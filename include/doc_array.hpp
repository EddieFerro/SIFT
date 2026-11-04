#pragma once

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/util.hpp>

namespace sift {

// DA[i] = tag of the suffix at SA position i of the forward syncmer parse, where a suffix's tag is the tag of the record its first syncmer came from.
// DA[Psi(i)] = DA[i] unless SA[i] is the last syncmer of its record
// tags(b, e) lists the distinct tags in SA interval [b, e) by scanning it, stopping once every tag has been seen; count(b, e) gives the number of rows with each tag.

// How records map to tags.
enum class TagBy { Record, PanSNHap, PanSNSample };

// record name -> tag string. PanSN names look like sample#hap#contig.
inline std::string tag_of_name(const std::string &name, TagBy by) {
    if (by == TagBy::Record) return name;
    const auto h1 = name.find('#');
    if (h1 == std::string::npos) return name; // not PanSN: fall back to the record
    if (by == TagBy::PanSNSample) return name.substr(0, h1);
    const auto h2 = name.find('#', h1 + 1);
    return h2 == std::string::npos ? name : name.substr(0, h2);
}

class DocArray {
    static constexpr char MAGIC[9] = "PFPITAG2";

    uint64_t n_ = 0; // SA length (parse length + 1 sentinel)
    uint64_t n_tags_ = 0;
    uint64_t s_ = 1; // sampling step along Psi (1 = every row stored)
    uint64_t runs_ = 0; // tag runs over all n rows (statistic, computed at build)
    sdsl::sd_vector<> sampled_;  // s > 1: sampled rows
    sdsl::sd_vector<>::rank_1_type sampled_rank_;
    sdsl::int_vector<> tag_; // s = 1: DA; s > 1: tags of sampled rows
    std::vector<std::string> tag_names_;

    // distinct-tag bookkeeping without clearing an n_tags-sized array per query
    struct Seen {
        std::vector<uint32_t> stamp;
        uint32_t gen = 0;
        void reset(size_t n) {
            if (stamp.size() < n) stamp.assign(n, 0), gen = 0;
            if (++gen == 0) { std::fill(stamp.begin(), stamp.end(), 0); gen = 1; }
        }
        bool insert(uint64_t t) { if (stamp[t] == gen) return false; stamp[t] = gen; return true; }
    };

public:
    // Psi steps taken by tag_at() on this thread
    static inline thread_local uint64_t psi_steps = 0;

    DocArray() = default;
    DocArray(const DocArray &) = delete; // sampled_rank_ points into sampled_
    DocArray &operator=(const DocArray &) = delete;

    uint64_t size() const { return n_; }
    uint64_t n_tags() const { return n_tags_; }
    uint64_t sample_step() const { return s_; }
    uint64_t n_samples() const { return tag_.size(); }
    uint64_t runs() const { return runs_; }
    const std::string &tag_name(uint64_t t) const { return tag_names_[t]; }
    uint64_t bytes() const {
        return sdsl::size_in_bytes(tag_) + (s_ > 1 ? sdsl::size_in_bytes(sampled_) + sdsl::size_in_bytes(sampled_rank_) : 0);
    }

    // sa: the forward SA indexable and with .size(); read once, in order
    // doc_first: first parse entry of each record, ascending, plus a final entry = parse length (PfspMeta::record_starts())
    // doc_names: record names, same order
    // s: sampling step (>= 1)
    template <class SA>
    void build(SA &&sa, const std::vector<uint64_t> &doc_first,
               const std::vector<std::string> &doc_names, TagBy by, uint64_t s = 1) {
        if (s == 0) throw std::runtime_error("DocArray: sampling step must be >= 1");
        const uint64_t n_par = doc_first.back();   // parse length
        n_ = sa.size();
        s_ = s;
        if (n_ != n_par + 1)
            throw std::runtime_error("DocArray: SA length " + std::to_string(n_) +
                                     " != parse length + 1 (" + std::to_string(n_par + 1) + ")");
        if (doc_names.size() + 1 != doc_first.size())
            throw std::runtime_error("DocArray: record starts and names disagree");
        // record -> tag id, in order of first appearance
        tag_names_.clear();
        std::vector<uint64_t> rec_tag(doc_names.size());
        {
            std::vector<std::pair<std::string, uint64_t>> seen;
            for (size_t r = 0; r < doc_names.size(); ++r) {
                const std::string t = tag_of_name(doc_names[r], by);
                auto it = std::find_if(seen.begin(), seen.end(), [&](auto &p) { return p.first == t; });
                if (it == seen.end()) { seen.push_back({t, tag_names_.size()}); tag_names_.push_back(t); rec_tag[r] = seen.back().second; }
                else rec_tag[r] = it->second;
            }
        }
        n_tags_ = tag_names_.size();
        const uint8_t width = sdsl::bits::hi(n_tags_) + 1; // values 0..n_tags_

        // records with no parse entries have empty ranges; upper_bound skips them
        auto rec_of = [&](uint64_t v) -> size_t {
            return std::upper_bound(doc_first.begin(), doc_first.end() - 1, v) - doc_first.begin() - 1;
        };

        sampled_ = sdsl::sd_vector<>();
        sampled_rank_ = sdsl::sd_vector<>::rank_1_type();
        runs_ = 0;
        uint64_t last = UINT64_MAX;
        if (s_ == 1) {
            tag_ = sdsl::int_vector<>(n_, 0, width);
            for (uint64_t i = 0; i < n_; ++i) {
                const uint64_t v = sa[i];
                const uint64_t t = v >= n_par ? n_tags_ : rec_tag[rec_of(v)];
                tag_[i] = t;
                if (t != last) { ++runs_; last = t; }
            }
            return;
        }
        uint64_t m = 1; // sentinel row
        for (size_t r = 0; r + 1 < doc_first.size(); ++r) m += (doc_first[r + 1] - doc_first[r] + s_ - 1) / s_;
        sdsl::sd_vector_builder sb(n_, m);
        tag_ = sdsl::int_vector<>(m, 0, width);
        uint64_t k = 0;
        for (uint64_t i = 0; i < n_; ++i) {
            const uint64_t v = sa[i];
            uint64_t t = n_tags_;
            bool keep = true; // the sentinel row is always sampled
            if (v < n_par) {
                const size_t r = rec_of(v);
                t = rec_tag[r];
                keep = (doc_first[r + 1] - 1 - v) % s_ == 0;
            }
            if (t != last) { ++runs_; last = t; }
            if (!keep) continue;
            if (k == m) throw std::runtime_error("DocArray: more sampled rows than expected (SA is not a permutation?)");
            sb.set(i);
            tag_[k++] = t;
        }
        if (k != m) throw std::runtime_error("DocArray: fewer sampled rows than expected (SA is not a permutation?)");
        sampled_ = sdsl::sd_vector<>(sb);
        sampled_rank_ = sdsl::sd_vector<>::rank_1_type(&sampled_);
    }

    // DA[i]. psi(i) must return the row of SA[i] + 1 (index_t::psi); it is not called when s == 1.
    template <class Psi> uint64_t tag_at(uint64_t i, const Psi &psi) const {
        if (s_ == 1) return tag_[i];
        uint64_t steps = 0;
        while (!sampled_[i]) {
            i = psi(i);
            if (++steps >= s_) throw std::logic_error("DocArray: Psi walk exceeded the sampling step (wrong Psi or index?)");
        }
        psi_steps += steps;
        return tag_[sampled_rank_(i)];
    }

    // Distinct tags occurring in SA interval [b, e). Appends to `out` (unsorted).
    template <class Psi> void tags(uint64_t b, uint64_t e, std::vector<uint64_t> &out, const Psi &psi) const {
        if (b >= e || e > n_) return;
        thread_local Seen seen;
        seen.reset(n_tags_ + 1);
        seen.insert(n_tags_); // never report the sentinel
        uint64_t found = 0;
        for (uint64_t i = b; i < e && found < n_tags_; ++i) {
            const uint64_t t = tag_at(i, psi);
            if (seen.insert(t)) { ++found; out.push_back(t); }
        }
    }

    // Rows of [b, e) per tag, as (tag, count) pairs sorted by tag.
    template <class Psi> void count(uint64_t b, uint64_t e, std::vector<std::pair<uint64_t, uint64_t>> &out, const Psi &psi) const {
        out.clear();
        if (b >= e || e > n_) return;
        thread_local std::vector<uint64_t> cnt;
        thread_local std::vector<uint64_t> touched;
        if (cnt.size() < n_tags_ + 1) cnt.assign(n_tags_ + 1, 0);
        touched.clear();
        for (uint64_t i = b; i < e; ++i) {
            const uint64_t t = tag_at(i, psi);
            if (t == n_tags_) continue;
            if (cnt[t]++ == 0) touched.push_back(t);
        }
        for (uint64_t t : touched) { out.push_back({t, cnt[t]}); cnt[t] = 0; }
        std::sort(out.begin(), out.end());
    }

    void save(const std::string &path) const {
        std::ofstream file(path, std::ios::binary);
        std::ostream &f = file;
        f.write(MAGIC, 8);
        sdsl::write_member(n_, f);
        sdsl::write_member(n_tags_, f);
        sdsl::write_member(s_, f);
        sdsl::write_member(runs_, f);
        tag_.serialize(f);
        if (s_ > 1) sampled_.serialize(f);
        const uint64_t nt = tag_names_.size();
        sdsl::write_member(nt, f);
        for (const auto &s : tag_names_) {
            const uint64_t l = s.size();
            sdsl::write_member(l, f);
            f.write(s.data(), (std::streamsize)l);
        }
        if (!f) throw std::runtime_error("error writing " + path);
    }
    void load(const std::string &path) {
        std::ifstream file(path, std::ios::binary);
        if (!file) throw std::runtime_error("cannot open " + path);
        std::istream &f = file;
        char magic[8];
        if (!f.read(magic, 8) || std::string(magic, 8) != std::string(MAGIC, 8))
            throw std::runtime_error(path + ": not a " + std::string(MAGIC, 8) +
                                     " tag array (rebuild it with this version of `build`)");
        sdsl::read_member(n_, f);
        sdsl::read_member(n_tags_, f);
        sdsl::read_member(s_, f);
        sdsl::read_member(runs_, f);
        tag_.load(f);
        sampled_ = sdsl::sd_vector<>();
        sampled_rank_ = sdsl::sd_vector<>::rank_1_type();
        if (s_ > 1) {
            sampled_.load(f);
            sampled_rank_ = sdsl::sd_vector<>::rank_1_type(&sampled_);
        }
        uint64_t nt = 0;
        sdsl::read_member(nt, f);
        tag_names_.resize(nt);
        for (auto &s : tag_names_) {
            uint64_t l = 0;
            sdsl::read_member(l, f);
            s.resize(l);
            f.read(s.data(), (std::streamsize)l);
        }
        if (!f) throw std::runtime_error(path + ": truncated tag array");
    }
};

} // namespace sift
