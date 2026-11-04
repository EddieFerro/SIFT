#pragma once

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "utils.hpp"
#include "config.hpp"

namespace sift {

    struct Meta {
        int K = 0, S = 0;
        uint64_t seed = 0, n_syncmers = 0;
        uint32_t n_distinct = 0;
        std::vector<RecordInfo> recs;

        Meta(const std::string& path) {
            std::ifstream file(path);
            if(!file) throw std::runtime_error("cannot open " + path);
            bool fmt = false;

            std::string line;

            while(std::getline(file, line)) {
                if(line.empty()) continue;
                if(line[0] == '#') {
                    if(line.rfind("#name\t", 0) == 0) continue;
                    auto tab = line.find('\t');
                    const std::string k = line.substr(1, tab - 1), v = line.substr(tab + 1);
                    if (k == "format") fmt = (v == "pfsp-syncmer-v1");
                    else if (k == "K") K = std::stoi(v);
                    else if (k == "S") S = std::stoi(v);
                    else if (k == "seed") seed = std::stoull(v);
                    else if (k == "n_syncmers") n_syncmers = std::stoull(v);
                    else if (k == "n_distinct") n_distinct = (uint32_t)std::stoul(v);
                    else if ((k == "rank_base" && v != "1") || (k == "gap_bytes" && v != "2") ||(k == "id_bytes" && v != "4"))
                        throw std::runtime_error(path + ": unsupported " + k + " " + v);
                    continue;
                }
                std::istringstream ss(line);
                RecordInfo r;
                std::getline(ss, r.name, '\t');
                ss >> r.length >> r.first_entry >> r.n_entries;
                if (ss.fail()) throw std::runtime_error(path + ": bad record row: " + line);
                recs.push_back(std::move(r));
            }
            if (!fmt) throw std::runtime_error(path + ": not a pfsp-syncmer-v1 .meta file");
        }


        std::vector<uint64_t> record_starts() const {
            std::vector<uint64_t> v;
            v.reserve(recs.size() + 1);
            for (const auto &r : recs) v.push_back(r.first_entry);
            v.push_back(n_syncmers);
            return v;
        }

        std::vector<std::string> record_names() const {
            std::vector<std::string> v;
            for (const auto &r : recs) v.push_back(r.name);
            return v;
        }
    };

class SyncmerDict {
    int K_ = 0;
    size_t PB_ = 0;
    std::vector<uint8_t> recs_;        // n * PB_, rank order
    std::vector<uint32_t> slots_;      // open addressing: rank + 1, 0 = empty
    uint64_t mask_ = 0;
    static constexpr uint64_t SEED = 0x9e3779b97f4a7c15ULL;

    uint64_t hash(const uint8_t *r) const { return murmurhash64(r, PB_, SEED); }

public:
    static constexpr uint32_t NONE = UINT32_MAX;

    void load(const std::string &dict_path, int K, uint32_t n_expected) {
        K_ = K;
        PB_ = (3 * (size_t)K + 7) / 8;
        FILE *f = std::fopen(dict_path.c_str(), "rb");
        if (!f) throw std::runtime_error("cannot open " + dict_path);
        recs_.resize((size_t)n_expected * PB_);
        if (std::fread(recs_.data(), 1, recs_.size(), f) != recs_.size() || std::fgetc(f) != EOF) {
            std::fclose(f);
            throw std::runtime_error(dict_path + ": size does not match n_distinct * record width");
        }
        std::fclose(f);
        uint64_t cap = 16;
        while (cap < 2ull * n_expected) cap <<= 1;   // load factor <= 0.5
        slots_.assign(cap, 0);
        mask_ = cap - 1;
        // MurmurHash3 reads 8-byte blocks: hash an aligned copy of each record
        // (records sit at PB_-byte strides, so in place they are misaligned)
        std::vector<uint64_t> aligned((PB_ + 7) / 8);
        uint8_t *tmp = reinterpret_cast<uint8_t *>(aligned.data());
        for (uint32_t r = 0; r < n_expected; ++r) {
            std::memcpy(tmp, record(r), PB_);
            uint64_t i = hash(tmp) & mask_;
            while (slots_[i]) i = (i + 1) & mask_;
            slots_[i] = r + 1;
        }
    }

    int K() const { return K_; }
    uint32_t size() const { return (uint32_t)(recs_.size() / PB_); }
    size_t width() const { return PB_; }
    const uint8_t *record(uint32_t rank) const { return &recs_[(size_t)rank * PB_]; }

    // 0-based rank of a packed k-mer, or NONE (`packed` should be 8-byte aligned)
    uint32_t find(const uint8_t *packed) const {
        for (uint64_t i = hash(packed) & mask_; slots_[i]; i = (i + 1) & mask_)
            if (std::memcmp(record(slots_[i] - 1), packed, PB_) == 0) return slots_[i] - 1;
        return NONE;
    }
};

inline uint64_t construct_reverse_parse(const std::string& parse_path, const std::string& rev_path) {
    FILE *in = std::fopen(parse_path.c_str(), "rb");
    if (!in) throw std::runtime_error("cannot open " + parse_path);
    std::fseek(in, 0, SEEK_END);
    const uint64_t n = (uint64_t)std::ftell(in) / sizeof(rank_t);
    FILE *out = std::fopen(rev_path.c_str(), "wb");
    if (!out) { std::fclose(in); throw std::runtime_error("cannot open " + rev_path); }
    constexpr uint64_t CHUNK = 1u << 20;
    std::vector<rank_t> buf(CHUNK);
    bool bad = false;
    for (uint64_t left = n; left > 0 && !bad;) {
        const uint64_t m = std::min(CHUNK, left);
        left -= m;
        std::fseek(in, (long)(left * sizeof(rank_t)), SEEK_SET);
        if (std::fread(buf.data(), sizeof(rank_t), m, in) != m) { bad = true; break; }
        std::reverse(buf.begin(), buf.begin() + m);
        bad |= std::fwrite(buf.data(), sizeof(rank_t), m, out) != m;
    }
    std::fclose(in);
    bad |= std::fclose(out) != 0;
    if (bad) throw std::runtime_error("error writing " + rev_path);
    return n;
}
 

} //namespace sift
