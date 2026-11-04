#pragma once

#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>
#include "efrank.hpp"

namespace sift {

class EFRLCSA {
public:
    struct range_type { uint64_t start, end; }; // half-open [start, end)

    // bwt[i] = T[SA[i]-1] (0 for the sentinel row); symbols 0..sigma-1
    void build(const std::vector<uint32_t>& bwt, uint32_t sigma) {
        n_ = bwt.size(); ef_.build(bwt, sigma);
    }

    uint64_t size() const { return n_; }
    range_type full_range() const { return {0, n_}; }
    bool empty(const range_type& r) const { return !(r.start < r.end); }

    // Backward search step: interval of c·P from the interval of P.
    range_type extend_left(const range_type& cur, uint64_t c) const {
        if (c >= ef_.sigma) return {0, 0};
        const uint64_t base = c * n_;
        const uint64_t b = ef_.KS.count_less(base); // first run of c
        const uint64_t e = ef_.KS.count_less(base + n_); // one past its last run
        if (b == e) return {0, 0}; // c does not occur
        const uint64_t Cc = ef_.FP.access(b); // C[c]
        if (cur.start == 0 && cur.end == n_) { // first step: block of F
            const uint64_t Cn = e < ef_.nruns ? ef_.FP.access(e) : n_;
            return {Cc, Cn};
        }
        return {Cc + rank_(b, base, c, cur.start), Cc + rank_(b, base, c, cur.end)};
    }

    // psi(i): row of SA[i] + 1; select on the BWT. Row i is in F-run j.
    uint64_t psi(uint64_t i) const {
        const uint64_t j = ef_.FP.count_less(i + 1) - 1;
        const uint64_t v = ef_.KS.access(j); // c * n + start_j
        return v % n_ + (i - ef_.FP.access(j));
    }

    uint64_t runs() const { return ef_.nruns; } // r of the BWT 
    uint64_t bytes() const { return ef_.bytes(); }
    void save(const std::string& path) const {
        FILE* f = std::fopen(path.c_str(), "wb"); if (!f) throw std::runtime_error("cannot write " + path);
        std::fwrite(MAGIC, 1, 8, f); std::fwrite(&n_, 8, 1, f); ef_.save(f); std::fclose(f);
    }
    void load(const std::string& path) {
        FILE* f = std::fopen(path.c_str(), "rb"); if (!f) throw std::runtime_error("cannot read " + path);
        char m[8]; if (std::fread(m, 1, 8, f) != 8 || std::string(m, 8) != std::string(MAGIC, 8))
            throw std::runtime_error(path + ": not a SIFT EF index");
        if (std::fread(&n_, 8, 1, f) != 1) throw std::runtime_error("read " + path);
        ef_.load(f); std::fclose(f);
    }

private:
    static constexpr char MAGIC[9] = "PFPIEF01";
    uint64_t n_ = 0;
    EFRank ef_;

    // occurrences of c in BWT[0, i), given b = index of c's first run
    uint64_t rank_(uint64_t b, uint64_t base, uint64_t c, uint64_t i) const {
        if (i == 0) return 0;
        const uint64_t a = ef_.KS.count_less(base + (i > n_ ? n_ : i));
        return a == b ? 0 : ef_.rank_from(b, a, (uint32_t)c, i);
    }
};

} // namespace sift
