#pragma once

#include <memory>
#include <string>
#include <vector>
#include "pfsp_ext.hpp"
#include "config.hpp"
#include "syncmer_parsing.hpp"

namespace sift {

struct PhraseToken {
    rank_t id; //phrase id
    uint32_t start; //beginning of phrase in read
    uint32_t end; //end of phrase in read

    PhraseToken(rank_t r, uint32_t s, uint32_t e) : id(r), start(s), end(e) {}
};

class ReadParser {
    SyncmerDict dict_;
    Meta meta_;

    

public:
    
    const Meta& meta() const { return meta_; }

    ReadParser(const std::string& prefix): meta_(prefix + ".meta") {
        dict_.load(prefix + ".dict", meta_.K, meta_.n_distinct);
    }

    void parse(const char *seq, size_t len, std::vector<PhraseToken>& read_parse) {
        read_parse.clear();
        const int K = meta_.K;

        if(len < (size_t)K) { 
            read_parse.emplace_back(0, 0, (uint32_t) len); 
            return;
        }

        thread_local std::unique_ptr<Syncmers> sm;
        sm = std::make_unique<Syncmers>(K, meta_.S, meta_.seed);
        thread_local std::vector<uint8_t> norm, key;

        norm.resize(len);
        for(size_t i = 0; i < len; i++) norm[i] = NORM.v[(uint8_t) seq[i]];
        key.resize(dict_.width());

        //Handle syncmers
        auto emit = [&](uint32_t p) {
            pack(norm.data() + p, key.data(), K);
            const uint32_t r = dict_.find(key.data());
            if (read_parse.empty() && p > 0) read_parse.emplace_back(0, 0, p);   // leading partial
            read_parse.emplace_back(r == SyncmerDict::NONE ? 0 : (rank_t)r + 1, p, p + (uint32_t)K);
        };

        sm->reset(norm.data());
        if(sm->is_syncmer()) emit(0);
        for(size_t i = (size_t) K; i < len; i++) {
            sm->shift(norm[i]);
            if(sm->is_syncmer()) emit((uint32_t)(i-K+1));
        }
        if (read_parse.empty()) { read_parse.emplace_back(0, 0, (uint32_t)len); return; }
        if (read_parse.back().end < len) read_parse.emplace_back(0, read_parse.back().end, (uint32_t)len);

    }
};



inline std::string reverse_complement(const char* seq, int len) {
    std::string rc(len, 'N');
    for(size_t i = 0; i < len; ++i) {
        switch(seq[len - 1 - i]) {
            case 'A': case 'a': rc[i] = 'T'; break;
            case 'T': case 't': rc[i] = 'A'; break;
            case 'C': case 'c': rc[i] = 'G'; break;
            case 'G': case 'g': rc[i] = 'C'; break;
            default: rc[i] = 'N'; break;
        }
    }
    return rc;
}

} //namespace sift
