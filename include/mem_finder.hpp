#pragma once
#include <algorithm>
#include <extended_rlcsa.hpp>
#include <read_parser.hpp>
#include "config.hpp"
#include <queue>


namespace sift {

using index_t = EFRLCSA;


int first_long_substring(const std::vector<PhraseToken>& tokens, int s, int L) {
    const uint32_t target = tokens[s].start + (uint32_t)L;
    auto it = std::lower_bound(tokens.begin() + s, tokens.end(), target, [](const PhraseToken& t, uint32_t v) { return t.end < v; });
    if(it == tokens.end()) return -1;
    return int(it - tokens.begin());
}


void find_pseudo_mems(const std::vector<PhraseToken>& read_parse,
    int L,
    const index_t& index_fwd, 
    const index_t& index_rev,
    std::vector<std::tuple<int, int, typename index_t::range_type>>& mems,
    int top_t = 0) 
{

    using range_t = typename index_t::range_type;
    using mem_t = std::tuple<int, int, range_t>;

    using heap_entry = std::pair<int, mem_t>;
    
    auto cmp = [](const heap_entry& a, const heap_entry& b) {
        return a.first > b.first;
    };

    std::priority_queue<heap_entry, std::vector<heap_entry>, decltype(cmp)> min_heap(cmp);

    // Helper: try to add a MEM to the heap
    auto try_add = [&](int out_l, int out_r, range_t range) {
        int span = read_parse[out_r].end - read_parse[out_l].start;
        mem_t mem{out_l, out_r, range};
        if(top_t == 0) {
            mems.emplace_back(mem);
        } else if((int)min_heap.size() < top_t) {
            min_heap.push({span, mem});
        } else if(span > min_heap.top().first) {
            min_heap.pop();
            min_heap.push({span, mem});
        }
    };

    const int read_len = (int) read_parse.size();

    int l = 0;
    while(l < read_len) {

        if(read_parse[l].id == 0) {
            l++;
            continue;
        }
        int r = l;
        while(r+1 < read_len && read_parse[r + 1].id != 0) r++;
        const int next_l = r + 1;

        const int run_span = read_parse[r].end - read_parse[l].start;
        if(run_span < L || (top_t > 0 && (int)min_heap.size() >= top_t && run_span <= min_heap.top().first)) {
            l = next_l;
            continue;
        }

        int s = l;  // current search start (absolute index into read_parse)

        while(s <= r) {
            // Gagie L-skip
            // Backward search for P[s..q] in index_fwd (right-to-left), where
            // [s..q] spans >= L characters.  If it fails at phrase k, P[k..q]
            // is absent so no MEM of length >= L starts in [s..k]; advance to k+1.
            int q = first_long_substring(read_parse, s, L);
            if(q < 0 || q > r) break;

            range_t interval = index_fwd.full_range();
            bool skipped = false;
            for(int k = q; k >= s; --k) {
                auto next = index_fwd.extend_left(interval, read_parse[k].id);
                if(index_fwd.empty(next)) { s = k + 1; skipped = true; break; }
                interval = next;
            }
            if(skipped) continue;

            // Forward extension
            // "search for P^rev in the index for T^rev": search P[s..] left-to-right
            // in index_rev to find the right end e of the MEM starting at s.
            range_t range_rev = index_rev.full_range();
            int e = s - 1;
            for(int j = s; j <= r; ++j) {
                auto next = index_rev.extend_left(range_rev, read_parse[j].id);
                if(index_rev.empty(next)) break;
                range_rev = next; e = j;
            }
            if(e < s) { s++; continue; }

            // Get the forward SA interval for MEM P[s..e] via backward search in index_fwd.
            range_t range_fwd = index_fwd.full_range();
            for(int k = e; k >= s; --k)
                range_fwd = index_fwd.extend_left(range_fwd, read_parse[k].id);

            // Left-maximality is guaranteed by BML structure.
            // Apply boundary zero folding for candidates with adjacent zeros.
            int out_l = (s > 0) ? s - 1 : s;
            int out_r = (e < read_len - 1) ? e + 1 : e;
            // mems.emplace_back(out_l, out_r, range_fwd);
            try_add(out_l, out_r, range_fwd);

            if(e == r) break;

            // Backward extension to find next start s2
            // "find s2 through backward extension by searching for P[0..e1+1] in
            // the index for T": search P[s..e+1] right-to-left in index_fwd.
            // Failure at phrase k means P[k..e+1] is absent; next start = k+1.
            range_t rb = index_fwd.full_range();
            int s2 = e + 1;
            for(int k = e + 1; k >= s; --k) {
                auto next = index_fwd.extend_left(rb, read_parse[k].id);
                if(index_fwd.empty(next)) { s2 = k + 1; break; }
                rb = next;
            }
            if(s2 <= s) s2 = s + 1;
            s = s2;
        }
        l = next_l;
    }
    while(!min_heap.empty()) {
        mems.push_back(std::get<1>(min_heap.top()));
        min_heap.pop();
    }
}


} //namespace sift
