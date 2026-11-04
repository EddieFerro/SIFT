#include <iostream>
#include <fstream>
#include <algorithm>
#include <filesystem>
#include <chrono>
#include <atomic>
#include <limits>
#include <malloc.h>
#include <sdsl/qsufsort.hpp>

#include "mem_finder.hpp"
#include "pfsp_ext.hpp"
#include "read_parser.hpp"
#include "doc_array.hpp"
#include "config.hpp"
#include "CLI11.hpp"
#include "zlib.h"
#include "kseq.h"
#include "omp.h"

KSEQ_INIT(gzFile, gzread)

// ── Timing ────────────────────────────────────────────────────────────────────
using Clock  = std::chrono::steady_clock;
using Micros = std::chrono::microseconds;
 
// Thread-safe per-step accumulator.
// Records total, count, min and max duration (in microseconds) across all reads
// processed by all threads.  Uses relaxed atomics — no ordering guarantees
// needed since each statistic is updated independently.
struct StepStats {
    std::atomic<long long> total_us{0};
    std::atomic<long long> count   {0};
    std::atomic<long long> min_us  {std::numeric_limits<long long>::max()};
    std::atomic<long long> max_us  {0};
 
    void record(long long us) {
        total_us.fetch_add(us, std::memory_order_relaxed);
        count.fetch_add(1,  std::memory_order_relaxed);
 
        // Atomic CAS loop for min
        long long old = min_us.load(std::memory_order_relaxed);
        while (us < old &&
               !min_us.compare_exchange_weak(old, us, std::memory_order_relaxed)) {}
 
        // Atomic CAS loop for max
        old = max_us.load(std::memory_order_relaxed);
        while (us > old &&
               !max_us.compare_exchange_weak(old, us, std::memory_order_relaxed)) {}
    }
 
    void print(const std::string& name) const {
        long long n = count.load();
        if (n == 0) return;
        double avg = static_cast<double>(total_us.load()) / n;
        std::cerr << "[timing] " << name
                  << "  total=" << total_us.load() / 1000 << "ms"
                  << "  avg="   << avg                    << "µs"
                  << "  min="   << min_us.load()          << "µs"
                  << "  max="   << max_us.load()          << "µs"
                  << "  n="     << n                      << "\n";
    }
};
// ─────────────────────────────────────────────────────────────────────────────


struct IndexParams {
    std::string pfsp_prefix;
    std::string temp_dir = DEFAULT_TEMP_PATH;
    std::string index_prefix;
    std::string tag_by = "pansn-hap";
    uint64_t tag_sample = 1;
    bool stats = false;
};

struct QueryParams {
    std::string reads;
    std::string pfsp_prefix;
    std::string index_prefix;
    std::string output;
    uint32_t top_t = DEFAULT_TOP_T;
    uint32_t threads = DEFAULT_THREADS;
    uint32_t L = DEFAULT_MEM_LENGTH;
    bool tags = true;
};


static sift::TagBy parse_tag_by(const std::string& s) {
    if (s == "record") return sift::TagBy::Record;
    if (s == "pansn-hap") return sift::TagBy::PanSNHap;
    if (s == "pansn-sample") return sift::TagBy::PanSNSample;
    throw std::runtime_error("unknown --tag-by " + s);
}


void process_reads(QueryParams& params, sift::index_t& fwd_index, sift::index_t& rev_index, sift::ReadParser& RP, const sift::DocArray* tags) {

    gzFile in;
    in = gzopen(params.reads.c_str(), "r");
    if(in == 0) {
        throw std::runtime_error("Error opening file: " + params.reads);
    }
    kseq_t* sequences = kseq_init(in);

    FILE* out = fopen(params.output.c_str(), "w");
    if (!out) {
        kseq_destroy(sequences); gzclose(in);
        throw std::runtime_error("cannot write " + params.output);
    }
    int fd = fileno(out);
    int64_t buffer = fpathconf(fd, _PC_REC_XFER_ALIGN);
    if (buffer <= 0) buffer = 4096;
    setvbuf(out, nullptr, _IOFBF, buffer);

    omp_set_num_threads(params.threads);

    std::atomic<size_t> total_fragments{0};
    int read_error = 0;   // kseq_read < -1: truncated quality string (-2) or read error (-3)

    // Per-step timing accumulators
    StepStats t_parse, t_pseudo, t_tags, t_write, t_io;
    std::atomic<uint64_t> total_psi_steps{0}, total_tagged{0};
    auto psi = [&fwd_index](uint64_t i) { return fwd_index.psi(i); };


    #pragma omp parallel 
    {
        thread_local static std::vector<sift::PhraseToken> read_parse;
        thread_local static std::vector<std::tuple<int, int, sift::index_t::range_type>> pseudo_mems;

        char* seq_text = nullptr;
        char* seq_name = nullptr;
        int64_t seq_len = 0;

        while(true) {

            auto t_io_start = Clock::now();
            #pragma omp critical(read_seq) 
            {
                seq_len = kseq_read(sequences);
                if(seq_len < -1) read_error = (int) seq_len;
                if(seq_len >= 0) {
                    if(params.threads > 1) {
                        seq_text = strdup(sequences->seq.s);
                        seq_name = strdup(sequences->name.s);
                    } else {
                        seq_text = sequences->seq.s;
                        seq_name = sequences->name.s;
                    }
                }
            }
            if(seq_len < 0) break;
            // empty record: nothing to query, keep going
            if(seq_len == 0) { 
                if(params.threads > 1) { free(seq_text); free(seq_name); }
                continue;
            }

            t_io.record(std::chrono::duration_cast<Micros>(Clock::now() - t_io_start).count());

            read_parse.clear();
            pseudo_mems.clear();

            {
                auto t0 = Clock::now();
                RP.parse(seq_text, (size_t) seq_len, read_parse);
                t_parse.record(std::chrono::duration_cast<Micros>(Clock::now() - t0).count());
            }

            {
                auto t0 = Clock::now();
                sift::find_pseudo_mems(read_parse, params.L, fwd_index, rev_index, pseudo_mems, params.top_t);
                t_pseudo.record(std::chrono::duration_cast<Micros>(Clock::now() - t0).count());
            }


            // RP.parse(std::string(seq_text, seq_len), read_parse);
            // sift::find_pseudo_mems(read_parse, candidate_mems, params.L, fwd_index, rev_index, pseudo_mems, params.top_t);

            // rc_parse.clear();
            // rc_pseudo_mems.clear();

            // RP.parse(std::string(sift::reverse_complement(seq_text, seq_len)), rc_parse);
            // sift::find_pseudo_mems(rc_parse, rc_candidate_mems, params.L, fwd_index, rev_index, rc_pseudo_mems);

            
            thread_local static std::vector<std::pair<int,int>> all_fragments;
            thread_local static std::vector<std::string> all_tags;
            thread_local static std::vector<uint64_t> tag_ids;
            all_fragments.clear();
            all_tags.clear();
            const auto t_tags0 = Clock::now();
            for(auto& [ml, mr, sa]: pseudo_mems) {
                all_fragments.push_back({read_parse[ml].start, read_parse[mr].end});
                std::string t;
                if (tags) {
                    tag_ids.clear();
                    tags->tags(sa.start, sa.end, tag_ids, psi);
                    total_tagged.fetch_add(1, std::memory_order_relaxed);
                    std::sort(tag_ids.begin(), tag_ids.end());
                    for (size_t k = 0; k < tag_ids.size(); ++k) {
                        if (k) t += ',';
                        t += tags->tag_name(tag_ids[k]);
                    }
                }
                all_tags.push_back(std::move(t));
            }
            if (tags) t_tags.record(std::chrono::duration_cast<Micros>(Clock::now() - t_tags0).count());

            // for(auto& [ml, mr, sa]: rc_pseudo_mems) {
            //     int fwd_start = seq_len - rc_parse[mr].end;
            //     int fwd_end = seq_len - rc_parse[ml].start;

            //     all_fragments.push_back({fwd_start, fwd_end});
            // }

            // if(params.top_t > 0 && all_fragments.size() > params.top_t) {
            //     std::sort(all_fragments.begin(), all_fragments.end(),
            //         [](const auto& a, const auto& b) {
            //             return (a.second - a.first) > (b.second - b.first);
            //         });
            //     all_fragments.resize(params.top_t);
            // }
            
            if(all_fragments.empty())
            {
                all_fragments.push_back({0, seq_len});
                all_tags.emplace_back();
            }
            thread_local static std::string out_buf;
            out_buf.clear();

            for (size_t f = 0; f < all_fragments.size(); ++f) {
                const auto [start, end] = all_fragments[f];
                // Write header
                out_buf += '>';
                out_buf += seq_name;
                out_buf += ':';
                out_buf += std::to_string(start + 1);
                out_buf += '-';
                out_buf += std::to_string(end);
                if (!all_tags[f].empty()) {
                    out_buf += " tags=";
                    out_buf += all_tags[f];
                }
                out_buf += '\n';
                // Append sequence
                out_buf.append(seq_text + start, end - start);
                out_buf += '\n';
            }

            {
                
                #pragma omp critical(write_mems)
                {
                    auto t0 = Clock::now();
                    total_fragments += all_fragments.size();
                    fwrite(out_buf.data(), 1, out_buf.size(), out);
                    t_write.record(std::chrono::duration_cast<Micros>(Clock::now() - t0).count());
                }
                
            }
            
            if(params.threads > 1) {
                free(seq_text);
                free(seq_name);
            }
            total_psi_steps.fetch_add(sift::DocArray::psi_steps, std::memory_order_relaxed);
            sift::DocArray::psi_steps = 0;
        }
    }

    // Summary 
    std::cerr << "\n[sift] Total fragments generated: " << total_fragments << "\n";
    t_parse.print("parse          ");
    t_pseudo.print("find_pseudo    ");
    if (tags) {
        t_tags.print("tags           ");
        std::cerr << "[sift] tags: s=" << tags->sample_step() << ", " << total_tagged << " tagged fragments, "
                  << total_psi_steps << " Psi steps ("
                  << (total_tagged ? (double)total_psi_steps / total_tagged : 0.0) << " per fragment)\n";
    }
    t_write.print("write_output   ");
    t_io.print("read_seq       ");

    kseq_destroy(sequences);
    gzclose(in);
    const bool write_failed = std::ferror(out) != 0;
    if (std::fclose(out) != 0 || write_failed) throw std::runtime_error("error writing " + params.output);
    if (read_error == -2) throw std::runtime_error(params.reads + ": truncated FASTQ quality string");
    if (read_error) throw std::runtime_error("error reading " + params.reads);


}


static void sa_bwt(const std::string& path, uint64_t len, uint32_t sigma,
                   sdsl::int_vector<>& SA, std::vector<uint32_t>& bwt) {
    std::vector<uint32_t> P(len);
    std::ifstream in(path, std::ios::binary);
    if (!in.read(reinterpret_cast<char*>(P.data()), len * 4)) throw std::runtime_error("cannot read " + path);
    const uint64_t n = len + 1;
    { sdsl::int_vector<> T(n, 0, sdsl::bits::hi(sigma) + 1);
      for (uint64_t i = 0; i < len; ++i) T[i] = P[i];
      sdsl::qsufsort::construct_sa(SA, T); }
    bwt.resize(n);
    for (uint64_t i = 0; i < n; ++i) bwt[i] = SA[i] ? P[SA[i] - 1] : 0;
}


static int run(int argc, char** argv) {
    //Parameter Objects
    IndexParams index_params;
    QueryParams query_params;
    
    
    //Argument Parsing
    CLI::App app{"SIFT: pangenome index on a PFSP syncmer parse; pseudo-MEM queries with haplotype tags"};
    app.require_subcommand(1);
    app.failure_message(CLI::FailureMessage::help);

    //Building Index subcommand and options
    auto build = app.add_subcommand("build", "Build the SIFT RLCSA indexes.");

    build->add_option("pfsp_prefix", index_params.pfsp_prefix, "PFSP output prefix (<prefix>.meta/.dict/.parse/.gaps).")->required();
    build->add_option("-o,--output", index_params.index_prefix, "Prefix for index files. [PREFIX].[fwd/rev].rlcsa");
    build->add_option("-t,--temp", index_params.temp_dir, "Unused: construction runs in memory (kept so existing scripts work).");
    build->add_option("--tag-by", index_params.tag_by, "Tag granularity: record | pansn-hap (sample#hap) | pansn-sample")
        ->check(CLI::IsMember({"record", "pansn-hap", "pansn-sample"}))->capture_default_str();
    build->add_option("--tag-sample", index_params.tag_sample, "Store the tag of every s-th syncmer of each record and recover the rest with at most s-1 Psi steps (1 = store every tag).")
        ->check(CLI::PositiveNumber)->capture_default_str();
    build->add_flag("--stats", index_params.stats, "Also report the number of BWT runs r of the parse (loads the parse into memory).");

    auto query = app.add_subcommand("query", "Uses the SIFT index to identify pseudo-mems.");

    query->add_option("input", query_params.reads, "Input reads file for SIFT querying. Only Fasta/Fastq files supported currently.")->check(CLI::ExistingFile)->required();
    query->add_option("-d,--pfsp", query_params.pfsp_prefix, "PFSP output prefix used to build the index (<prefix>.meta/.dict).")->required();
    query->add_option("-i,--index", query_params.index_prefix, "Prefix for index files. [PREFIX].[fwd/rev].rlcsa")->required();
    query->add_option("-o,--output", query_params.output, "Output file (fasta format).");
    query->add_flag("!--no-tags", query_params.tags, "Do not load the tag array or report tags.");
    query->add_option("--top-t", query_params.top_t, "Only top-t longest fragments are kept.")->check(CLI::PositiveNumber);
    query->add_option("-t,--threads", query_params.threads, "Number of threads.")->check(CLI::PositiveNumber);
    query->add_option("-L,--length", query_params.L, "Minimum MEM length.")->check(CLI::PositiveNumber);


    CLI11_PARSE(app, argc, argv);

    if(build->parsed()) {
        if(index_params.index_prefix.empty()) index_params.index_prefix = index_params.pfsp_prefix;
        const std::string& pf = index_params.pfsp_prefix;
        const std::string& ix = index_params.index_prefix;

        // Index input: the PFSP syncmer parse itself, ranks shifted to 1-based
        sift::Meta meta = sift::Meta(pf + ".meta");
        const std::string parse_path = pf + ".parse";
        const std::string rev_path = ix + ".parse.rev";
        if (sift::construct_reverse_parse(parse_path, rev_path) != meta.n_syncmers)
            throw std::runtime_error(parse_path + ": length does not match n_syncmers in .meta");
        std::cerr << "[sift] " << meta.n_syncmers << " syncmer occurrences over "
                  << meta.n_distinct << " distinct syncmers, " << meta.recs.size() << " records\n";


        // Forward index, then the tag array from its SA (still in the cache)
        const uint32_t sigma = (uint32_t)meta.n_distinct + 1;
        {
            sdsl::int_vector<> SA;
            std::vector<uint32_t> bwt;
            sa_bwt(parse_path, meta.n_syncmers, sigma, SA, bwt);
            sift::index_t fwd_index;
            fwd_index.build(bwt, sigma);
            bwt = std::vector<uint32_t>(); // free before the tag array
            fwd_index.save(ix + ".fwd" + SIFT_EXT);
            std::cerr << "[sift] fwd: n=" << fwd_index.size() << " r=" << fwd_index.runs()
                      << " n/r=" << (double)fwd_index.size() / fwd_index.runs()
                      << ", " << fwd_index.bytes() << " bytes\n";

            sift::DocArray tags;
            tags.build(SA, meta.record_starts(), meta.record_names(), parse_tag_by(index_params.tag_by), index_params.tag_sample);

            if (tags.size() != fwd_index.size())
                throw std::runtime_error("tag array size does not match the forward index");
            tags.save(ix + SIFT_TAGS_EXT);
            std::cerr << "[sift] tags: " << tags.n_tags() << " tags, n=" << tags.size()
                      << ", tag runs t=" << tags.runs() << ", s=" << tags.sample_step()
                      << ", " << tags.n_samples() << " stored tags, " << tags.bytes() << " bytes ("
                      << 8.0 * tags.bytes() / tags.size() << " bits/row)\n";

            // Psi must map row i to the row of SA[i] + 1: sampled tags depend on it
            {
                const uint64_t n = SA.size(), step = std::max<uint64_t>(1, n / 1000);
                for (uint64_t i = 0; i < n; i += step) {
                    const uint64_t want = (SA[i] + 1) % n, got = fwd_index.psi(i);
                    if (got >= n || SA[got] != want)
                        throw std::runtime_error("index_t::psi(" + std::to_string(i) + ") is not the row of SA[i] + 1");
                }
            }
            
            if (index_params.stats)
                std::cerr << "[sift] stats: n=" << fwd_index.size() << " r=" << fwd_index.runs()
                          << " n/r=" << (double)fwd_index.size() / fwd_index.runs() << " t=" << tags.runs()
                          << " t/n=" << (double)tags.runs() / fwd_index.size() << "\n";
        }
        malloc_trim(0);


        // Reverse index (no tags)
        {
            sdsl::int_vector<> SA;
            std::vector<uint32_t> bwt;
            sa_bwt(rev_path, meta.n_syncmers, sigma, SA, bwt);
            SA = sdsl::int_vector<>();
            sift::index_t rev_index;
            rev_index.build(bwt, sigma);
            rev_index.save(ix + ".rev" + SIFT_EXT);
            std::cerr << "[sift] rev: r=" << rev_index.runs() << ", " << rev_index.bytes() << " bytes\n";
        }

        std::filesystem::remove(rev_path);
    }

    if(query->parsed()) {
        StepStats t_load,t_dict;

        // Check existence of index files
        std::filesystem::path f = query_params.index_prefix + ".fwd" + SIFT_EXT;
        std::filesystem::path r = query_params.index_prefix + ".rev" + SIFT_EXT;
        // std::filesystem::path d = query_params.dictionary;

        if(!(std::filesystem::exists(f) && std::filesystem::exists(r))) throw std::runtime_error("Index files not found.");
        // if(d.extension().string() != vcfbwt::EXT::DICT) throw std::runtime_error("Incorrect type for dictionary file.");

        if(query_params.output.empty()) {
            std::filesystem::path o = query_params.reads;
            query_params.output = o.replace_extension("").string() + SIFT_MEMS_EXT;
        }

        auto t0 = Clock::now();
        // Load indexes
        sift::index_t fwd_index;
        // sdsl::load_from_file(fwd_index, f.string());
        fwd_index.load(f.string());

        sift::index_t rev_index;
        // sdsl::load_from_file(rev_index, r.string());
        rev_index.load(r.string());
        
        std::unique_ptr<sift::DocArray> tags;
        if (query_params.tags) {
            tags = std::make_unique<sift::DocArray>();
            tags->load(query_params.index_prefix + SIFT_TAGS_EXT);
            if (tags->size() != fwd_index.size())
                throw std::runtime_error("tag array does not match the forward index");
        }
        t_load.record(std::chrono::duration_cast<Micros>(Clock::now() - t0).count());
        t_load.print("load indexes   ");

        t0 = Clock::now();
        //Create readparser, which maintains dictionary
        sift::ReadParser RP(query_params.pfsp_prefix);
        t_dict.record(std::chrono::duration_cast<Micros>(Clock::now() - t0).count());
        t_dict.print("load dict   ");
        

        process_reads(query_params, fwd_index, rev_index, RP, tags.get());

        // D.load_dictionary(query_params.dictionary);
        // test_integration_mems(D, query_params.pfp_params, fwd_index, rev_index, query_params.reads, query_params.L);
    }

    return 0;
}

int main(int argc, char** argv) {
    try {
        return run(argc, argv);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "[sift] error: %s\n", e.what());
        return 1;
    }
}
