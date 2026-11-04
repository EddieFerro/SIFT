#ifndef config_hpp
#define config_hpp
#include <cstddef>
#include <cstdint>


static constexpr const char* SIFT_EXT = ".rlcsa";
static constexpr const char* SIFT_TAGS_EXT = ".tags";
static constexpr const char* SIFT_MEMS_EXT = ".sift";
static constexpr const char* DEFAULT_TEMP_PATH = "./temp";
static constexpr const uint32_t DEFAULT_THREADS = 8;
static constexpr const uint32_t DEFAULT_MEM_LENGTH = 25;
static constexpr const uint32_t DEFAULT_TOP_T = 0;

using rank_t = uint32_t;   // PFSP phrase ranks are 4 bytes (.meta: id_bytes 4)

#endif
