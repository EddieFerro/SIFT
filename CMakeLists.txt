cmake_minimum_required(VERSION 3.20)
project(sift LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# PFSP checkout: its include/ provides utils.hpp, syncmer_parsing.hpp and
# dictionary.hpp, so syncmer selection and packing are byte-identical to the
# parser that produced the input. Default: the git submodule external/PFSP,
# else a PFSP checkout next to this one.
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/PFSP/include/syncmer_parsing.hpp")
  set(PFSP_DEFAULT "${CMAKE_CURRENT_SOURCE_DIR}/external/PFSP")
else()
  set(PFSP_DEFAULT "${CMAKE_CURRENT_SOURCE_DIR}/../PFSP")
endif()
set(PFSP_DIR "${PFSP_DEFAULT}" CACHE PATH "Path to the PFSP source tree")

if(NOT EXISTS "${PFSP_DIR}/include/syncmer_parsing.hpp")
  message(FATAL_ERROR "PFSP not found at ${PFSP_DIR}; pass -DPFSP_DIR=/path/to/PFSP")
endif()

include(FetchContent)

# -------------------------------------------------------------------------
# sdsl-lite (int_vector, sd_vector, qsufsort): the fork sr-index uses, pinned
# -------------------------------------------------------------------------
FetchContent_Declare(sdsl
  GIT_REPOSITORY https://github.com/duscob/sdsl-lite.git
  GIT_TAG        55a9b35f8a46f886a4846da23e03ff454ed29071   # branch feature/io-cache
)

set(SDSL_ENABLE_TESTS OFF CACHE BOOL "")
set(SDSL_ENABLE_TUTORIALS OFF CACHE BOOL "")
set(SDSL_ENABLE_EXAMPLES OFF CACHE BOOL "")
FetchContent_MakeAvailable(sdsl)


# -------------------------------------------------------------------------
# Single-file deps: CLI11 (command line), and kseq and MurmurHash3 as in PFSP's build
# -------------------------------------------------------------------------

FetchContent_Declare(cli11
  URL      https://github.com/CLIUtils/CLI11/releases/download/v2.6.2/CLI11.hpp
  DOWNLOAD_NO_EXTRACT TRUE
)


FetchContent_Declare(kseq
  URL      https://raw.githubusercontent.com/lh3/seqtk/94e707082d39b0a038f234df676e32d9802c0dc7/kseq.h
  DOWNLOAD_NO_EXTRACT TRUE
)
set(SMHASHER_RAW https://raw.githubusercontent.com/aappleby/smhasher/07bb4de10a63e8cc2e1724865454eba635742383/src)
FetchContent_Declare(murmur3_h
  URL      ${SMHASHER_RAW}/MurmurHash3.h
  DOWNLOAD_NO_EXTRACT TRUE
)
FetchContent_Declare(murmur3_cpp
  URL      ${SMHASHER_RAW}/MurmurHash3.cpp
  DOWNLOAD_NO_EXTRACT TRUE
)
FetchContent_MakeAvailable(cli11 kseq murmur3_h murmur3_cpp)

set(THIRD_PARTY_DIR ${CMAKE_BINARY_DIR}/third_party)
file(MAKE_DIRECTORY ${THIRD_PARTY_DIR})
file(COPY ${cli11_SOURCE_DIR}/CLI11.hpp
          ${kseq_SOURCE_DIR}/kseq.h
          ${murmur3_h_SOURCE_DIR}/MurmurHash3.h
          ${murmur3_cpp_SOURCE_DIR}/MurmurHash3.cpp
     DESTINATION ${THIRD_PARTY_DIR})

add_library(murmurhash3 STATIC ${THIRD_PARTY_DIR}/MurmurHash3.cpp)
target_include_directories(murmurhash3 PUBLIC ${THIRD_PARTY_DIR})
target_compile_options(murmurhash3 PRIVATE -w)

find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)

# Common include paths: SIFT first, then PFSP, then third-party files
add_library(sift_common INTERFACE)
target_include_directories(sift_common INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${PFSP_DIR}/include
  ${THIRD_PARTY_DIR}
  ${sdsl_SOURCE_DIR}/include
)

target_link_libraries(sift_common INTERFACE murmurhash3 sdsl ZLIB::ZLIB)

# -------------------------------------------------------------------------
# sift
# -------------------------------------------------------------------------
add_executable(sift src/main.cpp)
target_link_libraries(sift PRIVATE sift_common OpenMP::OpenMP_CXX)

