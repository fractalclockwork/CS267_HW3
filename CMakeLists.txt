cmake_minimum_required(VERSION 3.14)
project(hw3 LANGUAGES CXX)

find_package(UPCXX REQUIRED)

# Group number
set(GROUP_NAME "None" CACHE STRING "Your group name as it appears on bCourses (no spaces)")

# Packaging system
if (NOT ${GROUP_NAME} STREQUAL None)
    set(CPACK_GENERATOR TGZ)
    set(CPACK_PACKAGE_FILE_NAME "cs267${GROUP_NAME}_hw3")
    install(FILES kmer_hash.cpp hash_map.hpp DESTINATION .)
    install(FILES ${CPACK_PACKAGE_FILE_NAME}.pdf DESTINATION .)
    include(CPack)
endif ()

# Build the kmer_hash executable
add_executable(kmer_hash_19 kmer_hash.cpp)
target_link_libraries(kmer_hash_19 PRIVATE UPCXX::upcxx)
target_compile_definitions(kmer_hash_19 PRIVATE "KMER_LEN=19")

# Build the kmer_hash executable
add_executable(kmer_hash_51 kmer_hash.cpp)
target_link_libraries(kmer_hash_51 PRIVATE UPCXX::upcxx)
target_compile_definitions(kmer_hash_51 PRIVATE "KMER_LEN=51")

# Copy the job scripts
#configure_file(job-perlmutter-starter job-perlmutter-starter COPYONLY)
#configure_file(check_it.sh check_it.sh COPYONLY)
#configure_file(work_it.sh work_it.sh COPYONLY)

# Glob Copy the job scripts
file(GLOB SH_FILES "${CMAKE_SOURCE_DIR}/scripts/*.sh")
file(GLOB JOB_FILES "${CMAKE_SOURCE_DIR}/scripts/job-*")

foreach(JOB_FILE ${JOB_FILES})
  configure_file(${JOB_FILE} ${CMAKE_BINARY_DIR} COPYONLY)
endforeach()

foreach(SH_FILE ${SH_FILES})
  configure_file(${SH_FILE} ${CMAKE_BINARY_DIR} COPYONLY)
endforeach()