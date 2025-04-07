cmake_minimum_required(VERSION 3.10)
project(DistributedHashMapTest CXX)

find_package(UPCXX REQUIRED)

add_executable(distributed_hashmap_test
  distributed_hashmap_test.cpp
  distributed_hashmap.hpp
  distributed_hash.hpp
)

target_link_libraries(distributed_hashmap_test PUBLIC UPCXX::upcxx)

configure_file(run_it.sh run_it.sh COPYONLY)