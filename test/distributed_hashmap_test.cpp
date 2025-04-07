#include <upcxx/upcxx.hpp>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include "distributed_hashmap.hpp"

//======================================================
// Main: Test driver for the Distributed HashMap
//======================================================
int main() {
  upcxx::init();
  
  int rank = upcxx::rank_me();
  int world_size = upcxx::rank_n();
  
  // Instantiate our DistributedHashMap for int keys and std::string values.
  DistributedHashMap<int, std::string> dhm(rank, world_size);
  
  // Generate 100 test entries.
  std::vector<std::pair<int, std::string>> entries;
  for (int i = 0; i < 100; i++) {
    entries.push_back({i, "Value_" + std::to_string(i)});
  }
  
  // Perform a batch insertion.
  dhm.batch_insert(entries).wait();
  
  upcxx::barrier();
  
  // From rank 0, perform lookups and print the results.
  if (rank == 0) {
    for (int i = 0; i < 100; i++) {
      std::string value = dhm.find(i).wait();
      std::cout << "Key: " << i << " Value: " << value << std::endl;
    }
  }
  
  upcxx::finalize();
  return 0;
}
