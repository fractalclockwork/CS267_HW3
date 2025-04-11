#pragma once
#include <upcxx/upcxx.hpp>
#include <unordered_map>
#include <string>
#include <vector>
#include <utility>
#include "kmer_t.hpp"

// DistributedHashMap is a nontrivial implementation that partitions 
// the key-space by having each rank “own” a portion of the hash space.
// Instead of issuing an RPC per insertion, we batch remote updates.
class DistributedHashMap {
public:
  // Type aliases for clarity.
  using local_map_type = std::unordered_map<std::string, kmer_pair>;
  using kv_pair = std::pair<std::string, kmer_pair>;
  using batch_type = std::vector<kv_pair>;
  
private:
  // Each rank holds a local copy, wrapped in a UPC++ dist_object.
  upcxx::dist_object<local_map_type> local_map;
  size_t table_size_;
  int rank_id_;
  int world_size_;

  // Partition function: each rank owns keys whose hash falls in its interval.
  // (For simplicity, we still use key % world_size.)
  int get_target_rank(const std::string &key) const {
    return std::hash<std::string>{}(key) % world_size_;
  }

  // Local insertion: add a key-value pair to the local hash table.
  void insert_locally(const std::string &key, const kmer_pair &value) {
    (*local_map)[key] = value;
  }

  // Remote insertion: insert a batch of updates via a single RPC.
  void insert_batch_remote(int target_rank, const batch_type &batch) {
    upcxx::rpc(target_rank,
      [](upcxx::dist_object<local_map_type> &lmap, const batch_type &batch) {
        for (const auto &entry : batch) {
          (*lmap)[entry.first] = entry.second;
        }
      },
      local_map, batch).wait();
  }

public:
  // Constructor. Each rank initializes its local hash table.
  DistributedHashMap(size_t table_size, int rank_id, int world_size)
      : table_size_(table_size), rank_id_(rank_id), world_size_(world_size),
        local_map({}) {}

  // Batch insert: partition input items by owner and update with one RPC per target.
  void insert_all(const std::vector<kmer_pair> &items) {
    // Partition batch: map target_rank -> vector of key-value pairs.
    std::unordered_map<int, batch_type> batches;
    for (const auto &item : items) {
      std::string key = item.kmer_str();
      int target = get_target_rank(key);
      batches[target].push_back({key, item});
    }
    // Issue batch updates to each target.
    for (const auto &pair : batches) {
      int target = pair.first;
      const batch_type &batch = pair.second;
      if (target == rank_id_) {
        // Insert locally.
        for (const auto &entry : batch) {
          insert_locally(entry.first, entry.second);
        }
      }
      else {
        // Insert remotely in one call.
        insert_batch_remote(target, batch);
      }
    }
    // Synchronize so that all updates are visible.
    upcxx::barrier();
  }

  // Distributed find: look up a key on the owning rank.
  bool find(const std::string &key, kmer_pair &result) {
    int target = get_target_rank(key);
    if(target == rank_id_){
      auto it = local_map->find(key);
      if(it != local_map->end()){
        result = it->second;
        return true;
      }
      return false;
    }
    else {
      auto fut = upcxx::rpc(target,
         [](upcxx::dist_object<local_map_type> &lmap, const std::string &key) -> upcxx::future<kmer_pair> {
             auto it = lmap->find(key);
             return upcxx::make_future((it != lmap->end()) ? it->second : kmer_pair());
         },
         local_map, key);
      kmer_pair found = fut.wait();
      if(!found.kmer_str().empty()){
         result = found;
         return true;
      }
      return false;
    }
  }

  // Process requests and synchronize.
  void process_requests() {
    upcxx::progress(upcxx::progress_level::user);
    upcxx::barrier();
  }
};
