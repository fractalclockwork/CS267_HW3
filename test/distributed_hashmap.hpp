#pragma once

#include <upcxx/upcxx.hpp>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include "distributed_hash.hpp"

//======================================================
// DistributedHashMap: higher-level abstraction (with batching and caching)
//======================================================
template<typename Key, typename Value>
class DistributedHashMap {
private:
  // The underlying distributed hash stored in a dist_object.
  upcxx::dist_object< DistributedHash<Key,Value> > d_hash;
  std::unordered_map<Key, Value> local_cache; // local cache for lookups
  int rank_id_;
  int world_size_;
  
  int get_target_rank(const Key &key) const {
    return std::hash<Key>{}(key) % world_size_;
  }
  
public:
  DistributedHashMap(int rank_id, int world_size)
    : rank_id_(rank_id),
      world_size_(world_size),
      d_hash( DistributedHash<Key,Value>(rank_id, world_size) )
  { }
  
  // Batch insert:
  // Group a collection of (Key,Value) pairs by their target rank.
  // For remote groups, perform a single RPC passing all entries.
  upcxx::future<> batch_insert(const std::vector<std::pair<Key, Value>> &entries) {
    std::unordered_map<int, std::vector<std::pair<Key,Value>>> groups;
    for(const auto &entry : entries) {
      int target = get_target_rank(entry.first);
      groups[target].push_back(entry);
      local_cache[entry.first] = entry.second; // update local cache immediately
    }
    
    std::vector<upcxx::future<>> rpc_futures;
    for(auto &grp : groups) {
      int target = grp.first;
      auto &vec = grp.second;
      if(target == rank_id_) {
         // Local insertion.
         for(auto &kv : vec) {
           d_hash->insert(kv.first, kv.second);
         }
      } else {
         // For remote insertion, obtain the global pointer to the underlying DistributedHash
         // Here we use global_ptr() method which is available on dist_object.
         /*
         rpc_futures.push_back(
           upcxx::rpc(target,
             // Remote lambda: receives a global pointer to DistributedHash
             [](upcxx::global_ptr< DistributedHash<Key,Value> > remote_ptr,
                std::vector<std::pair<Key,Value>> batch)
             {
               DistributedHash<Key,Value>* ptr = remote_ptr.local();
               for(auto &kv : batch) {
                 ptr->insert(kv.first, kv.second);
               }
             },
             d_hash.global_ptr(), std::move(vec)
           )
         );
         */
        rpc_futures.push_back(
            upcxx::rpc(target,
                // Pass `d_hash` directly as a dist_object instead of converting it to a global_ptr.
                [](upcxx::dist_object<DistributedHash<Key, Value>> &remote_hash,
                   std::vector<std::pair<Key, Value>> batch) {
                    for (const auto &kv : batch) {
                        remote_hash->insert(kv.first, kv.second);
                    }
                },
                d_hash, std::move(vec) // Correctly passing the dist_object directly
            )
        );

      }
    }
    // Wait (blocking) on all remote RPCs.
    for(auto &f : rpc_futures) {
      f.wait();
    }
    return upcxx::make_future();
  }
  
  // Single insert wrapper.
  void insert(const Key &key, const Value &value) {
    batch_insert({{key, value}}).wait();
  }
  
  // Lookup operation.
  // First checks the local cache; if absent, delegates to the underlying distributed hash.
  upcxx::future<Value> find(const Key &key) {
    auto it = local_cache.find(key);
    if(it != local_cache.end())
      return upcxx::make_future(it->second);
    return d_hash->find(key).then([this, key](Value res) {
      if(!(res == Value{}))
         local_cache[key] = res;
      return res;
    });
  }
  
  // Process pending UPC++ RPC requests.
  void process_requests() {
    d_hash->process_requests();
  }
};

