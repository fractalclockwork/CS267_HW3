#pragma once
#include <upcxx/upcxx.hpp>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>

//======================================================
// DistributedHash: the low-level, templated distributed hash table
//======================================================
template<typename Key, typename Value>
class DistributedHash {
public:
  using HashTable = std::unordered_map<Key, Value>;
  using DistMap = upcxx::dist_object<HashTable>;
  
private:
  DistMap local_data; // local storage on each rank
  int rank_id_;
  int world_size_;
  
public:
  DistributedHash(int rank_id, int world_size)
    : rank_id_(rank_id),
      world_size_(world_size),
      local_data(HashTable{}) // Construct with an empty hash table
  { }
  
  // Insert a key/value pair.
  // If the key hashes to the local rank, insert it locally;
  // otherwise, use an RPC (fire-and-forget) remote insertion.
  void insert(const Key &key, const Value &value) {
    int target = get_target_rank(key);
    if(target == rank_id_) {
      local_data->emplace(key, value);
    } else {
      upcxx::rpc_ff(target,
         [](upcxx::dist_object<HashTable> &dobj, const Key &key, const Value &value) {
           dobj->emplace(key, value);
         },
         local_data, key, value);
    }
  }
  
  // Lookup: if the key exists locally, return it immediately; else, perform a remote RPC lookup.
  upcxx::future<Value> find(const Key &key) {
    int target = get_target_rank(key);
    if(target == rank_id_) {
      auto it = local_data->find(key);
      if(it != local_data->end())
         return upcxx::make_future(it->second);
      else
         return upcxx::make_future(Value{}); // default-constructed value
    } else {
      return upcxx::rpc(target,
         [](upcxx::dist_object<HashTable>& dobj, const Key &key) -> Value {
           auto it = dobj->find(key);
           return (it != dobj->end()) ? it->second : Value{};
         },
         local_data, key);
    }
  }
  
  // Process pending UPC++ RPC requests.
  void process_requests() {
    upcxx::progress(upcxx::progress_level::user);
  }
  
private:
  int get_target_rank(const Key &key) const {
    return std::hash<Key>{}(key) % world_size_;
  }
};

