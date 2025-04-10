#pragma once

#include <upcxx/upcxx.hpp>
#include <upcxx/dist_object.hpp>
#include "kmer_t.hpp"
#include <vector>
#include <stdexcept>
#include <cstdint>

//------------------------------------------------------------------------------
// DistributedHashMap
//
// A simple synchronous distributed hash map implementing linear probing,
// using UPC++ global memory and PGAS. Each rank allocates a local block
// of entries. Global-to-local mapping is performed via block decomposition.
// Remote accesses are done synchronously using rget/rput and rpc calls.
//------------------------------------------------------------------------------
template<typename Key, typename Value>
class DistributedHashMap {
private:
  // Each slot (entry) in the table.
  struct Entry {
    bool occupied;
    Value value;
    Entry() : occupied(false), value() {}
  };

  size_t global_size_;      // Total number of slots in the table.
  size_t local_size_;       // Number of slots allocated on this rank.
  int world_size_;          // Total number of ranks.
  int rank_;                // My rank.

  // Local table allocated on this rank.
  upcxx::global_ptr<Entry> local_table_;

  // Global vector of pointers to each rank's local table.
  std::vector<upcxx::global_ptr<Entry>> table_ptrs_;

  // Maps a global slot index to its owner and local index via block decomposition.
  void global_to_local(uint64_t slot, int &owner, uint64_t &local_idx) const {
      size_t base = global_size_ / world_size_;
      size_t rem = global_size_ % world_size_;
      uint64_t start = 0;
      for (int r = 0; r < world_size_; r++) {
          size_t block_size = (static_cast<size_t>(r) < rem) ? base + 1 : base;
          if (slot < start + block_size) {
              owner = r;
              local_idx = slot - start;
              return;
          }
          start += block_size;
      }
      throw std::runtime_error("global_to_local: slot index out of range");
  }

public:
  // Constructor: Allocate the distributed table.
  DistributedHashMap(size_t table_size)
      : global_size_(table_size)
  {
      world_size_ = upcxx::rank_n();
      rank_ = upcxx::rank_me();

      // Determine local block size using block decomposition.
      size_t base = global_size_ / world_size_;
      size_t rem  = global_size_ % world_size_;
      local_size_ = (static_cast<size_t>(rank_) < rem) ? base + 1 : base;

      // Allocate local memory for entries.
      local_table_ = upcxx::new_array<Entry>(local_size_);
      // Initialize each entry to "unoccupied."
      for (size_t i = 0; i < local_size_; i++) {
          upcxx::rput(Entry{}, local_table_ + i).wait();
      }

      // Use a dist_object to publish each rank’s local_table_ pointer.
      upcxx::dist_object<upcxx::global_ptr<Entry>> table_obj(local_table_);
      upcxx::barrier(); // Ensure all ranks have created their dist_object.

      // Gather each rank's local pointer.
      table_ptrs_.resize(world_size_);
      for (int r = 0; r < world_size_; r++) {
          table_ptrs_[r] = table_obj.fetch(r).wait();
      }
      upcxx::barrier(); // Ensure all pointers are gathered.
  }

  // Returns the global number of slots.
  size_t size() const { return global_size_; }

  // Synchronous insert: Attempts to insert (key, value) using linear probing.
  // Returns true if insertion succeeds; false if the table is full.
  bool insert(const Key &key, const Value &value) {
      uint64_t hash_val = key.hash();
      uint64_t probe = 0;
      while (probe < global_size_) {
          uint64_t slot = (hash_val + probe) % global_size_;
          int owner;
          uint64_t local_idx;
          global_to_local(slot, owner, local_idx);
          bool success = false;
          if (owner == rank_) {
              // Local insertion.
              Entry entry = upcxx::rget(local_table_ + local_idx).wait();
              if (!entry.occupied) {
                  entry.occupied = true;
                  entry.value = value;
                  upcxx::rput(entry, local_table_ + local_idx).wait();
                  success = true;
              }
          } else {
              // Remote insertion via RPC.
              success = upcxx::rpc(owner,
                  [=](upcxx::global_ptr<Entry> remote_table, uint64_t idx, const Value &val) -> bool {
                      Entry entry = upcxx::rget(remote_table + idx).wait();
                      if (!entry.occupied) {
                          entry.occupied = true;
                          entry.value = val;
                          upcxx::rput(entry, remote_table + idx).wait();
                          return true;
                      }
                      return false;
                  },
                  table_ptrs_[owner], local_idx, value).wait();
          }
          if (success)
              return true;
          ++probe;
      }
      return false;  // Table is full.
  }

  // Synchronous find: Searches for the given key.
  // If found, sets result and returns true; otherwise returns false.
  bool find(const Key &key, Value &result) {
      uint64_t hash_val = key.hash();
      uint64_t probe = 0;
      while (probe < global_size_) {
          uint64_t slot = (hash_val + probe) % global_size_;
          int owner;
          uint64_t local_idx;
          global_to_local(slot, owner, local_idx);
          Entry entry;
          if (owner == rank_) {
              entry = upcxx::rget(local_table_ + local_idx).wait();
          } else {
              entry = upcxx::rpc(owner,
                  [=](upcxx::global_ptr<Entry> remote_table, uint64_t idx) -> Entry {
                      return upcxx::rget(remote_table + idx).wait();
                  },
                  table_ptrs_[owner], local_idx).wait();
          }
          if (entry.occupied && entry.value.kmer == key) {
              result = entry.value;
              return true;
          }
          ++probe;
      }
      return false;
  }
};
