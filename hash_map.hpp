#pragma once

#include "kmer_t.hpp"
#include "butil.hpp"
#include <upcxx/upcxx.hpp>

class HashMap {
public:
    // Global arrays for each rank's local portion.
    upcxx::global_ptr<kmer_pair> global_data;
    upcxx::global_ptr<int> global_used;
    // Atomic domain to reserve slots.
    upcxx::atomic_domain<int> atomic_used;
    // Local size (number of slots per rank) and total global table size.
    size_t my_size;
    size_t total_size; 

    // Constructor: 'size' is the total table size across ranks.
    HashMap(size_t size);

    // Distributed insertion with linear probing.
    bool insert(const kmer_pair& kmer);
    // Distributed lookup with linear probing.
    bool find(const pkmer_t& key_kmer, kmer_pair& val_kmer);

    // Attempts to reserve a slot (using atomics); returns true if successful.
    bool request_slot(uint64_t slot);
};

HashMap::HashMap(size_t size)
    : atomic_used({upcxx::atomic_op::compare_exchange})
{
    // Distribute the table evenly across ranks.
    size_t num_ranks = upcxx::rank_n();
    my_size = size / num_ranks;
    total_size = my_size * num_ranks;

    // Allocate local portions.
    global_data = upcxx::new_array<kmer_pair>(my_size);
    global_used = upcxx::new_array<int>(my_size);

    // Initialize local arrays (to avoid reading uninitialized memory).
    for (size_t i = 0; i < my_size; i++) {
        upcxx::rput(kmer_pair(), global_data + i).wait();
        upcxx::rput(0, global_used + i).wait();
    }
}

bool HashMap::request_slot(uint64_t slot) {
    int expected = 0;
    int desired = 1;
    // Reserve the slot atomically. The compare_exchange operation returns the *old* value.
    return atomic_used.compare_exchange(global_used + slot, expected, desired,
                                          std::memory_order_relaxed).wait() == expected;
}

bool HashMap::insert(const kmer_pair& kmer) {
    uint64_t hash = kmer.hash();
    // Compute the initial global slot.
    uint64_t global_slot = hash % total_size;
    // Try up to total_size probes.
    for (size_t i = 0; i < total_size; i++) {
        uint64_t probe_slot = (global_slot + i) % total_size;
        int target_rank = static_cast<int>(probe_slot / my_size);
        uint64_t local_slot = probe_slot % my_size;

        // Perform an RPC to the rank that “owns” this slot.
        bool success = upcxx::rpc(target_rank,
            [this](uint64_t slot_index, kmer_pair km) -> bool {
                // Try to reserve the slot atomically.
                if (request_slot(slot_index)) {
                    // Write the k-mer into the reserved slot.
                    upcxx::rput(km, global_data + slot_index).wait();
                    return true;
                }
                return false;
            },
            local_slot, kmer).wait();

        if (success) {
            BUtil::print("Rank %d: Inserted k-mer %s into global slot %lu (target_rank=%d, local_slot=%lu)\n",
                         upcxx::rank_me(), kmer.kmer_str().c_str(), probe_slot, target_rank, local_slot);
            return true;
        } 
        else {
            // For verbose debugging, print each collision.
            BUtil::print("Rank %d: Slot %lu already in use (probe %lu) for k-mer %s\n",
                         upcxx::rank_me(), probe_slot, i, kmer.kmer_str().c_str());
        }
    }
    // No free slot found after a complete probe.
    BUtil::print("Rank %d: Failed to insert k-mer %s after probing entire table.\n",
                 upcxx::rank_me(), kmer.kmer_str().c_str());
    return false;
}

bool HashMap::find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {
    uint64_t hash = key_kmer.hash();
    uint64_t global_slot = hash % total_size;
    // Linear probe over the entire table.
    for (size_t i = 0; i < total_size; i++) {
        uint64_t probe_slot = (global_slot + i) % total_size;
        int target_rank = static_cast<int>(probe_slot / my_size);
        uint64_t local_slot = probe_slot % my_size;

        kmer_pair retrieved = upcxx::rpc(target_rank,
            [this](uint64_t slot_index) -> kmer_pair {
                return upcxx::rget(global_data + slot_index).wait();
            },
            local_slot).wait();

        // If the slot is empty, then the key is not present.
        if (retrieved.kmer_str().empty()) {
            BUtil::print("Rank %d: Find: Empty slot encountered at global slot %lu while looking for k-mer %s.\n",
                         upcxx::rank_me(), probe_slot, key_kmer.get().c_str());
            return false;
        }
        // If the retrieved k-mer matches the query, we have found the key.
        if (retrieved.kmer == key_kmer) {
            BUtil::print("Rank %d: Found k-mer %s at global slot %lu (probe %lu).\n",
                         upcxx::rank_me(), key_kmer.get().c_str(), probe_slot, i);
            val_kmer = retrieved;
            return true;
        } else {
            // For debugging, print when a slot is occupied by a different k-mer.
            BUtil::print("Rank %d: Probe %lu at global slot %lu: Found %s rather than %s.\n",
                         upcxx::rank_me(), i, probe_slot,
                         retrieved.kmer_str().c_str(), key_kmer.get().c_str());
        }
    }
    BUtil::print("Rank %d: K-mer %s not found after probing entire table.\n",
                 upcxx::rank_me(), key_kmer.get().c_str());
    return false;
}
