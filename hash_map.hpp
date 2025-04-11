#pragma once

#include "kmer_t.hpp"
#include <upcxx/upcxx.hpp>

class HashMap {
public:
    upcxx::global_ptr<kmer_pair> global_data;
    upcxx::global_ptr<int> global_used;
    upcxx::atomic_domain<int> atomic_used;
    size_t my_size;

    HashMap(size_t size);
    
    bool insert(const kmer_pair& kmer);
    bool find(const pkmer_t& key_kmer, kmer_pair& val_kmer);
    bool request_slot(uint64_t slot);
    bool slot_used(uint64_t slot);
};

HashMap::HashMap(size_t size)
    : atomic_used({upcxx::atomic_op::compare_exchange}) {
    my_size = size / upcxx::rank_n();
    global_data = upcxx::new_array<kmer_pair>(my_size);
    global_used = upcxx::new_array<int>(my_size);

    // Initialize memory to avoid segmentation faults
    for (size_t i = 0; i < my_size; i++) {
        upcxx::rput(kmer_pair(), global_data + i).wait();
        upcxx::rput(0, global_used + i).wait();
    }
}

bool HashMap::request_slot(uint64_t slot) {
    int expected = 0;
    int desired = 1;
    return atomic_used.compare_exchange(global_used + slot, expected, desired, std::memory_order_relaxed).wait() == expected;
}

bool HashMap::insert(const kmer_pair& kmer) {
    uint64_t hash = kmer.hash();
    int target_rank = hash % upcxx::rank_n();
    uint64_t local_slot = hash % my_size;

    upcxx::barrier(); // Ensure memory consistency before insertion

    return upcxx::rpc(target_rank, [this](uint64_t slot_index, kmer_pair km) {
        upcxx::rput(km, global_data + slot_index).wait();
        return true;
    }, local_slot, kmer).wait();
}

bool HashMap::find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {
    uint64_t hash = key_kmer.hash();
    int target_rank = hash % upcxx::rank_n();
    uint64_t local_slot = hash % my_size;

    upcxx::barrier(); // Synchronize before lookup

    val_kmer = upcxx::rpc(target_rank, [this](uint64_t slot_index) {
        return upcxx::rget(global_data + slot_index).wait();
    }, local_slot).wait();

    if (val_kmer.kmer_str().empty()) {
        return false;  // Handle lookup failure gracefully
    }

    return (val_kmer.kmer == key_kmer);
}
