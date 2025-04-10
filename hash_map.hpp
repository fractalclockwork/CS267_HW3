#pragma once

#include <cstring>          // Needed so that memcmp is available (required by pkmer_t.hpp)
#include "kmer_t.hpp"
#include <vector>
#include <cstdint>
#include <stdexcept>

//------------------------------------------------------------------------------
// HashMap
// This is the serial hash map implementation.
// It preserves the functionality of the original starter code and has 
// been reorganized for clarity in preparation for the DHM implementation.
//------------------------------------------------------------------------------
struct HashMap {
    size_t table_size;             // Global number of slots
    std::vector<kmer_pair> data;   // Storage for kmer_pair values
    std::vector<int> used;         // Occupancy flags (0 = free, 1 = used)

    explicit HashMap(size_t size);

    size_t size() const noexcept;

    // Inserts a kmer_pair into the table.
    // Returns true if insertion succeeded (false if full).
    bool insert(const kmer_pair &kmer);

    // Searches for a key (of type pkmer_t) in the table.
    // If found, sets val_kmer and returns true.
    bool find(const pkmer_t &key_kmer, kmer_pair &val_kmer);

    // Writes a kmer_pair value into a given slot.
    void write_slot(uint64_t slot, const kmer_pair &kmer);

    // Reads and returns the kmer_pair from a given slot.
    kmer_pair read_slot(uint64_t slot);

    // Attempts to mark a slot as used; returns true if the slot was free.
    bool request_slot(uint64_t slot);

    // Returns true if the slot is marked as used.
    bool slot_used(uint64_t slot);
};

//------------------------------------------------------------------------------
// Method Definitions
//------------------------------------------------------------------------------

HashMap::HashMap(size_t size)
    : table_size(size), data(size), used(size, 0) {}

size_t HashMap::size() const noexcept {
    return table_size;
}

bool HashMap::insert(const kmer_pair &kmer) {
    uint64_t hash_val = kmer.hash();
    uint64_t probe = 0;
    while (probe < size()) {
        uint64_t slot = (hash_val + probe++) % size();
        if (request_slot(slot)) {
            write_slot(slot, kmer);
            return true;
        }
    }
    return false;  // Table is full.
}

bool HashMap::find(const pkmer_t &key_kmer, kmer_pair &val_kmer) {
    uint64_t hash_val = key_kmer.hash();
    uint64_t probe = 0;
    while (probe < size()) {
        uint64_t slot = (hash_val + probe++) % size();
        if (slot_used(slot)) {
            kmer_pair candidate = read_slot(slot);
            if (candidate.kmer == key_kmer) {
                val_kmer = candidate;
                return true;
            }
        }
    }
    return false;  // Not found.
}

void HashMap::write_slot(uint64_t slot, const kmer_pair &kmer) {
    data[slot] = kmer;
}

kmer_pair HashMap::read_slot(uint64_t slot) {
    return data[slot];
}

bool HashMap::request_slot(uint64_t slot) {
    if (used[slot] != 0)
        return false;
    used[slot] = 1;
    return true;
}

bool HashMap::slot_used(uint64_t slot) {
    return used[slot] != 0;
}
