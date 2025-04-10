#pragma once

#include <cstring>          // For memcmp()
#include "kmer_t.hpp"
#include <vector>
#include <cstdint>
#include <stdexcept>

//------------------------------------------------------------------------------
// HashMap
// A serial hash map implementation that stores kmer_pair values.
struct HashMap {
    size_t table_size;               // Total number of slots.
    std::vector<kmer_pair> data;     // Storage for kmer_pair values.
    std::vector<int> used;           // Occupancy flags (0 = free, 1 = used)

    explicit HashMap(size_t size);

    size_t size() const noexcept;

    // Inserts a kmer_pair into the table.
    // Returns true if insertion succeeded, false otherwise.
    bool insert(const kmer_pair &kmer);

    // Searches for a k-mer (key) in the table.
    // On success, sets val_kmer and returns true.
    bool find(const pkmer_t &key, kmer_pair &val_kmer) const;

    // Writes the kmer_pair into the given slot.
    void write_slot(uint64_t slot, const kmer_pair &kmer);

    // Reads and returns the kmer_pair from the given slot.
    kmer_pair read_slot(uint64_t slot) const;

    // Attempts to claim a slot; returns true on success.
    bool request_slot(uint64_t slot);

    // Returns true if the specified slot is already used.
    bool slot_used(uint64_t slot) const;
};

//------------------------------------------------------------------------------
// Method definitions
//------------------------------------------------------------------------------

HashMap::HashMap(size_t size)
    : table_size(size), data(size), used(size, 0) {}

size_t HashMap::size() const noexcept {
    return table_size;
}

void HashMap::write_slot(uint64_t slot, const kmer_pair &kmer) {
    data[slot] = kmer;
}

// Marked const.
kmer_pair HashMap::read_slot(uint64_t slot) const {
    return data[slot];
}

bool HashMap::request_slot(uint64_t slot) {
    if (used[slot] != 0)
        return false;
    used[slot] = 1;
    return true;
}

// Marked const.
bool HashMap::slot_used(uint64_t slot) const {
    return used[slot] != 0;
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

bool HashMap::find(const pkmer_t &key, kmer_pair &val_kmer) const {
    uint64_t hash_val = key.hash();
    uint64_t probe = 0;
    while (probe < size()) {
        uint64_t slot = (hash_val + probe++) % size();
        if (slot_used(slot)) {
            kmer_pair candidate = read_slot(slot);
            if (candidate.kmer == key) {  // Assumes operator== (using memcmp) is correctly defined.
                val_kmer = candidate;
                return true;
            }
        }
    }
    return false;  // Not found.
}
