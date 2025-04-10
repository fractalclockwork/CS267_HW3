#pragma once

// We add <cstring> here so that pkmer_t.hpp (which we cannot modify)
// finds the declaration for memcmp.
#include <cstring>
#include "kmer_t.hpp"
#include <vector>
#include <cstdint>
#include <stdexcept>

//------------------------------------------------------------------------------
// SerialHashMap
//
// This class implements the original serial hash table behavior.
// It encapsulates the core functionality and collision resolution
// (using linear probing) in one place. It also provides a stub for global
// memory mapping (to be overridden in DHM).
//------------------------------------------------------------------------------
class SerialHashMap {
public:
    // Constructor takes total table size.
    explicit SerialHashMap(size_t size);

    // Returns the number of slots.
    size_t size() const noexcept;

    // Inserts kmer_pair into the table. Returns true if insertion succeeded.
    bool insert(const kmer_pair &kmer);

    // Searches for a k-mer by key. If found, sets val_kmer and returns true.
    bool find(const pkmer_t &key_kmer, kmer_pair &val_kmer);

protected:
    size_t table_size;                // Total number of slots.
    std::vector<kmer_pair> data;        // Storage for kmer_pair values.
    std::vector<int> used;            // Occupancy flags (0 = free, 1 = used).

    // Helper: Computes the next slot using linear probing.
    inline uint64_t next_slot(uint64_t hash_val, uint64_t probe) const {
        return (hash_val + probe) % table_size;
    }

    // For future DHM: Map a global slot index to a local index.
    // In serial mode, the mapping is trivial: owner = 0 and local_index = slot.
    virtual void global_to_local(uint64_t slot, int &owner, uint64_t &local_index) const {
        owner = 0;
        local_index = slot;
    }

    // Writes a kmer_pair into a given slot.
    void write_slot(uint64_t slot, const kmer_pair &kmer);

    // Reads the kmer_pair from a given slot.
    kmer_pair read_slot(uint64_t slot);

    // Attempts to mark a slot as used. Returns true if the slot was free.
    bool request_slot(uint64_t slot);

    // Returns true if the slot is marked as used.
    bool slot_used(uint64_t slot);
};

//------------------------------------------------------------------------------
// Method Definitions
//------------------------------------------------------------------------------

SerialHashMap::SerialHashMap(size_t size)
    : table_size(size), data(size), used(size, 0) {}

size_t SerialHashMap::size() const noexcept {
    return table_size;
}

void SerialHashMap::write_slot(uint64_t slot, const kmer_pair &kmer) {
    data[slot] = kmer;
}

kmer_pair SerialHashMap::read_slot(uint64_t slot) {
    return data[slot];
}

bool SerialHashMap::request_slot(uint64_t slot) {
    if (used[slot] != 0)
        return false;
    used[slot] = 1;
    return true;
}

bool SerialHashMap::slot_used(uint64_t slot) {
    return used[slot] != 0;
}

bool SerialHashMap::insert(const kmer_pair &kmer) {
    uint64_t hash_val = kmer.hash();
    uint64_t probe = 0;
    while (probe < size()) {
        uint64_t slot = next_slot(hash_val, probe);
        if (request_slot(slot)) {
            write_slot(slot, kmer);
            return true;
        }
        ++probe;
    }
    return false;  // Hash table is full.
}

bool SerialHashMap::find(const pkmer_t &key_kmer, kmer_pair &val_kmer) {
    uint64_t hash_val = key_kmer.hash();
    uint64_t probe = 0;
    while (probe < size()) {
        uint64_t slot = next_slot(hash_val, probe);
        if (slot_used(slot)) {
            kmer_pair candidate = read_slot(slot);
            if (candidate.kmer == key_kmer) {
                val_kmer = candidate;
                return true;
            }
        }
        ++probe;
    }
    return false;
}
