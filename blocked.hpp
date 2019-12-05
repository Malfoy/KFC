#ifndef BLOCKEDBLOOM_HPP
#define BLOCKEDBLOOM_HPP

#include <cstdint>
#include <climits>
#include <cstdlib>
#include <cstring>

#include <memory>
#include <exception>
#include <iostream>

#include "MurmurHash3.h"
#undef assert // FIXME upstream
#include <gatbl/common.hpp>

/**
 * Type traits for mutating integers of small bitwidth
 * The bitwidth can be as small as 1bit and as large as the bitwidth of the word representation, but must be constant.
 * This is the generic implementation with bit masking.
 */
template<size_t _bit_width, typename _word_t = uint64_t, typename = void>
struct subword_traits {
	using word_t = _word_t;
	static constexpr size_t subword_bits = _bit_width;
	static constexpr size_t word_bits = CHAR_BIT * sizeof(word_t);
	static constexpr size_t nsubwords = word_bits / subword_bits;
	static_assert(subword_bits < word_bits, "Subword should be smaller than the word");
	static_assert(word_bits % subword_bits == 0, "Subword bitwidth must divide the word width");

	/// A reference to a constant bit slice
	struct const_reference {
		const_reference(const word_t& word, size_t subword_idx)
		  : m_word(word)
		  , m_offset(static_cast<uint8_t>(subword_bits * subword_idx)) {
			assume(subword_idx < nsubwords, "Access out of bound of word");
		}

	  public:
		operator word_t() const { return (m_word >> m_offset) & m_mask; }

	  protected:
		const word_t& m_word;
		const uint8_t m_offset;
	};

	/// A reference to a bit slice
	struct reference {
		reference(word_t& word, size_t subword_idx)
		  : m_word(word)
		  , m_offset(static_cast<uint8_t>(subword_bits * subword_idx)) {
			assume(subword_idx < nsubwords, "Access out of bound of word");
		}

		reference& operator=(word_t val) {
			assume(val <= m_mask, "Value %llu is not representable in %llu bits", val, word_bits);
			m_word &= ~(m_mask << m_offset);
			m_word |= val << m_offset;
			return *this;
		}

		operator word_t() const { return (m_word >> m_offset) & m_mask; }

	  protected:
		word_t& m_word;
		const uint8_t m_offset;
	};

	static reference at(word_t& word, size_t idx) { return {word, idx}; }
	static const_reference at(const word_t& word, size_t idx) { return {word, idx}; }

  private:
	static constexpr word_t m_mask = (word_t(1) << subword_bits) - 1;
};

/// The trivial case where subword == word
template<size_t _bit_width, typename _word_t>
struct subword_traits<_bit_width, _word_t, std::enable_if<_bit_width == sizeof(_word_t) * CHAR_BIT>> {
	using word_t = _word_t;
	using reference = word_t&;
	using const_reference = const word_t&;

	static reference at(word_t& word, size_t idx) { return word; }
	static const_reference at(const word_t& word, size_t idx) { return word; }

	static constexpr size_t word_bits = CHAR_BIT * sizeof(word_t);
	static constexpr size_t subword_bits = word_bits;
	static constexpr size_t nsubwords = 1;
};

// TODO: 8 bit version

#ifndef __cpp_lib_hardware_interference_size
static constexpr size_t hardware_constructive_interference_size = 64;
#endif

/// A cacheline of integers mutated through a subword_traits
template<typename _subword_traits, size_t block_bytes = hardware_constructive_interference_size>
struct Block : public _subword_traits {
	using typename _subword_traits::const_reference;
	using typename _subword_traits::reference;
	using typename _subword_traits::word_t;
	static constexpr size_t nwords = block_bytes / sizeof(word_t);
	static constexpr size_t nsubwords = _subword_traits::nsubwords * nwords;

	const_reference at(size_t idx) const {
		assume(idx <= nsubwords, "Access out of bound in block %llu > %llu", idx, nwords);
		return const_reference(m_arr[idx / _subword_traits::nsubwords], idx % _subword_traits::nsubwords);
	}
	reference at(size_t idx) {
		assume(idx <= nsubwords, "Access out of bound in block %llu > %llu", idx, nwords);
		return reference(m_arr[idx / _subword_traits::nsubwords], idx % _subword_traits::nsubwords);
	}

  private:
	word_t m_arr[nwords] = {0};
};

// Before C++17 new doesn't correctly align over-aligned type, so we use posix_memalign and a free() based deallocator
template<typename T>
struct FreeDelete {
	void operator()(T* ptr) const {
		ptr->~T();
		free(ptr);
	}
};

template<typename T>
using aligned_array = std::unique_ptr<T[], FreeDelete<T>>;

template<typename T>
aligned_array<T> make_unique_aligned(size_t n) {
	void* void_ptr;
	if (posix_memalign(&void_ptr, alignof(T), n * sizeof(T)) != 0) {
		throw std::bad_alloc();
	}
	auto ptr = aligned_array<T>(static_cast<T*>(void_ptr));
	// Initialize the items with placement new and default constructor
	for (unsigned i = 0; i < n; i++) {
		new (&ptr[i]) T{};
	}

	return ptr;
}

/**
 * @brief A vector of overaligned items
 */
template<typename T>
class AlignedVector {
  public:
	using element_type = T;

	AlignedVector(uint64_t n)
	  : m_size(n) {
		m_arr = make_unique_aligned<element_type>(n);
	}

	size_t size() const { return m_size; }
	size_t nbytes() const { return sizeof(element_type) * m_size; }

	const element_type& at(size_t i) const {
		assume(i < m_size, "indice out of range: %llu >= %llu", i, m_size);
		return m_arr[i];
	}
	element_type& at(size_t i) {
		assume(i < m_size, "indice out of range: %llu >= %llu", i, m_size);
		return m_arr[i];
	}

  private:
	aligned_array<element_type> m_arr;
	const size_t m_size;
};

template<uint8_t counter_bits = 2, uint8_t hash_funcs = 3>
class BlockedCMS {
	using subword_trait = subword_traits<counter_bits>;
	using block_t = Block<subword_trait>;
	using array_t = AlignedVector<block_t>;
	using word_t = typename subword_trait::word_t;
    using hash_t = __uint128_t;

	static constexpr size_t ncounter_per_block = block_t::nsubwords;

  public:
	static constexpr word_t max_count = (word_t(1) << counter_bits) - 1;
	static constexpr size_t block_bytes = sizeof(block_t);

	/// @param log2size: log2 of the number of block
	BlockedCMS(uint8_t log2size)
	  : m_arr(size_t(1) << log2size)
      , m_log2size(log2size) {
        assume(log2size * hash_funcs <= CHAR_BIT * sizeof(hash_t), "not enough bits in hash type");
    }

	word_t aproximate(void* data, size_t len) {
        hash_t hash0;
		MurmurHash3_x64_128(data, static_cast<int>(len), 0, &hash0);
        this->aproximate(hash0);
    }

    word_t aproximate(hash_t hash) {
        block_t& block = m_arr.at(static_cast<size_t>(hash) % (size_t(1) << m_log2size));
        hash >>= m_log2size;

		auto min_val = ~word_t(0);
		for (unsigned i = 0; i < hash_funcs; i++) {
			size_t counter_idx = hash % ncounter_per_block;
			hash /= ncounter_per_block;

			word_t val = block.at(counter_idx * counter_bits, counter_bits);
			min_val = val < min_val ? val : min_val;
		}

		return min_val;
	}

    /**
     * Increment the counter
     * @returns the previous value
     */
    word_t add(void* data, size_t len) {
        hash_t hash0;
        MurmurHash3_x64_128(data, static_cast<int>(len), 0, &hash0);

        return this->add(hash0);
    }

    /**
     * Increment the counter
     * @returns the previous value
     */
    word_t add(hash_t hash0) {
        block_t& block = m_arr.at(static_cast<size_t>(hash0) % (size_t(1) << m_log2size));
        hash0 >>= m_log2size;

        hash_t hash = hash0;
        auto min_val = ~word_t(0);
        for (unsigned i = 0; i < hash_funcs; i++) {
            word_t val = block.at(hash % ncounter_per_block);
            hash /= ncounter_per_block;

            min_val = val < min_val ? val : min_val;
        }

        // Max count achieved by minimum : reset all associated counters
        if (min_val == max_count) {
            hash = hash0;
            for (unsigned i = 0; i < hash_funcs; i++) {
                auto ref = block.at(hash % ncounter_per_block);
                hash /= ncounter_per_block;
                ref = 0;
            }
        } else { // Increment the counters equal to the minimum
            hash = hash0;
            for (unsigned i = 0; i < hash_funcs; i++) {
                auto ref = block.at(hash % ncounter_per_block);
                hash /= ncounter_per_block;

                word_t val = ref;
                if (val == min_val) ref = val + 1;
                // else if(val > min_val) ref = val - 1;
            }
        }

        return min_val;
    }

    std::array<size_t, max_count + 1> histogram() const {
		std::array<size_t, max_count + 1> res{};
		for (size_t i = 0; i < m_arr.size(); i++) {
			auto block = m_arr.at(i);
			for (size_t j = 0; j < ncounter_per_block; j++) {
				word_t val = block.at(j);
				res[val]++;
			}
		}
		return res;
	}

    size_t size() const { return m_arr.size() * ncounter_per_block; }

    friend std::ostream& operator<<(std::ostream& out, const BlockedCMS& that) {
        std::array<size_t, max_count + 1> hist = that.histogram();
        size_t ncounters = that.size();
		for (size_t i = 0; i <= max_count; i++)
            out << "\tP(x=" << i << ") =\t" << hist[i] << '/' << ncounters << "\t= " << double(hist[i]) / ncounters << std::endl;
		return out;
	}

  private:
	array_t m_arr;
	const uint8_t m_log2size;
};

template<uint8_t hash_funcs = 3>
class BlockedBloom {
	using subword_trait = subword_traits<1>;
	using block_t = Block<subword_trait>;
	using array_t = AlignedVector<block_t>;
	using word_t = typename subword_trait::word_t;
    using hash_t = __uint128_t;

	static constexpr size_t nbits_per_block = block_t::nsubwords;

  public:
	static constexpr size_t block_bytes = sizeof(block_t);

	BlockedBloom(uint8_t log2size)
	  : m_arr(size_t(1) << log2size)
	  , m_bits_set(0)
      , m_log2size(log2size) {
        assume(log2size * hash_funcs <= CHAR_BIT * sizeof(hash_t), "not enough bits in hash type");
    }

    /**
     * Add an element
     * @returns the previous value
     */
    word_t add(void* data, size_t len) {
        hash_t hash0;
        MurmurHash3_x64_128(data, static_cast<int>(len), 0, &hash0);

        return this->add(hash0);
    }

	/**
	 * Add an element
	 * @returns the previous value
	 */
    bool add(hash_t hash0) {
		block_t& block = m_arr.at(static_cast<size_t>(hash0) % (size_t(1) << m_log2size));
		hash0 >>= m_log2size;

        hash_t hash = hash0;
		bool was_present = true;
		for (unsigned i = 0; i < hash_funcs; i++) {
			auto ref = block.at(hash % nbits_per_block);
			hash /= nbits_per_block;

			bool was_set = ref;
			was_present &= was_set;
			m_bits_set += !was_set;
			ref = 1;
		}

		if (was_present) {
			hash = hash0;
			for (unsigned i = 0; i < hash_funcs; i++) {
				block.at(hash % nbits_per_block) = 0;
				hash /= nbits_per_block;
			}
			assume(m_bits_set >= hash_funcs, "Negative number of bits ?");
			m_bits_set -= hash_funcs;
		}

		return was_present;
	}

	bool possiblyContains(void* data, size_t len) {
        hash_t hash0;
		MurmurHash3_x64_128(data, static_cast<int>(len), 0, &hash0);
        return this->possiblyContains(hash0);
    }

    bool possiblyContains(hash_t hash0) {
		block_t& block = m_arr.at(static_cast<size_t>(hash0) % (size_t(1) << m_log2size));
		hash0 >>= m_log2size;

        hash_t hash = hash0;
		bool present = true;
		for (unsigned i = 0; i < hash_funcs; i++) {
			present &= block.at(hash % nbits_per_block);
			hash /= nbits_per_block;
		}

		return present;
	}

    size_t size() const { return m_arr.size() * nbits_per_block; }

	void reset() {
		m_bits_set = 0;
		memset(&m_arr.at(0), 0, block_bytes << m_log2size);
	}

	bool add_resetting(void* data, size_t len, double reset_ratio) {
		bool was_present = this->add(data, len);
		if (not was_present && this->nbBitsSet() >= reset_ratio * this->size()) {
			this->reset();
		}
		return was_present;
	}

	size_t nbBitsSet() const { return m_bits_set; }

    friend std::ostream& operator<<(std::ostream& out, const BlockedBloom& that) {
        size_t set = that.nbBitsSet();
        size_t size = that.size();
        out << "\tP(x=1) =\t" << set << '/' << size << "\t= " << double(set) / size << std::endl;
        return out;
    }

  private:
	array_t m_arr;
	size_t m_bits_set;
	const uint8_t m_log2size;
};

template<uint8_t num_blooms = 4, uint8_t hash_funcs = 3>
class LayeredBlockedBloom {
  public:
	using bloom_t = BlockedBloom<hash_funcs>;
	static constexpr size_t block_bytes = bloom_t::block_bytes;

	LayeredBlockedBloom(std::array<size_t, num_blooms> sizes, double reset_ratio = 0.65)
	  : m_blooms(sizes)
	  , m_reset_ratio(reset_ratio) {}

	bool insert(void* data, size_t len) {
		for (auto& level : m_blooms) {
			if (!level.add_resetting(data, len, this->m_reset_ratio)) {
				return false;
			}
		}
		return true;
	}

  private:
	std::array<bloom_t, num_blooms> m_blooms;
	double m_reset_ratio;
};

#endif // BLOCKEDBLOOM_HPP
