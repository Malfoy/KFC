#include <climits>
#include <cstdint>
#include <memory>

class BitSet {
  public:
	using word_type = uint64_t;
	static constexpr size_t word_bitwidth = sizeof(word_type) * CHAR_BIT;

	BitSet(uint64_t nb_bits);

	bool get(uint64_t pos) const;
	void set(uint64_t pos);
	bool get_and_set(uint64_t pos);

	void reset();

	uint64_t size() const;

  private:
	std::unique_ptr<uint64_t[]> m_bits;
	uint64_t m_nb_bits;
	uint64_t m_int_size;
};
