#include <climits>
#include <cstdint>
#include <memory>

#define WORD_TYPE uint64_t
#define WORD_SIZE (sizeof(WORD_TYPE) * CHAR_BIT)
/* reads bit p from e */
#define bitget(e, p) ((((e)[(p) / WORD_SIZE] >> ((p) % WORD_SIZE))) & 1)
/* sets bit p in e */
#define bitset(e, p) ((e)[(p) / WORD_SIZE] |= ((WORD_TYPE)1 << ((p) % WORD_SIZE)))

class BitSet {
  public:
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
