#include <cstdint>

#define WORD_TYPE uint64_t
#define WORD_SIZE (sizeof(WORD_TYPE) * 8)
/* reads bit p from e */
#define bitget(e, p) ((((e)[(p) / WORD_SIZE] >> ((p) % WORD_SIZE))) & 1)
/* sets bit p in e */
#define bitset(e, p) ((e)[(p) / WORD_SIZE] |= ((WORD_TYPE)1 << ((p) % WORD_SIZE)))

class BitSet
{
  public:
	BitSet(uint64_t nb_bits);
	~BitSet();

	bool get(uint64_t pos) const;
	void set(uint64_t pos);

	void reset();

	uint64_t size() const;

  private:
	uint64_t* m_bits;
	uint64_t m_nb_bits;
	uint64_t m_int_size;
};
