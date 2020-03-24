#include <iostream>
#include <cstdint>

#ifndef SKC_H
#define SKC_H

class SKC {
public:
	__uint128_t sk;
	uint8_t size;
	uint8_t counts[16] = {0};
	uint8_t minimizer_idx;

	/** Construct a superkmer from one kmer and the minimizer position.
	* @param kmer The unsigned int 64 used to represent the binary kmer.
	* @param mini_idx The minimizer position in the kmer.
	*/
	SKC(const uint64_t, const uint8_t);
	/** Look for the minimizer between the SKC and both of the fwd and rev kmer.
	* Call the right compact function if needed.
	*/
	bool add_kmer(const kmer_full& kmer);

	friend std::ostream & operator << (std::ostream& out, const SKC& skc);


private:
	bool compact_right(const uint64_t);
	bool compact_left(const uint64_t);
	bool is_present(uint64_t, uint64_t);
};

#endif
