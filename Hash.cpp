#include "Hash.hpp"
#include <gatbl/kmer.hpp>

using namespace std;

std::array<uint64_t, 2> hash64(const uint8_t* data, std::size_t len) {
	assume(len <= sizeof(uint64_t), "Can only up to 64bit integers"); // pretty brutal
	uint64_t in = 0;
	assume(in == 0, "wtf");
	memcpy(&in, data, len);

	std::array<uint64_t, 2> out;
	out[0] = gatbl::ReversibleHash()(in);
	out[1] = gatbl::ReversibleHash()(in + 0xDEADBEEF);
	return out;
}

uint64_t nthHash(uint8_t n, uint64_t hashA, uint64_t hashB, uint64_t filterSize) {
	return (hashA + n * hashB) % filterSize;
}
