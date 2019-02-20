#include <array>
#include <cstdint>

#ifndef __HASH__
#	define __HASH__

std::array<uint64_t, 2>
hash64(const uint8_t* data, std::size_t len);

uint64_t
nthHash(uint8_t n, uint64_t hashA, uint64_t hashB, uint64_t filterSize);

#endif
