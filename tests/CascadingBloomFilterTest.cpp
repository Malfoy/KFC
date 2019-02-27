#include "CascadingBloomFilter.hpp"
#include "lest.hpp"
#include <array>
#include <iostream>
#include <vector>

using namespace std;

// clang-format off
const lest::test module[] = {
	CASE("Setting and querying BloomFilter") {
		SETUP("A CascadingBloomFilter with 256 bytes composed of 3 BF") {
			const unsigned byte_length = 256;
			CascadingBloomFilter cbf = CascadingBloomFilter(byte_length, 3, .5);

			EXPECT(cbf.size() == byte_length);

			
			SECTION("Setting and querying sizes") {
				uint8_t val[] = {0};

				// must contain 0s only
				vector<uint64_t> sizes = cbf.filter_sizes();
				EXPECT(sizes[0] == 0);
				EXPECT(sizes[1] == byte_length*8/2);
				EXPECT(sizes[2] == 0);
				EXPECT(sizes[3] == byte_length*8/4);
				EXPECT(sizes[4] == 0);
				EXPECT(sizes[5] == byte_length*8/4);

				cbf.insert(val, 1);
				// Must contain 1s only in the 1st BF
				sizes = cbf.filter_sizes();
				EXPECT(sizes[0] == 2);
				EXPECT(sizes[2] == 0);
				EXPECT(sizes[4] == 0);

				cbf.insert(val, 1);
				//                 in the 1st & 2nd BFs
				sizes = cbf.filter_sizes();
				EXPECT(sizes[0] == 2);
				EXPECT(sizes[2] == 2);
				EXPECT(sizes[4] == 0);

				cbf.insert(val, 1);
				//                 in the 1st & 2nd & 3rd BFs
				sizes = cbf.filter_sizes();
				EXPECT(sizes[0] == 2);
				EXPECT(sizes[2] == 2);
				EXPECT(sizes[4] == 2);

				cbf.insert(val, 1);
				// BFs must not change
				sizes = cbf.filter_sizes();
				EXPECT(sizes[0] == 2);
				EXPECT(sizes[1] == byte_length*8/2);
				EXPECT(sizes[2] == 2);
				EXPECT(sizes[3] == byte_length*8/4);
				EXPECT(sizes[4] == 2);
				EXPECT(sizes[5] == byte_length*8/4);
			}

		// SECTION("Reset when setting 50% different bits") {
		// 		uint8_t val[] = {0};
		// 		unsigned i = 0;
		// 		while (bf.nbBitsSet() < (unsigned) (byte_length*8/2-1)) {
		// 			vxal[0] = i++;
		// 			bf.add(val, 1);
		// 		}
		// 		uint64_t current_nb_bits_set = bf.nbBitsSet();
		// 		// Next inserted bit should induce a reset
		// 		// But next inserted bit may fall on an existing set bit (proba: 1/2) .
		// 		while (bf.nbBitsSet() == current_nb_bits_set) {
		// 			val[0] = i++;
		// 			bf.add(val, 1);
		// 		}

		// 		// Now the BF should have been reset
		// 		EXPECT(bf.nbBitsSet() == (unsigned)0);
		// 		for (unsigned i = 0; i < 255; i++) {
		// 			val[0] = i;
		// 			EXPECT(! bf.possiblyContains(val, 1));
		// 		}
		// }
		
		// SECTION("No reset when modifying same bits") {

		// }
		}

	}
};
extern lest::tests & specification();

MODULE( specification(), module )
