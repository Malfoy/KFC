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
				EXPECT(sizes[0] == 0u);
				EXPECT(sizes[1] == byte_length*8/2);
				EXPECT(sizes[2] == 0u);
				EXPECT(sizes[3] == byte_length*8/4);
				EXPECT(sizes[4] == 0u);
				EXPECT(sizes[5] == byte_length*8/4);

				cbf.insert(val, 1);
				// Must contain 1s only in the 1st BF
				sizes = cbf.filter_sizes();
				EXPECT(sizes[0] == 2u);
				EXPECT(sizes[2] == 0u);
				EXPECT(sizes[4] == 0u);

				cbf.insert(val, 1);
				//                 in the 1st & 2nd BFs
				sizes = cbf.filter_sizes();
				EXPECT(sizes[0] == 2u);
				EXPECT(sizes[2] == 2u);
				EXPECT(sizes[4] == 0u);

				cbf.insert(val, 1);
				//                 in the 1st & 2nd & 3rd BFs
				sizes = cbf.filter_sizes();
				EXPECT(sizes[0] == 2u);
				EXPECT(sizes[2] == 2u);
				EXPECT(sizes[4] == 2u);

				cbf.insert(val, 1);
				// BFs must not change
				sizes = cbf.filter_sizes();
				EXPECT(sizes[0] == 2u);
				EXPECT(sizes[1] == byte_length*8/2);
				EXPECT(sizes[2] == 2u);
				EXPECT(sizes[3] == byte_length*8/4);
				EXPECT(sizes[4] == 2u);
				EXPECT(sizes[5] == byte_length*8/4);
			}

		}

	}
};
extern lest::tests & specification();

MODULE( specification(), module )
