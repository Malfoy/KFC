#include "lest.hpp"
#include "BitSet.hpp"

// clang-format off
using namespace std;
const lest::test module[] = {

	CASE( "A BitSet can be built and queried" "[bitset]" ) {

		SETUP("A BitSet with 100 bits") {
			const unsigned length = 100;
			BitSet bits(length);

			EXPECT(bits.size() == length);

			SECTION("At initialization all bits are zero") {
				for (unsigned i = 0; i < length; i++) {
					EXPECT(! bits.get(i));
				}
			}

			SECTION("Properly (un)sets bits") {
				bits.set(0);
				EXPECT(bits.get(0));

				bits.set(15);
				EXPECT(bits.get(15));

				bits.set(BitSet::word_bitwidth-1);
				EXPECT(bits.get(BitSet::word_bitwidth-1));

				bits.set(99);
				EXPECT(bits.get(99));

				for (unsigned i = 0; i < length; i++) {
					if (i != 0 && i != 15 && i != (BitSet::word_bitwidth - 1) && i != 99) {
						EXPECT(! bits.get(i));
					}
				}

				bits.reset();
				for (unsigned i = 0; i < length; i++) {
					EXPECT(! bits.get(i));
				}
			}

			SECTION("Setting bits several times") {
				for (unsigned i = 0 ; i < length; i++) {
					bits.set(i);
					EXPECT(bits.get(i));
					bits.set(i);
					EXPECT(bits.get(i));
				}
				for (unsigned i = 0 ; i < length; i++) {
					EXPECT(bits.get(i));
				}
			}
		}
	}
};

extern lest::tests & specification();

MODULE( specification(), module )
