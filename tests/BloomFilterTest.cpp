#include "lest.hpp"
#include "BloomFilter.hpp"
#include <array>

const lest::test module[] = {
  CASE("Setting and querying BloomFilter") {
    SETUP("A BloomFilter with 100 bits") {
      const uint length = 100;
      BloomFilter bf(length, 1, .5);
      EXPECT(bf.size() == length);

      
      SECTION("Setting and querying") {
        uint8_t val[] = {0};
        bf.add(val, 1);
        EXPECT(bf.possiblyContains(val, 1));
        EXPECT(bf.nbBitsSet() == (uint)1);

        val[0] = 1;
        EXPECT(! bf.possiblyContains(val, 1));

        val[0] = 255;
        bf.add(val, 1);
        EXPECT(bf.possiblyContains(val, 1));
        // This will depend on hash functions but the test is likely to succeed
        EXPECT(bf.nbBitsSet() == (uint)2);
        val[0] = 0;
        EXPECT(bf.possiblyContains(val, 1));
        val[0] = 0;
        EXPECT(bf.possiblyContains(val, 1));
      }

    SECTION("Reset when setting 50% different bits") {
        uint8_t val[] = {0};
        uint i = 0;
        while (bf.nbBitsSet() < (uint) length/2-1) {
          val[0] = i++;
          bf.add(val, 1);
        }
        uint64_t current_nb_bits_set = bf.nbBitsSet();
        // Next inserted bit should induce a reset
        // But next inserted bit may fall on an existing set bit (proba: 1/2) .
        while (bf.nbBitsSet() == current_nb_bits_set) {
          val[0] = i++;
          bf.add(val, 1);
        }

        // Now the BF should have been reset
        EXPECT(bf.nbBitsSet() == (uint)0);
        for (uint i = 0; i < 255; i++) {
          val[0] = i;
          EXPECT(! bf.possiblyContains(val, 1));
        }
    }
    
    // SECTION("No reset when modifying same bits") {

    // }
    }

  }
};
extern lest::tests & specification();

MODULE( specification(), module )
