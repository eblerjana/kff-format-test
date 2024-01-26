#include <iostream>
#include "kff_io.hpp"
#include <cassert>
#include <sstream>  

using namespace std;

uint8_t uint8_packing(std::string sequence);
/* Encode the sequence into an array of uint8_t packed sequence slices.
 * The encoded sequences are organised in big endian order.
 */
void encode_sequence(std::string sequence, uint8_t * encoded) {
        size_t size = sequence.length();
        // Encode the truncated first 8 bits sequence
        size_t remnant = size % 4;
        if (remnant > 0) {
                encoded[0] = uint8_packing(sequence.substr(0, remnant));
                encoded += 1;
        }

        // Encode all the 8 bits packed
        size_t nb_uint_needed = size / 4;
        for (size_t i=0 ; i<nb_uint_needed ; i++) {
                encoded[i] = uint8_packing(sequence.substr(remnant + 4*i, 4));
                // encoded[i] = uint8_packing(sequence + remnant + (i<<2), 4);
        }
}

/* Transform a char * sequence into a uint8_t 2-bits/nucl
 * Encoding ACTG
 * Size must be <= 4
 */
uint8_t uint8_packing(std::string sequence) {
        size_t size = sequence.length();
        assert(size <= 4);

        uint8_t val = 0;
        for (size_t i=0 ; i<size ; i++) {
                val <<= 2;
                val += (sequence[i] >> 1) & 0b11;
        }

        return val;
}


void uint8_unpacking(uint8_t packed, char * decoded, size_t size);
string decode_sequence(uint8_t * encoded, size_t size) {
        stringstream ss;
        char tmp_chars[4] = {0, 0, 0, 0};

        // Decode the truncated first compacted 8 bits
        size_t remnant = size % 4;
        if (remnant > 0) {
                uint8_unpacking(encoded[0], tmp_chars, remnant);
                for (size_t i=0 ; i<remnant ; i++) {
                        ss << tmp_chars[i];
                }
                encoded += 1;
        }

        // Decode all the 8 bits packed
        size_t nb_uint_used = size / 4;
        for (size_t i=0 ; i<nb_uint_used ; i++) {
                uint8_unpacking(encoded[i], tmp_chars, 4);
                for (size_t i=0 ; i<4 ; i++) {
                        ss << tmp_chars[i];
                }
        }

        return ss.str();
}

char const_nucleotides[4] = {'A', 'C', 'T', 'G'};
void uint8_unpacking(uint8_t packed, char * decoded, size_t size) {
        assert(size <= 4);

        size_t offset = 4 - size;
        for (size_t i=0 ; i<size ; i++) {
                decoded[i] = const_nucleotides[(packed >> ((3-i-offset) * 2)) & 0b11];
        }
}




int main(int argc, char* argv[]) {

        Kff_reader * reader = new Kff_reader("test.kff");
	uint64_t k = reader->k;
        uint64_t i = 0;
        while (reader->has_next()) {
                uint8_t * kmer;
                uint8_t * data;
                reader->next_kmer(kmer, data);
                cout << (i++) << " " << decode_sequence(kmer, k) << " " << (uint)*data << endl;
        }

        delete reader;
	return 0;
};
