#include "murmurhash32.hpp"
#include "murmurhash64.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <seqan/seq_io.h>
#include <fstream>
#include <cstdint>
#include <getopt.h>
using namespace std;

uint64_t canonical_kmer(uint64_t kmer, uint k = 31)
{
    uint64_t reverse = 0;
    uint64_t b_kmer = kmer;

    kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
    kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
    kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
    kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
    kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
    reverse = (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - (k << 1));

    return (b_kmer < reverse) ? b_kmer : reverse;
}

class HyperLogLog {
private:
    int b;                // número de bits para índice
    int m;                // número de registros = 2^b
    vector<uint8_t> M;    // contadores
    double alpha_m;       // constante de corrección

    static double get_alpha(int m) {
        if (m == 16) return 0.673;
        else if (m == 32) return 0.697;
        else if (m == 64) return 0.709;
        else return 0.7213 / (1.0 + 1.079 / m);
    }

public:
    explicit HyperLogLog(int b_bits = 14) : b(b_bits), m(1 << b_bits), M(m, 0) {
        alpha_m = get_alpha(m);
    }

    

    void add(uint64_t x) {
    uint64_t hash = murmurhash64A(&x, sizeof(x), 0);

    int idx = hash >> (64 - b);

    // w = los bits restantes 
    uint64_t w = hash << b;

    int rho;
    if (w == 0) {
        rho = (64 - b) + 1;  
    } else {
        rho = __builtin_clzll(w) + 1;
    }

    int max_rho = (64 - b) + 1;
    if (rho > max_rho) rho = max_rho;

    if (rho > M[idx]) M[idx] = (uint8_t)rho;
}


    double estimate() const {
        double Z = 0.0;
        for (auto v : M)
            Z += pow(2.0, -double(v));

        Z = 1.0 / Z;

        double E = alpha_m * m * m * Z;

        if (E <= 2.5 * m) {
            int V = count(M.begin(), M.end(), 0);
            if (V > 0)
                E = m * log((double)m / V);
        }

        return E;
    }
    void merge_from(const HyperLogLog &other) {
        for (size_t i = 0; i < M.size(); i++) {
            M[i] = std::max(M[i], other.M[i]);
        }
    }

};

//  Procesamiento de los kmers
void process_kmers_hll(HyperLogLog &hll, const string &filename, uint k) {
    seqan::SeqFileIn seqFileIn;
    if (!open(seqFileIn, filename.c_str())) {
        cerr << "ERROR: Could not open file " << filename << endl;
        return;
    }

    seqan::CharString id;
    seqan::IupacString seq;

    while (!atEnd(seqFileIn)) {
        try {
            seqan::readRecord(id, seq, seqFileIn);
        } catch (seqan::ParseError &) {
            break;
        }

        uint64_t kmer = 0;
        uint bases = 0;

        for (size_t i = 0; i < length(seq); i++) {
            uint8_t two_bit = 0;

            switch (char(seq[i])) {
                case 'A': case 'a': two_bit = 0; break;
                case 'C': case 'c': two_bit = 1; break;
                case 'G': case 'g': two_bit = 2; break;
                case 'T': case 't': two_bit = 3; break;
                default: 
                    kmer = 0; bases = 0;
                    continue;
            }

            kmer = ((kmer << 2) | two_bit) & ((1ULL << (k*2)) - 1);
            bases++;

            if (bases == k) {
                uint64_t canonical = canonical_kmer(kmer, k);
                hll.add(canonical);
                bases--;
            }
        }
    }

    close(seqFileIn);
}

//  Procesamiento de los kmers con minimizers
void process_kmers_minimizers(HyperLogLog &hll, const string &filename, uint k, uint w) {
    seqan::SeqFileIn seqFileIn;
    if (!open(seqFileIn, filename.c_str())) {
        cerr << "ERROR: Could not open file " << filename << endl;
        return;
    }

    seqan::CharString id;
    seqan::IupacString seq;

    while (!atEnd(seqFileIn)) {
        try {
            seqan::readRecord(id, seq, seqFileIn);
        } catch (seqan::ParseError &) {
            break;
        }

        vector<uint64_t> kmer_hashes;
        uint64_t kmer = 0;
        uint bases = 0;

        for (size_t i = 0; i < length(seq); i++) {
            uint8_t two_bit = 0;
            switch (char(seq[i])) {
                case 'A': case 'a': two_bit = 0; break;
                case 'C': case 'c': two_bit = 1; break;
                case 'G': case 'g': two_bit = 2; break;
                case 'T': case 't': two_bit = 3; break;
                default: kmer = 0; bases = 0; kmer_hashes.clear(); continue;
            }

            kmer = ((kmer << 2) | two_bit) & ((1ULL << (k*2)) - 1);
            bases++;

            if (bases == k) {
                uint64_t canonical = canonical_kmer(kmer, k);
                uint64_t hash = murmurhash64A(&canonical, sizeof(canonical), 0);
                kmer_hashes.push_back(hash);

                if (kmer_hashes.size() == w) {
                    uint64_t min_hash = *min_element(kmer_hashes.begin(), kmer_hashes.end());
                    hll.add(min_hash);
                    kmer_hashes.erase(kmer_hashes.begin()); // deslizar ventana
                }

                bases--;
            }
        }
    }

    close(seqFileIn);
}

