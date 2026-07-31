// RandomGen.C

#include "TRandomGen.h"
#include "TRandom.h"
#include "TStopwatch.h"

#include <vector>
#include <iostream>

void RandomGen( uint seed = 0 )
{

    // ===================================================
    // TRandom scalar generation
    // ===================================================

    TRandom3 rng2(seed);
    double randomnumber;
    for (size_t i=0; i<5; ++i) {
        randomnumber = rng2.Rndm();
        std::cout << "Random number: " << i << " : "<< randomnumber << std::endl;
    }

    const size_t N = 10000000;

    std::vector<Double_t> v(N);

    TStopwatch sw_scalar2;
    sw_scalar2.Start();

    for (size_t i = 0; i < N; ++i)
        v[i] = rng2.Rndm();

    sw_scalar2.Stop();

    double checksum_scalar2 = 0.0;
    for (auto x : v)
        checksum_scalar2 += x;

    std::cout << "TRandom3 scalar:            "
              << sw_scalar2.RealTime()
              << " s"
              << "   checksum = "
              << checksum_scalar2
              << std::endl;


    // ===================================================
    // TRandom bulk generation
    // ===================================================

    TStopwatch sw_array2;
    sw_array2.Start();

    rng2.RndmArray(N, v.data());

    sw_array2.Stop();

    double checksum_array2 = 0.0;
    for (auto x : v)
        checksum_array2 += x;

    std::cout << "TRandom3 RndmArray:         "
              << sw_array2.RealTime()
              << " s"
              << "   checksum = "
              << checksum_array2
              << std::endl;


    // ===================================================
    // Speedup
    // ===================================================

    std::cout << "Speedup (TRandom3): "
              << sw_scalar2.RealTime() / sw_array2.RealTime()
              << "x" << std::endl;
}