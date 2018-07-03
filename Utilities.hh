
#ifndef UTILITIES
#define UTILITIES 1

#include <cmath>
#include <iostream>
#include <fstream>

namespace Utilities
{
    double fwhm2sigma = 1.0 / ( 2.0*std::sqrt(2.0*std::log(2.)) );

    bool Exists(const char * filename) {
        std::ifstream infile(filename);
        return infile.good();
    }
}

#endif
