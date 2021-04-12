#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <math.h>
namespace constants
{
    constexpr double avo{6.022e23};
    constexpr double pi{4.0 * atan(1.0)};
    constexpr double hbar{1.054571800e-34}; // Js
    constexpr double c0{299792458.0}; // ms
    constexpr double echarge{1.6021766028e-19}; //C
    constexpr double emass{9.10938356e-31}; // kg
    constexpr double epsilon0{8.854187817e-12}; // Fm^-1
    constexpr double bohr{4 * constants::pi * constants::epsilon0 * pow(constants::hbar,2) /
                          (constants::emass * pow(constants::echarge,2) ) };
    constexpr double bohrAng{bohr * 1e10}; // for convenience in angstroms
    constexpr double gravity{9.812}; // ms^-2i
    constexpr double hartree{pow(hbar, 2) / (emass * pow(bohr, 2))};
    constexpr double boltz{1.386485279e-23}; // J K^-1
    constexpr double boltzHar{boltz / hartree}; // Eh K^-1
}
#endif
