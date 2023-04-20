#ifndef BECPP_EVOLUTION_H
#define BECPP_EVOLUTION_H

#include "data.h"
#include "wavefunction.h"
#include <complex>

constexpr std::complex<double> I{0, 1};

void fourierStep(Wavefunction1D& wfn, const Parameters& params);
void fourierStep(Wavefunction2D& wfn, const Parameters& params);
void fourierStep(Wavefunction3D& wfn, const Parameters& params);

void interactionStep(Wavefunction1D& wfn, const Parameters& params);
void interactionStep(Wavefunction2D& wfn, const Parameters& params);
void interactionStep(Wavefunction3D& wfn, const Parameters& params);

double calculateAtomNum(const Wavefunction1D& wfn);
double calculateAtomNum(const Wavefunction2D& wfn);
double calculateAtomNum(const Wavefunction3D& wfn);

void renormaliseAtomNum(Wavefunction1D& wfn);
void renormaliseAtomNum(Wavefunction2D& wfn);
void renormaliseAtomNum(Wavefunction3D& wfn);

#endif//BECPP_EVOLUTION_H
