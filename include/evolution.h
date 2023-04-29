#ifndef BECPP_EVOLUTION_H
#define BECPP_EVOLUTION_H

#include "data.h"
#include "wavefunction.h"
#include <complex>

constexpr std::complex<double> I{0, 1};

/** Computes the Fourier step of the evolution.
 *
 * Computes the Fourier subsystem of the evolution equations for a 1D system.
 *
 * @param wfn The 1D wavefunction object.
 * @param params Struct containing the parameters of the system.
 */
void fourierStep(Wavefunction1D& wfn, const Parameters& params);

/** Computes the Fourier step of the evolution.
 *
 * Computes the Fourier subsystem of the evolution equations for a 2D system.
 *
 * @param wfn The 2D wavefunction object.
 * @param params Struct containing the parameters of the system.
 */
void fourierStep(Wavefunction2D& wfn, const Parameters& params);

/** Computes the Fourier step of the evolution.
 *
 * Computes the Fourier subsystem of the evolution equations for a 3D system.
 *
 * @param wfn The 3D wavefunction object.
 * @param params Struct containing the parameters of the system.
 */
void fourierStep(Wavefunction3D& wfn, const Parameters& params);

/** Computes the non-linear step of the evolution.
 *
 * Computes the non-linear subsystem of the evolution equations for a 1D system.
 *
 * @param wfn The 1D wavefunction object.
 * @param params Struct containing the parameters of the system.
 */
void interactionStep(Wavefunction1D& wfn, const Parameters& params);

/** Computes the non-linear step of the evolution.
 *
 * Computes the non-linear subsystem of the evolution equations for a 2D system.
 *
 * @param wfn The 2D wavefunction object.
 * @param params Struct containing the parameters of the system.
 */
void interactionStep(Wavefunction2D& wfn, const Parameters& params);

/** Computes the non-linear step of the evolution.
 *
 * Computes the non-linear subsystem of the evolution equations for a 3D system.
 *
 * @param wfn The 3D wavefunction object.
 * @param params Struct containing the parameters of the system.
 */
void interactionStep(Wavefunction3D& wfn, const Parameters& params);

/** Calculates the atom number of the wavefunction.
 *
 * @param wfn The 1D wavefunction object.
 */
double calculateAtomNum(const Wavefunction1D& wfn);

/** Calculates the atom number of the wavefunction.
 *
 * @param wfn The 2D wavefunction object.
 */
double calculateAtomNum(const Wavefunction2D& wfn);

/** Calculates the atom number of the wavefunction.
 *
 * @param wfn The 2D wavefunction object.
 */
double calculateAtomNum(const Wavefunction3D& wfn);

/** Calculates the atom number of the wavefunction.
 *
 * @param wfn The 3D wavefunction object.
 */

/** Renormalises the atom number of the system.
 *
 * Renormalises the atom number of the system if it has been altered, typically
 * when using imaginary time evolution.
 *
 * @param wfn The 1D wavefunction object.
 */
void renormaliseAtomNum(Wavefunction1D& wfn);

/** Renormalises the atom number of the system.
 *
 * Renormalises the atom number of the system if it has been altered, typically
 * when using imaginary time evolution.
 *
 * @param wfn The 2D wavefunction object.
 */
void renormaliseAtomNum(Wavefunction2D& wfn);

/** Renormalises the atom number of the system.
 *
 * Renormalises the atom number of the system if it has been altered, typically
 * when using imaginary time evolution.
 *
 * @param wfn The 3D wavefunction object.
 */
void renormaliseAtomNum(Wavefunction3D& wfn);

#endif  // BECPP_EVOLUTION_H
