#ifndef BECPP_WAVEFUNCTION_H
#define BECPP_WAVEFUNCTION_H

#include "constants.h"
#include "fftw3.h"
#include "grid.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using complexVector_t = std::vector<std::complex<double>>;

struct FFTPlans {
  fftw_plan plan_forward{};
  fftw_plan plan_backward{};
};

/** 1D wave function class.
 *
 * This class encapsulates all the details of the wave function of the system.
 * It handles both the position space and Fourier space arrays, as well as
 * associated functions to operate on those arrays. It is the fundamental object
 * of the library.
 */
class Wavefunction1D {
 private:
  const Grid1D& m_grid;
  FFTPlans m_plans{};
  complexVector_t m_component{};
  complexVector_t m_fourierComponent{};
  double m_atomNumber{};

  void createFFTPlans(const Grid1D& grid);
  void destroyFFTPlans() const;
  void updateAtomNumber();

 public:
  /** Constructs the wave function object from the associated grid object.
   *
   * @param grid The 1D grid object of the system.
   */
  explicit Wavefunction1D(const Grid1D& grid);

  /** Destructor that automatically cleans up stored FFT plans.
   *
   */
  ~Wavefunction1D();

  /** Returns a reference to the grid object of the system.
   *
   * This is particularly useful when we need access to the underlying grid
   * variables, such as grid points, grid spacing etc.
   */
  [[nodiscard]] const Grid1D& grid() const;

  /** Returns a reference to the position space vector.
   */
  [[nodiscard]] complexVector_t& component();

  /** Returns a reference to the Fourier space vector.
   */
  [[nodiscard]] complexVector_t& fourierComponent();

  /** Returns a vector of the density of the system.
   */
  [[nodiscard]] std::vector<double> density() const;

  /** Calculates and returns the current atom number of the system.
   */
  [[nodiscard]] double atomNumber() const;

  /** Performs a forward FFT.
   *
   * Uses the pre-built FFT plans to compute the forward fast fourier transform.
   * In particular, the Fourier space vector will updated with new values based
   * on the position space vector.
   */
  void fft() const;

  /** Performs an inverse FFT.
   *
   * Uses the pre-built FFT plans to compute the inverse fast fourier transform.
   * In particular, the position space vector will updated with new values based
   * on the Fourier space vector.
   */
  void ifft();

  /** Sets the position space vector to the inputted vector.
   *
   * Takes in any complexVector_t array and sets it to the position space
   * vector of the Wavefunction object.
   * Note: the size of the input vector must be the same as the size of the
   * numerical grid.
   *
   * @param component A complexVector_t containing the wave function state.
   */
  void setComponent(complexVector_t& component);
};

/** 2D wave function class.
 *
 * This class encapsulates all the details of the wave function of the system.
 * It handles both the position space and Fourier space arrays, as well as
 * associated functions to operate on those arrays. It is the fundamental object
 * of the library.
 */

class Wavefunction2D {
 private:
  const Grid2D& m_grid;
  FFTPlans m_plans{};
  complexVector_t m_component{};
  complexVector_t m_fourierComponent{};
  double m_atomNumber{};

  void createFFTPlans(const Grid2D& grid);
  void destroyFFTPlans() const;
  void updateAtomNumber();

 public:
  /** Constructs the wave function object from the associated grid object.
   *
   * @param grid The 2D grid object of the system.
   */
  explicit Wavefunction2D(const Grid2D& grid);

  /** Destructor that automatically cleans up stored FFT plans.
   *
   */
  ~Wavefunction2D();

  /** Returns a reference to the grid object of the system.
   *
   * This is particularly useful when we need access to the underlying grid
   * variables, such as grid points, grid spacing etc.
   */
  [[nodiscard]] const Grid2D& grid() const;

  /** Returns a reference to the position space vector. Note that the return
   * type is still a 1D vector, but should be treated as a 2D vector.
   */
  [[nodiscard]] complexVector_t& component();

  /** Returns a reference to the Fourier space vector. Note that the return type
   * is still a 1D vector, but should be treated as a 2D vector.
   */
  [[nodiscard]] complexVector_t& fourierComponent();

  /** Returns a vector of the density of the system. Note that the return type
   * is still a 1D vector, but should be treated as a 2D vector.
   */
  [[nodiscard]] std::vector<double> density() const;

  /** Calculates and returns the current atom number of the system.
   */
  [[nodiscard]] double atomNumber() const;

  /** Performs a forward FFT.
   *
   * Uses the pre-built FFT plans to compute the forward fast fourier transform.
   * In particular, the Fourier space vector will updated with new values based
   * on the position space vector.
   */
  void fft() const;

  /** Performs an inverse FFT.
   *
   * Uses the pre-built FFT plans to compute the inverse fast fourier transform.
   * In particular, the position space vector will updated with new values based
   * on the Fourier space vector.
   */
  void ifft();

  /** Sets the position space vector to the inputted vector.
   *
   * Takes in any complexVector_t array and sets it to the position space
   * vector of the Wavefunction object.
   * Note: the size of the input vector must be the same as the size of the
   * numerical grid.
   *
   * @param component A complexVector_t containing the wave function state.
   */
  void setComponent(complexVector_t& component);
};

/** 3D wave function class.
 *
 * This class encapsulates all the details of the wave function of the system.
 * It handles both the position space and Fourier space arrays, as well as
 * associated functions to operate on those arrays. It is the fundamental object
 * of the library.
 */
class Wavefunction3D {
 private:
  const Grid3D& m_grid;
  FFTPlans m_plans{};
  complexVector_t m_component{};
  complexVector_t m_fourierComponent{};
  double m_atomNumber{};

  void createFFTPlans(const Grid3D& grid);
  void destroyFFTPlans() const;
  void updateAtomNumber();

 public:
  /** Constructs the wave function object from the associated grid object.
   *
   * @param grid The 3D grid object of the system.
   */
  explicit Wavefunction3D(const Grid3D& grid);

  /** Destructor that automatically cleans up stored FFT plans.
   *
   */
  ~Wavefunction3D();

  /** Returns a reference to the grid object of the system.
   *
   * This is particularly useful when we need access to the underlying grid
   * variables, such as grid points, grid spacing etc.
   */
  [[nodiscard]] const Grid3D& grid() const;

  /** Returns a reference to the position space vector. Note that the return
   * type is still a 1D vector, but should be treated as a 3D vector.
   */
  [[nodiscard]] complexVector_t& component();

  /** Returns a reference to the Fourier space vector. Note that the return
   * type is still a 1D vector, but should be treated as a 3D vector.
   */
  [[nodiscard]] complexVector_t& fourierComponent();

  /** Returns a vector of the density of the system. Note that the return type
   * is still a 1D vector, but should be treated as a 3D vector.
   */
  [[nodiscard]] std::vector<double> density() const;

  /** Calculates and returns the current atom number of the system.
   */
  [[nodiscard]] double atomNumber() const;

  /** Performs a forward FFT.
   *
   * Uses the pre-built FFT plans to compute the forward fast fourier transform.
   * In particular, the Fourier space vector will updated with new values based
   * on the position space vector.
   */
  void fft() const;

  /** Performs an inverse FFT.
   *
   * Uses the pre-built FFT plans to compute the inverse fast fourier transform.
   * In particular, the position space vector will updated with new values based
   * on the Fourier space vector.
   */
  void ifft();

  /** Sets the position space vector to the inputted vector.
   *
   * Takes in any complexVector_t array and sets it to the position space
   * vector of the Wavefunction object.
   * Note: the size of the input vector must be the same as the size of the
   * numerical grid.
   *
   * @param component A complexVector_t containing the wave function state.
   */
  void setComponent(complexVector_t& component);
};

#endif  // BECPP_WAVEFUNCTION_H
