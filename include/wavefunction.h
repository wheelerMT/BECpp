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

class Wavefunction1D {
 private:
  const Grid1D &m_grid;
  FFTPlans m_plans{};
  complexVector_t m_component{};
  complexVector_t m_fourierComponent{};
  double m_atomNumber{};

  void createFFTPlans(const Grid1D &grid);
  void destroyFFTPlans() const;
  void updateAtomNumber();

 public:
  explicit Wavefunction1D(const Grid1D &grid);
  ~Wavefunction1D();

  [[nodiscard]] const Grid1D &grid() const;
  [[nodiscard]] complexVector_t &component();
  [[nodiscard]] complexVector_t &fourierComponent();
  [[nodiscard]] std::vector<double> density() const;
  [[nodiscard]] double atomNumber() const;

  void fft() const;
  void ifft();
  void setComponent(complexVector_t &component);
};

class Wavefunction2D {
 private:
  const Grid2D &m_grid;
  FFTPlans m_plans{};
  complexVector_t m_component{};
  complexVector_t m_fourierComponent{};
  double m_atomNumber{};

  void createFFTPlans(const Grid2D &grid);
  void destroyFFTPlans() const;
  void updateAtomNumber();

 public:
  explicit Wavefunction2D(const Grid2D &grid);
  ~Wavefunction2D();

  [[nodiscard]] const Grid2D &grid() const;
  [[nodiscard]] complexVector_t &component();
  [[nodiscard]] complexVector_t &fourierComponent();
  [[nodiscard]] std::vector<double> density() const;
  [[nodiscard]] double atomNumber() const;

  void fft() const;
  void ifft();
  void setComponent(complexVector_t &component);
};

class Wavefunction3D {
 private:
  const Grid3D &m_grid;
  FFTPlans m_plans{};
  complexVector_t m_component{};
  complexVector_t m_fourierComponent{};
  double m_atomNumber{};

  void createFFTPlans(const Grid3D &grid);
  void destroyFFTPlans() const;
  void updateAtomNumber();

 public:
  explicit Wavefunction3D(const Grid3D &grid);
  ~Wavefunction3D();

  [[nodiscard]] const Grid3D &grid() const;
  [[nodiscard]] complexVector_t &component();
  [[nodiscard]] complexVector_t &fourierComponent();
  [[nodiscard]] std::vector<double> density() const;
  [[nodiscard]] double atomNumber() const;

  void fft() const;
  void ifft();
  void setComponent(complexVector_t &component);
};

#endif  // BECPP_WAVEFUNCTION_H
