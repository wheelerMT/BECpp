#ifndef BECPP_GRID_H
#define BECPP_GRID_H

#include <tuple>
#include <utility>
#include <vector>

/**
 * Abstract Grid class, used to provide a framework for the different
 * dimensionality Grid classes.
 */
class Grid {
 protected:
  virtual void constructGridParams() = 0;
  virtual void constructMesh() = 0;
  virtual ~Grid() = default;
};

/** Struct containing the 1D grid meshes.
 */
struct Mesh1D {
  std::vector<double> xMesh{};
  std::vector<double> xFourierMesh{};
  std::vector<double> wavenumber{};
};

/** The 1D Grid class of the system.
 *
 * This classs encapsulates all of the details of the numerical grid,
 * from the grid points to grid spacings etc. It also contains useful functions
 * for operating and retrieving the underlying grid data.
 */
class Grid1D : public Grid {
 private:
  void constructGridParams() override;
  void constructMesh() override;

  const unsigned int m_gridPoints{};
  double m_gridSpacing{};
  double m_fourierGridSpacing{};
  double m_length{};
  Mesh1D m_mesh{};

 public:
  Grid1D(unsigned int nx, double dx);
  ~Grid1D() override = default;

  /** Returns the shape of the grid.
   */
  [[nodiscard]] unsigned int shape() const;

  /** Returns the grid spacing of the position space numerical grid.
   */
  [[nodiscard]] double gridSpacing() const;

  /** Returns the grid spacing of the Fourier space numerical grid.
   */
  [[nodiscard]] double fourierGridSpacing() const;

  /** Returns the length of the position space numerical grid. This is
   * calculated as the number of grid points x grid spacing.
   */
  [[nodiscard]] double gridLength() const;

  /** Returns a reference to the xMesh of the numerical grid.
   */
  [[nodiscard]] std::vector<double>& xMesh();

  /** Returns a reference to the wavenumber mesh of the numerical grid.
   */
  [[nodiscard]] std::vector<double>& wavenumber();

  friend class Wavefunction;
};

/** Struct containing the 2D grid meshes.
 */
struct Mesh2D {
  std::vector<double> xMesh{};
  std::vector<double> yMesh{};
  std::vector<double> xFourierMesh{};
  std::vector<double> yFourierMesh{};
  std::vector<double> wavenumber{};
};

/** The 2D Grid class of the system.
 *
 * This classs encapsulates all of the details of the numerical grid,
 * from the grid points to grid spacings etc. It also contains useful functions
 * for operating and retrieving the underlying grid data.
 */
class Grid2D : public Grid {
 private:
  void constructGridParams() override;
  void constructMesh() override;

  const std::tuple<unsigned int, unsigned int> m_gridPoints{};
  const std::tuple<double, double> m_gridSpacing{};
  std::tuple<double, double> m_fourierGridSpacing{};
  std::tuple<double, double> m_gridLength{};
  Mesh2D m_mesh{};

 public:
  Grid2D(std::tuple<unsigned int, unsigned int> points,
         std::tuple<double, double> gridSpacing);
  ~Grid2D() override = default;

  /** Returns the shape of the grid as a tuple. The first and second entries of
   * the tuple correspond to xPoints and yPoints, respectively.
   */
  [[nodiscard]] std::tuple<unsigned int, unsigned int> shape() const;

  /** Returns the position space grid spacings as a tuple. The first and second
   * entries of the tuple correspond to xGridSpacing and yGridSpacing,
   * respectively.
   */
  [[nodiscard]] std::tuple<double, double> gridSpacing() const;

  /** Returns the Fourier space grid spacings as a tuple. The first and second
   * entries of the tuple correspond to xFourierGridSpacing and
   * yFourierGridSpacing, respectively.
   */
  [[nodiscard]] std::tuple<double, double> fourierGridSpacing() const;

  /** Returns the lengths of the position space numerical grid as a tuple. The
   * first and second entries of the tuple correspond to xLength and yLength,
   * respectively.
   */
  [[nodiscard]] std::tuple<double, double> gridLength() const;

  /** Returns a reference to the xMesh of the numerical grid.
   */
  [[nodiscard]] std::vector<double>& xMesh();

  /** Returns a reference to the yMesh of the numerical grid.
   */
  [[nodiscard]] std::vector<double>& yMesh();

  /** Returns a reference to the wavenumber mesh of the numerical grid.
   */
  [[nodiscard]] std::vector<double>& wavenumber();

  friend class Wavefunction;
};

/** Struct containing the 3D grid meshes.
 */
struct Mesh3D {
  std::vector<double> xMesh{};
  std::vector<double> yMesh{};
  std::vector<double> zMesh{};
  std::vector<double> xFourierMesh{};
  std::vector<double> yFourierMesh{};
  std::vector<double> zFourierMesh{};
  std::vector<double> wavenumber{};
};

/** The 3D Grid class of the system.
 *
 * This classs encapsulates all of the details of the numerical grid,
 * from the grid points to grid spacings etc. It also contains useful functions
 * for operating and retrieving the underlying grid data.
 */
class Grid3D : public Grid {
 private:
  void constructGridParams() override;
  void constructMesh() override;

  const std::tuple<unsigned int, unsigned int, unsigned int> m_gridPoints{};
  const std::tuple<double, double, double> m_gridSpacing{};
  std::tuple<double, double, double> m_fourierGridSpacing{};
  std::tuple<double, double, double> m_gridLength{};
  Mesh3D m_mesh{};

 public:
  Grid3D(std::tuple<unsigned int, unsigned int, unsigned int> points,
         std::tuple<double, double, double> gridSpacing);
  ~Grid3D() override = default;

  /** Returns the shape of the grid as a tuple. The first, second and third
   * entries of the tuple correspond to xPoints, yPoints and zPoints,
   * respectively.
   */
  [[nodiscard]] std::tuple<unsigned int, unsigned int, unsigned int>
  shape() const;

  /** Returns the position space grid spacings as a tuple. The first, second and
   * third entries of the tuple correspond to xGridSpacing, yGridSpacing and
   * zGridSpacing, respectively.
   */
  [[nodiscard]] std::tuple<double, double, double> gridSpacing() const;

  /** Returns the Fourier space grid spacings as a tuple. The first, second and
   * third entries of the tuple correspond to xFourierGridSpacing,
   * yFourierGridSpacing and zFourierGridSpacing, respectively.
   */
  [[nodiscard]] std::tuple<double, double, double> fourierGridSpacing() const;

  /** Returns the lengths of the position space grids as a tuple. The first,
   * second and third entries of the tuple correspond to xLength, yLength and
   * zLength, respectively.
   */
  [[nodiscard]] std::tuple<double, double, double> gridLength() const;

  /** Returns a reference to the xMesh of the numerical grid.
   */
  [[nodiscard]] std::vector<double>& xMesh();

  /** Returns a reference to the yMesh of the numerical grid.
   */
  [[nodiscard]] std::vector<double>& yMesh();

  /** Returns a reference to the zMesh of the numerical grid.
   */
  [[nodiscard]] std::vector<double>& zMesh();

  /** Returns a reference to the wavenumber mesh of the numerical grid.
   */
  [[nodiscard]] std::vector<double>& wavenumber();

  friend class Wavefunction;
};

#endif  // BECPP_GRID_H
