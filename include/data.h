#ifndef BECPP_DATA_H
#define BECPP_DATA_H

#include "grid.h"
#include "highfive/H5DataSet.hpp"
#include "highfive/H5DataSpace.hpp"
#include "highfive/H5File.hpp"
#include "wavefunction.h"
#include <string>
#include <vector>

/** Struct containing all the parameters of the system.
 */
struct Parameters {
  double intStrength{};             ///< Interaction strength
  std::vector<double> trap{};       ///< Trapping potential
  int numTimeSteps{};               ///< Number of time steps in the simulation
  std::complex<double> timeStep{};  ///< Time step increment
  double currentTime{};             ///< Current time of the simulation
};

/** DataManager class that handles all the details of the save system of BEC++.
 * It automatically creates the appropriate datasets upon construction of the
 * object, and saves the initial details of the parameters and numerical grid.
 * It also provides a function for saving the current wavefunction data to a
 * file.
 */
class DataManager1D {
 private:
  unsigned int m_saveIndex{0};

  void saveParameters(const Parameters& params, const Grid1D& grid);
  void generateWavefunctionDatasets(const Grid1D& grid);

 public:
  /** Constructs the DataManager object. It automatically saves and creates the
   * datasets for the wave function, numerical grid, and parameters.
   *
   * @param filename The desired name of the file.
   * @param params The struct containing the system parameters.
   * @param grid The 1D grid object of the system.
   */
  DataManager1D(const std::string& filename, const Parameters& params,
                const Grid1D& grid);

  /** Saves the current wave function data to the file.
   *
   * @param wfn The Wavefunction object of the system.
   */
  void saveWavefunctionData(Wavefunction1D& wfn);

  std::string filename;  ///< Filename of the .hdf5 file
  HighFive::File file;   ///< Reference to the underlying .hdf5 file.
};

/** DataManager class that handles all the details of the save system of BEC++.
 * It automatically creates the appropriate datasets upon construction of the
 * object, and saves the initial details of the parameters and numerical grid.
 * It also provides a function for saving the current wavefunction data to a
 * file.
 */
class DataManager2D {
 private:
  unsigned int m_saveIndex{0};

  void saveParameters(const Parameters& params, const Grid2D& grid);
  void generateWavefunctionDatasets(const Grid2D& grid);

 public:
  /** Constructs the DataManager object. It automatically saves and creates the
   * datasets for the wave function, numerical grid, and parameters.
   *
   * @param filename The desired name of the file.
   * @param params The struct containing the system parameters.
   * @param grid The 2D grid object of the system.
   */
  DataManager2D(const std::string& filename, const Parameters& params,
                const Grid2D& grid);

  /** Saves the current wave function data to the file.
   *
   * @param wfn The Wavefunction object of the system.
   */
  void saveWavefunctionData(Wavefunction2D& wfn);

  std::string filename;  ///< Filename of the .hdf5 file

  HighFive::File file;  ///< Reference to the underlying .hdf5 file.
};

/** DataManager class that handles all the details of the save system of BEC++.
 * It automatically creates the appropriate datasets upon construction of the
 * object, and saves the initial details of the parameters and numerical grid.
 * It also provides a function for saving the current wavefunction data to a
 * file.
 */
class DataManager3D {
 private:
  unsigned int m_saveIndex{0};

  void saveParameters(const Parameters& params, const Grid3D& grid);
  void generateWavefunctionDatasets(const Grid3D& grid);

 public:
  /** Constructs the DataManager object. It automatically saves and creates the
   * datasets for the wave function, numerical grid, and parameters.
   *
   * @param filename The desired name of the file.
   * @param params The struct containing the system parameters.
   * @param grid The 3D grid object of the system.
   */
  DataManager3D(const std::string& filename, const Parameters& params,
                const Grid3D& grid);

  /** Saves the current wave function data to the file.
   *
   * @param wfn The Wavefunction object of the system.
   */
  void saveWavefunctionData(Wavefunction3D& wfn);

  std::string filename;  ///< Filename of the .hdf5 file

  HighFive::File file;  ///< Reference to the underlying .hdf5 file.
};

#endif  // BECPP_DATA_H
