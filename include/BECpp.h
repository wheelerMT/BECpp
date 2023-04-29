#ifndef BECPP_H
#define BECPP_H

#include "data.h"
#include "evolution.h"
#include "grid.h"
#include "wavefunction.h"

/** \mainpage Welcome to BEC++!
 *
 * BEC++ is a parallelised Gross-Pitaevskii equation solver for use in simulating the
 * dynamics of scalar Bose-Einstein condensate systems.
 * It offers an easy-to-use API that allows users to not worry about often tedious underlying numerical
 * details, and instead focus on the results of the simulations.
 * BEC++ is parallelised using OpenMP if available, making simulation times short.
 *
 * A complete API reference is provided in this documentation, and examples can be found in the examples folder on the
 * GitHub.
 * Here we provide a brief overview on how to get up and running with BEC++.
 *
 * \section getting_started Getting started
 * BEC++ is a static library, and should be built with CMake. The underlying library binary should then be linked
 * to your C++ program.
 *
 * The typical workflow of BEC++ is as follows:
 * 1. Create the numerical grid.
 * 2. Specify the parameters of the system.
 * 3. Construct the initial wave function object, and specify the initial state.
 * 4. If using, set up the DataManager class to handle correct saving of the parameters and wave function data.
 * 5. Inside a loop, call the evolution functions, and, if using a DataManager, save the wave function at your desired intervals.
 *
 * See the examples folder of the main GitHub repo for C++ code using the library.
 */
#endif  // BECPP_H
