<p align="center"><img src="docs/BEC++.png" alt="logo" ></p>

<h4 align="center">A parallelised Gross-Pitaevskii equation solver.</h4>

---

# <u> What is BEC++?</u>
BEC++ is a fast, adaptable Gross-Pitaevskii equation solver for scalar Bose-Einstein condensate systems.
It offers an easy-to-use interface allowing users to quickly get up and running in simulating the dynamics of these systems.

# <u> What can BEC++ do?</u>
BEC++ provides solvers for the Gross-Pitaevskii equation that are parallelised using OpenMP.
It provides support for 1D, 2D and 3D systems, with an intuite API.
It revolves around a `Wavefunction` class that handles the details of the wave function in one, seamless interface.

For examples of what BEC++ can do, see the [examples](examples/) folder.

## What about other BEC systems?
Currently, BEC++ only supports scalar BEC systems.
However, for users interested in simulating the dynamics of other BEC systems such as spinors, see my other library: [PyGPE](https://github.com/wheelerMT/pygpe/).
