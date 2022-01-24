### <img src="https://github.com/wheelerMT/BECpp/blob/master/docs/vids/becpp.gif">

---

# <u> What is BEC++?</u>
BEC++ is aiming to be a fast, adaptable Gross-Pitaevskii equation solver for Bose-Einstein condensate systems.
The project is very much a work in progress. 

# <u> What can BEC++ do?</u>
BEC++ provides solvers for the Gross-Pitaevskii equation that are solved using either multi-threaded CPU with OpenMP,
or GPU-accelerated using CUDA.

I plan to add support for
- Scalar BECs;
- Spin-1 BECs;
- Spin-2 BECs;
- Two-component BECs.

Initially this will be a 2D codebase, eventually expanding to support 1D and 3D grids.
The codebase will also contain a myriad of useful functions for generating initial states and phase profiles that
support vortices.

# <u> Current development </u>
I am currently developing the 1D version of the codebase.
