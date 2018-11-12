# High-Tune: RenDeRer

This program is a part of the [High-Tune](http://www.umr-cnrm.fr/high-tune/)
project: it illustrates the implementation of efficient radiative transfer
Monte-Carlo algorithms in cloudy atmospheres.

This program implements a rendering algorithm that computes the radiance in the
spectral interval [380, 780] nanometres that reaches an image through a pinhole
camera. The rendered scene is at least composed of an 1D atmosphere along the Z
axis. Optionally, one can add 3D data describing the cloud properties and/or a
geometry describing the ground with a lambertian reflectivity. The clouds and
the ground, can be both infinitely repeated along the X and Y axis.

In addition of shared memory parallelism, htrdr supports the [*M*essage
*P*assing *I*nterface](https://www.mpi-forum.org/) specification to
parallelise its computations in a distribute memory environment; the HTRDR
binary can be run either directly or through a MPI process launcher like
`mpirun`.

## How to build

This program is compatible GNU/Linux 64-bits. It relies on the
[CMake](http://www.cmake.org) and the
[RCMake](https://gitlab.com/vaplv/rcmake/) packages to build.  It also depends
on the
[HTCP](https://gitlab.com/meso-star/htcp/),
[HTGOP](https://gitlab.com/meso-star/htgop/),
[HTMIE](https://gitlab.com/meso-star/htmie/),
[RSys](https://gitlab.com/vaplv/rsys/),
[Star-3D](https://gitlab.com/meso-star/star-3d/),
[Star-3DAW](https://gitlab.com/meso-star/star-3daw/),
[Star-SF](https://gitlab.com/meso-star/star-sf/),
[Star-SP](https://gitlab.com/meso-star/stat-sp/) and
[Star-VX](https://gitlab.com/meso-star/star-vx/) libraries and on
[OpenMP](http://www.openmp.org) 1.2 and the
[MPI](https://www.mpi-forum.org/) 2.0 specification to parallelize its
computations.

First ensure that CMake is installed on your system. Then install the RCMake
package as well as the aforementioned prerequisites. Finally generate the
project from the `cmake/CMakeLists.txt` file by appending to the
`CMAKE_PREFIX_PATH` variable the install directories of its dependencies. The
resulting project can be edited, built, tested and installed as any CMake
project. Refer to the [CMake](https://cmake.org/documentation) for further
informations on CMake.

## License

htrdr is a free software copyright (C) 2018 Centre National de la Recherche
Scientifique, Université Paul Sabatier <contact-edstar@laplace.univ-tlse.fr>,
[|Meso|Star](http://www.meso-star.com) <contact@meso-star.com>. It is released
under the GPL v3+ license: GNU GPL version 3 or later. You are welcome to
redistribute it under certain conditions; refer to the COPYING file for
details.

