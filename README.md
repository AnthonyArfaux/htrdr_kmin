# High-Tune: RenDeRer

This program is a part of the [High-Tune](http://www.umr-cnrm.fr/high-tune/)
project: it illustrates the implementation of efficient radiative transfer
Monte-Carlo algorithms in cloudy atmospheres.

htrdr is an image renderer in the visible part of the spectrum, for scenes
composed of an atmospheric gas mixture, clouds, and a ground. It uses spectral
data that should be provided for the pressure and temperature atmospheric
vertical profile defined along the Z axis, the liquid water content in
suspension within the clouds that is a result of Large Eddy Simulation
computations, and the optical properties of water droplets that have been
obtained from a Mie code. The user also has to provide: the characteristics of
the simulated camera, the sensor definition, and the position of the sun. It is
also possible to provide a geometry representing the ground. Both, the clouds
and the ground, can be infinitely repeated along the X and Y axis.

htrdr evaluates the intensity incoming on each pixel of the sensor array. The
underlying algorithm is based on a Monte-Carlo method: it consists in
simulating a given number of optical paths originating from the camera,
directed into the atmosphere, taking into account light absorption and
scattering phenomena. The computation is performed over the whole visible part
of the spectrum, for the three components of the CIE 1931 XYZ colorimetric
space that are subsequently recombined in order to obtain the final color for
each pixel, and finally the whole image of the scene as seen from the required
observation position.

In addition of shared memory parallelism, htrdr supports the [*M*essage
*P*assing *I*nterface](https://www.mpi-forum.org/) specification to
parallelise its computations in a distribute memory environment; the htrdr
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
[MPI](https://www.mpi-forum.org/) specification to parallelize its
computations.

First ensure that CMake is installed on your system. Then install the RCMake
package as well as the aforementioned prerequisites. Finally generate the
project from the `cmake/CMakeLists.txt` file by appending to the
`CMAKE_PREFIX_PATH` variable the install directories of its dependencies. The
resulting project can be edited, built, tested and installed as any CMake
project. Refer to the [CMake](https://cmake.org/documentation) for further
informations on CMake.

## Release notes

### Version 0.0.4

- Fix the computation of the surface scattering: there was a bug in how Russian
  roulette was implemented at surface scattering leading to an underestimation
  of the surface reflection.
- Update the thread allocation policy: by default, the number of threads is now
  defined as the maximum between the number of processors detected by OpenMP
  and the number of threads defined by the `OMP_NUM_THREADS` environment
  variable. This variable can be used to counteract the number of processors
  detected by OpenMP that can be lower than the real number of processors of
  the system.

### Version 0.0.3

- Fix compilation on systems with a GNU C Library whose version is less than
  2.19.
- Fix a possible invalid memory access to cloud data leading to segmentation
  faults.

## License

Copyright (C) 2018-2019 Centre National de la Recherche
Scientifique (CNRS), [|Meso|Star](http://www.meso-star.com)
<contact@meso-star.com>, Université Paul Sabatier
<contact-edstar@laplace.univ-tlse.fr>. htrdr is free software released under
the GPL v3+ license: GNU GPL version 3 or later. You are welcome to
redistribute it under certain conditions; refer to the COPYING file for
details.

