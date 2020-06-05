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
[RCMake](https://gitlab.com/vaplv/rcmake/) packages to build. It also depends
on the
[AW](https://gitlab.com/vaplv/loader_aw/#tab-readme),
[HTSky](https://gitlab.com/meso-star/htsky/),
[MruMtl](https://gitlab.com/meso-star/mrumtl/),
[RSys](https://gitlab.com/vaplv/rsys/),
[Star-3D](https://gitlab.com/meso-star/star-3d/),
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

### Version 0.5

#### New feature

Add support of shortwave integration with respect to the Planck function for a
reference temperature whose default value is the blackbody temperature of the
sun. Actually this is the counterpart of the longwave integration introduced by
the "infrared rendering" in the 0.4 version. The main difference is that the
source of radiation is the sun rather than the medium and its boundaries.

The option `-l` that enabled the infrared rendering is now replaced by the new
`-s` option that controls the spectral integration that can be CIE XYZ (i.e.
regular image rendering), longwave or shortwave.

#### Fixes

- Fix the returned sun radiance: the precomputed per spectral band solar
  incoming flux is removed and the sun radiance is now retrieved by directly
  evaluating the monochromatic Planck for the blackbody temperature of the sun.
- Fix CIE XYZ spectral integration: the pdf used to sample the CIE tristimulus
  values was not correctly handled in the Monte-Carlo weight.
- Fix the longwave spectral integration: the Monte-Carlo weight was wrong
  leading to overestimated temperatures.

### Version 0.4

#### New features

- Add support of infrared rendering: when defined, the new `-l` option setups
  the range of longwave into which the rendering is performed. In infrared
  rendering, each pixel stores the radiance per pixel and its associated
  brightness temperature. Spectral integration is done with respect to the
  Planck function for a reference temperature of 290 K.
- The ground geometry can now have several materials whose data vary over the
  spectrum. These materials are listed in a new
  [htrdr-materials](https://gitlab.com/meso-star/htrdr/-/blob/master/doc/htrdr-materials.5.txt)
  file where each materials is defined by a name and a file storing its spectral
  data with respect to the
  [MruMtl](https://gitlab.com/meso-star/mrumtl/-/blob/master/doc/mrumtl.5.txt)
  fileformat. A material is mapped to a part of the OBJ geometry by using the
  `usemtl` directive of the
  [OBJ fileformat](https://gitlab.com/meso-star/htrdr/-/blob/master/doc/htrdr-obj.5.txt).
- Improve the sampling of the spectral dimension: the per wavelength
  realisation is now precisely sampled rather than arbitrarly fixed to the
  center of the sampled spectral band. Consequently, high
  resolution data defined per wavelength (e.g. Mie's properties and the
  reflectivity of the materials) are now fully taken into account.

#### Fixes

- Fix a deadlock when `htrdr` was run through MPI.
- Fix a memory leak: the output file was not closed on exit.

### Version 0.3

- Add the `-O` option that defines the file where the sky data are cached. If
  the file does not exist, the sky data structures are built from scratch and
  serialized into this new file. If this file exists, these data structures are
  directly read from it, leading to a huge speed up of the `htrdr`
  pre-processing step. Note that if the provided file exists but is filled with
  data that do not match the submitted HTGOP, HTCP and HTMie files, an error is
  detected and the program stops.
- Rely on the HTSky library to manage the sky data. This library handles the
  code previously defined into the `htrdr_sky.<c|h>` files. The HTCP, HTGOP,
  HTMie libraries are thus no more dependencies of `htrdr` since only the
  `htrdr_sky` files used them.

### Version 0.2

- Add the `-b` option that controls the BRDF of the ground geometry.
- Make optional the use of a ground geometry (option `-g`).
- Make optional the definition of the optical properties of water droplets
  (option `-m`) when no cloud field is used.

### Version 0.1

- Add the `-V` option that fixes the maximum definition of the octrees used to
  partitioned the radiative properties of the clouds.
- Add a per pixel estimation of the per radiative path computation time.

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

Copyright (C) 2018, 2019, 2020 [|Meso|Star>](http://www.meso-star.com)
<contact@meso-star.com>. Copyright (C) 2018, 2019 Centre National de la
Recherche Scientifique (CNRS), Université Paul Sabatier
<contact-edstar@laplace.univ-tlse.fr>. htrdr is free software released under
the GPL v3+ license: GNU GPL version 3 or later. You are welcome to
redistribute it under certain conditions; refer to the COPYING file for
details.

