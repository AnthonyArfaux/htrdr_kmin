# `htrdr`

`htrdr` evaluates the intensity at any position (probe) of the scene, in any
direction, in the presence of surfaces and an absorbing and diffusing
semi-transparent medium, both for radiation sources that are internal to the
medium (longwave) or external to the medium (shortwave). The intensity is
calculated using the *Monte-Carlo* method: a number of optical paths are
simulated backward, from the probe position and into the medium. Various
algorithms are used, depending on the specificities of the nature and shape of
the radiation source.

Applications are theoretically possible to any configuration. However, it all
eventually comes down to the possibility of using the physical data of
interest, in their most common formats, in each scientific community. `htrdr`
is currently suitable for two main application fields:

1. *Atmospheric radiative transfer*: the clear-sky atmosphere is vertically
   stratified, cloud thermodynamic data is provided on a regular 3D rectangular
   grid, and surface optical properties can be provided for an arbitrary number
   of materials. Internal radiation and solar radiation are taken into account.

2. *Combustion* processes: thermodynamic data is provided at the nodes of an
   unstructured tetrahedral mesh, while surface properties can still be
   provided for various materials. The radiation source is only external: a
   monochromatic laser sheet illuminates the inside of the combustion chamber
   for diagnostic purposes.

Since any observable radiative transfer is expressed as an integral of the
intensity, and since there is a strict equivalence between the integral to be
solved and the underlying Monte-Carlo algorithm (each integral results in the
sampling of a random variable), the algorithms that calculate the radiance are
used for computing various quantities:

- *Images* on a camera sensor, in a given field of view. For combustion
  applications, only monochromatic images are supported. In atmospheres, both
  visible and infrared images are possible: CIE colorimetry is used for visible
  images, while an infrared image is in fact a temperature map of luminosity,
  over the required spectral interval.

- *Flux density maps*, on a grid of sensors, integrated over an entire
  hemisphere. In the case of combustion chambers, only monochromatic flux maps
  can be calculated, while spectrally integrated flux density maps (both on the
  visible part of the spectrum and on the infrared) are possible for
  atmospheric applications.

## How to build

This program is compatible GNU/Linux 64-bits. It relies on the
[CMake](http://www.cmake.org) and the
[RCMake](https://gitlab.com/vaplv/rcmake/) packages to build. It also depends
on the
[AW](https://gitlab.com/vaplv/loader_aw/#tab-readme),
[MruMtl](https://gitlab.com/meso-star/mrumtl/),
[RSys](https://gitlab.com/vaplv/rsys/),
[Star-3D](https://gitlab.com/meso-star/star-3d/),
[Star-SF](https://gitlab.com/meso-star/star-sf/),
[Star-SP](https://gitlab.com/meso-star/stat-sp/) and
[Star-VX](https://gitlab.com/meso-star/star-vx/) libraries and on
[OpenMP](http://www.openmp.org) 1.2 and the
[MPI](https://www.mpi-forum.org/) 2 specification to parallelize its
computations.

`htrdr` finally depends on the [HTSky](https://gitlab.com/meso-star/htsky/)
library if the `HTRDR_BUILD_ATMOSPHERE` option is set and on
[AtrSTM](https://gitlab.com/meso-star/atrstm/) when `HTRDR_BUILD_COMBUSTION` is
set. These options enable/disable the build of the atmospheric part and the
combustion part of `htrdr`. By default, both options are activated.

To build `htrdr`, first ensure that CMake is installed on your system. Then
install the RCMake package as well as the aforementioned prerequisites. Finally
generate the project from the `cmake/CMakeLists.txt` file by appending to the
`CMAKE_PREFIX_PATH` variable the install directories of its dependencies. The
resulting project can be edited, built, tested and installed as any CMake
project. Refer to the [CMake](https://cmake.org/documentation) for further
informations on CMake.

## Release notes

### Version 0.8.1

Sets the required version of Star-SampPling to 0.12. This version fixes
compilation errors with gcc 11 but introduces API breaks.

### Version 0.8

- Adds support for a thin lens camera model and an orthographic camera model
  for combustion and atmosphere modes.
- Updates the size of a tile from 32x32 pixels to 8x8 pixels. A tile is a
  block of pixels rendered by a thread. However, a size of 32x32 pixels could
  be too large when rendering on several dozen threads: the image definition
  could be insufficient to give tiles to all threads. 
- Fixes the calculation of shortwave radiance by `htrdr-combustion` and the
  calculation of longwave radiance by `htrdr-atmosphere`. At each scattering
  position, the range of the traced ray could be incorrect.

### Version 0.7

#### Adds the simulation of radiative transfer in combustion media

The new `htrdr-combustion` command performs radiative transfer computations in
a scene representing a semi-transparent medium enlightened by a laser sheet. It
uses Monte-Carlo to calculate a monochromatic image of the medium or the
radiative flux density. Both computations are performed in the visible at a
given frequency.

The medium data are defined on the vertices of an unstructured tetrahedral mesh
that may be surrounded by a triangular surface mesh representing the inner
limits of the combustion chamber.

#### Updates the `htrdr` command

The previous `htrdr` command is renamed to `htrdr-atmosphere`. `htrdr` becomes
a proxy for the `htrdr-atmosphere` command or the `htrdr-combustion` command:
calling `htrdr` with the `<atmosphere|combustion>` options is equivalent to
directly calling the `htrdr-<atmosphere|combustion>` commands.

#### Miscellaneous

- Major update of the entire codebase to add multiple applications to `htrdr`:
  It was originally designed to handle atmospheric applications only.
- Always displays the number of processes and the number of threads: previously
  they were only printed on multi-node executions.
- Fixed auto intersection issue on surfaces not facing the sun.
- Fixed writing of pixel data: assumed pixel layout could be wrong.

### Version 0.6.1

- Fix the self-intersection issue in shortwave computations introduced by
  the 0.6 version.

### Version 0.6

#### New features

- Add support of flux map computation for both shortwave and longwave. The flux
  is computed for the part of the flux map lying outside any geometry. The new
  command line option `-p` defines the rectangle in the scene onto which the
  flux is going to be integrated. The flux map resolution and the realisations
  count per pixel is controlled by the `-i` option.
- Add support of thin materials, i.e. materials without geometric thickness as
  for instance the leaves of the trees.
- Add the temperature property to the materials and used it as the limit
  condition during longwave computations. Previously, the surface temperatures
  were fetched from the atmosphere at the given surface position.
- Add the `-n` option to fix the name of the material defining the atmosphere.

#### Fix

- In shortwave, fix how direct contribution is handled for purely specular
  BRDF.

### Version 0.5.1

- Fix the `undefined strtok_r symbol` issue exhibited by some GCC versions that
  leads to memory corruption and segmentation fault when parsing the ground
  interfaces.
- Fix typos in the man pages.

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

## Copyright notice

Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique  
Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux  
Copyright (C) 2022-2023 Institut de Physique du Globe de Paris  
Copyright (C) 2018-2023 [|Méso|Star>](http://www.meso-star.com) (<contact@meso-star.com>)  
Copyright (C) 2022-2023 Université de Reims Champagne-Ardenne  
Copyright (C) 2022-2023 Université de Versaille Saint-Quentin  
Copyright (C) 2018-2019, 2022-2023 Université Paul Sabatier

## License

`htrdr` is free software released under the GPL v3+ license: GNU GPL version 3 or
later. You are welcome to redistribute it under certain conditions; refer to
the COPYING file for details.

