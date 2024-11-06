# `htrdr`

`htrdr` evaluates the intensity at any position (probe) of the scene, in
any direction, in the presence of *surfaces* and an *absorbing and
diffusing semi-transparent medium*, for both *internal* (longwave) or
*external* (shortwave) *radiation sources*. The intensity is calculated
using the *Monte Carlo* method: a number of optical paths are simulated
backward, from the probe position and into the medium. Various
algorithms are used, depending on the specificities of the nature and
shape of the radiation source.

Applications are theoretically possible to any configuration. However,
it all eventually comes down to the possibility of using the physical
data of interest, in their most common formats, in each scientific
community. `htrdr` is currently suitable for three main application
fields:

1. *Atmospheric radiative transfer:* a clear-sky atmosphere is
   vertically stratified, neglecting Earth sphericity, and described in
   terms of absorption coefficients as a function of height and spectral
   quadrature point as per a correlated-k model. Cloud physical
   properties are provided on a 3D rectangular grid. Surface geometrical
   and optical properties can be provided for an arbitrary number of
   geometries. Internal radiation and solar radiation are taken into
   account.

2. *Combustion* processes: thermodynamic data is provided at the nodes
   of an unstructured tetrahedral mesh, while surface properties can
   still be provided for various materials. The radiation source is only
   external: a monochromatic laser sheet illuminates the inside of the
   combustion chamber for diagnostic purposes.

3. *Planetary science*: takes into account the geometry of a "ground" of
   arbitrary shape, described by a triangular mesh, with the possibility
   of using an arbitrary number of materials. The radiative properties
   of a gas mixture must be provided on a tetrahedral mesh, using the
   k-distribution spectral model. The radiative properties of an
   arbitrary number of aerosol and hydrometeores can also be provided on
   their individual tetrahedral mesh.  Calculations can be made for both
   internal and external radiation sources.  In the case of an external
   source, a sphere of arbitrary size and position is used. This sphere
   can radiate as a Planck source at a specified brightness temperature,
   or be associated with a high-resolution radiance spectrum.

Since any radiative transfer observable is expressed as an integral of
the intensity, and since there is a strict equivalence between the
integral to be solved and the underlying Monte Carlo algorithm (each
integral results in the sampling of a random variable), the algorithms
that calculate the radiance are used for computing various quantities:

- *Images* on a camera sensor, in a given field of view. For combustion
  applications, only monochromatic images are supported. In atmospheres,
  spectral integration is also possible, both for solar and thermal
  images: CIE colorimetry is used for solar images, while thermal images
  are in fact brightness temperature maps, obtained from the incoming
  radiative flux over a specified spectral interval.

- *Flux density maps*, on a sensor grid, integrated over an entire
  hemisphere. In the case of combustion chambers, only monochromatic
  flux maps can be calculated, while spectrally integrated flux density
  maps (both on the visible part of the spectrum and on the infrared)
  are possible for atmospheric applications.

## Requirements

- C compiler with OpenMP support
- POSIX make
- pkg-config
- Message Passing Interface (MPI)
- [mandoc](https://mandoc.bsd.lv)
- [AW](https://gitlab.com/vaplv/loader_aw)
- [Astoria: Semi-Transparent Medium](https://gitlab.com/meso-star/atrstm)
  (optional; required by COMBUSTION)
- [High-Tune: Sky](https://gitlab.com/meso-star/htsky)
  (optional; required by ATMOSHPERE)
- [ModRadUrb: MaTeriaL](https://gitlab.com/meso-star/mrumtl)
- [Rad-Net ATMopshere](https://gitlab.com/meso-star/rnatm)
  (optional; required by PLANETS)
- [Rad-Net GRounD](https://gitlab.com/meso-star/rngrd)
  (optional; required by PLANETS)
- [RSys](https://gitlab.com/vaplv/rsys)
- [Star 3D](https://gitlab.com/meso-star/star-3d)
- [Star Buffer](https://gitlab.com/meso-star/star-buffer)
  (optional; required by PLANETS)
- [Star Camera](https://gitlab.com/meso-star/star-camera)
- [Star Scattering Functions](https://gitlab.com/meso-star/star-sf)
- [Star SamPling](https://gitlab.com/meso-star/star-sp)
- [Star VoXel](https://gitlab.com/meso-star/star-vx)

## Installation

Edit config.mk as needed, then run:

    make clean install

## Release notes

### Version 0.10

#### Use POSIX make as a build system

The build procedure is now written in POSIX make instead of CMake.
In addition to the features already provided by its CMake alternative,
the Makefile supports the use of static libraries and provides an
uninstall target. In any case, the main motivation behind its writing is
to use a good old well-established standard with simple features,
available on all UNIX systems, thus simplifying its portability and
support while being much lighter

#### Proof-reading and editing manual pages

Write the man pages directly in mdoc's roff macros, instead of using the
scdoc markup language as a source for man pages.

Unlike writing manuals with man's roff macros, and even more so with
scdoc, mdoc macros take care of layout, font handling and all the other
typesetting details which, by construction, guarantee the consistency of
all manuals without leaving the responsibility to the individual author.
This also facilitates translation into other formats and documentation
tools. These are the main reasons for writing manual pages with mdoc
macros.

A complete re-reading of the manual pages was carried out during the
translation into mdoc, with several corrections and rewrites to make the
manual clearer.

#### Bug fixes

- Fix the construction of the planck/CIE cumulative: the memory space
  required could be prohibitive, leading to a shortage of memory space.
- Update error handling in `htrdr_ran_wlen_planck_create` to avoid
  freeing up the same memory space several times, and thus causing a
  crash.
- All dependencies have been updated. In particular, the `htmi` and
  `rnatm` libraries have been updated to version 0.1. They fixe several
  problems with the atmosphere and planeto commands. See their release
  note for more information.

### Version 0.9.2

- Update the `rnatm` library to version 0.0.1. This versions fixes
  several bugs when different atmospheric components do not have the
  same volumetric meshes.
- Display an error message when parsing unknown arguments to the
  `htrdr-planeto` command.
- Fix `htrdr-planeto` man page: there was an error in the `csplit`
  command provided as an example to extract octrees from the output of
  `htrdr-planeto` when the -d option is used.

### Version 0.9.1

- Fix invalid read/write memory access when ray tracing the ground in
  `htrdr-atmopshere`.
- Fix compilation warning detected by GCC 12.
- Fix `htrdr-planeto` man page (-S option): the unit of the radius and
  the distance from the source is not the meter but the kilometer.
- Fix `htrdr` man page: replaced long options with short options.
- Reference and install the rnrl fileformat man page.
- Proofreading the README and man pages: correcting typos, spelling and
  formatting errors and turns of phrase

### Version 0.9

#### Adds radiative transfer simulation in 3D planetary atmospheres

The new `htrdr-planeto` command simulates radiative transfer in
planetology context, i.e. in the 3D atmosphere of a telluric planet.
Both infrared and visible computations are supported. `htrdr-planeto` is
actually a renderer that calculates an image for a given observation
position. Its internal rendering algorithm is based on Monte Carlo
integration, which consists for each pixel in simulating a given number
of optical paths from the sensor, taking into account the phenomena of
light absorption and scattering.

The planet's ground can be any set of triangles with BRDFs and
temperatures defined per triangle. The atmosphere is composed of a gas
mixture and a potentially empty set of aerosols. Both can have arbitrary
tetrahedral meshes with per-node radiative properties.

#### Miscellaneous

- Use scdoc rather than asciidoc as file format for man sources.
- Update all dependencies. More notably, use
  [MruMtl 0.1](https://gitlab.com/meso-star/mrumtl/-/tree/0.1)
  which introduces API breaks.
- Add the discrete wavelength distribution currently used in
  `htrdr-planeto` only.

### Version 0.8.1

Sets the required version of Star-SampPling to 0.12. This version fixes
compilation errors with gcc 11 but introduces API breaks.

### Version 0.8

- Adds support for a thin lens camera model and an orthographic camera
  model for combustion and atmosphere modes.
- Updates the size of a tile from 32x32 pixels to 8x8 pixels. A tile is
  a block of pixels rendered by a thread. However, a size of 32x32
  pixels could be too large when rendering on several dozen threads: the
  image definition could be insufficient to give tiles to all threads.
- Fixes the calculation of shortwave radiance by `htrdr-combustion` and
  the calculation of longwave radiance by `htrdr-atmosphere`. At each
  scattering position, the range of the traced ray could be incorrect.

### Version 0.7

#### Adds the simulation of radiative transfer in combustion media

The new `htrdr-combustion` command performs radiative transfer
computations in a scene representing a semi-transparent medium
enlightened by a laser sheet. It uses Monte Carlo to calculate a
monochromatic image of the medium or the radiative flux density. Both
computations are performed in the visible at a given frequency.

The medium data are defined on the vertices of an unstructured
tetrahedral mesh that may be surrounded by a triangular surface mesh
representing the inner limits of the combustion chamber.

#### Updates the `htrdr` command

The previous `htrdr` command is renamed to `htrdr-atmosphere`. `htrdr`
becomes a proxy for the `htrdr-atmosphere` command or the
`htrdr-combustion` command: calling `htrdr` with the
`<atmosphere|combustion>` options is equivalent to directly calling the
`htrdr-<atmosphere|combustion>` commands.

#### Miscellaneous

- Major update of the entire codebase to add multiple applications to
  `htrdr`: It was originally designed to handle atmospheric applications
  only.
- Always displays the number of processes and the number of threads:
  previously they were only printed on multi-node executions.
- Fixed auto intersection issue on surfaces not facing the sun.
- Fixed writing of pixel data: assumed pixel layout could be wrong.

### Version 0.6.1

- Fix the self-intersection issue in shortwave computations introduced
  by the 0.6 version.

### Version 0.6

#### New features

- Add support of flux map computation for both shortwave and longwave.
  The flux is computed for the part of the flux map lying outside any
  geometry. The new command line option `-p` defines the rectangle in
  the scene onto which the flux is going to be integrated. The flux map
  resolution and the realisations count per pixel is controlled by the
  `-i` option.
- Add support of thin materials, i.e. materials without geometric
  thickness as for instance the leaves of the trees.
- Add the temperature property to the materials and used it as the limit
  condition during longwave computations. Previously, the surface
  temperatures were fetched from the atmosphere at the given surface
  position.
- Add the `-n` option to fix the name of the material defining the
  atmosphere.

#### Fix

- In shortwave, fix how direct contribution is handled for purely
  specular BRDF.

### Version 0.5.1

- Fix the `undefined strtok_r symbol` issue exhibited by some GCC
  versions that leads to memory corruption and segmentation fault when
  parsing the ground interfaces.
- Fix typos in the man pages.

### Version 0.5

#### New feature

Add support of shortwave integration with respect to the Planck function
for a reference temperature whose default value is the blackbody
temperature of the sun. Actually this is the counterpart of the longwave
integration introduced by the "infrared rendering" in the 0.4 version.
The main difference is that the source of radiation is the sun rather
than the medium and its boundaries.

The option `-l` that enabled the infrared rendering is now replaced by
the new `-s` option that controls the spectral integration that can be
CIE XYZ (i.e.  regular image rendering), longwave or shortwave.

#### Fixes

- Fix the returned sun radiance: the precomputed per spectral band solar
  incoming flux is removed and the sun radiance is now retrieved by
  directly evaluating the monochromatic Planck for the blackbody
  temperature of the sun.
- Fix CIE XYZ spectral integration: the pdf used to sample the CIE
  tristimulus values was not correctly handled in the Monte Carlo
  weight.
- Fix the longwave spectral integration: the Monte Carlo weight was
  wrong leading to overestimated temperatures.

### Version 0.4

#### New features

- Add support of infrared rendering: when defined, the new `-l` option
  setups the range of longwave into which the rendering is performed. In
  infrared rendering, each pixel stores the radiance per pixel and its
  associated brightness temperature. Spectral integration is done with
  respect to the Planck function for a reference temperature of 290 K.
- The ground geometry can now have several materials whose data vary
  over the spectrum. These materials are listed in a new
  [htrdr-materials](https://gitlab.com/meso-star/htrdr/-/blob/master/doc/htrdr-materials.5.txt)
  file where each materials is defined by a name and a file storing its
  spectral data with respect to the
  [MruMtl](https://gitlab.com/meso-star/mrumtl/-/blob/master/doc/mrumtl.5.txt)
  fileformat. A material is mapped to a part of the OBJ geometry by
  using the `usemtl` directive of the
  [OBJ fileformat](https://gitlab.com/meso-star/htrdr/-/blob/master/doc/htrdr-obj.5.txt).
- Improve the sampling of the spectral dimension: the per wavelength
  realisation is now precisely sampled rather than arbitrarly fixed to
  the center of the sampled spectral band. Consequently, high resolution
  data defined per wavelength (e.g. Mie's properties and the
  reflectivity of the materials) are now fully taken into account.

#### Fixes

- Fix a deadlock when `htrdr` was run through MPI.
- Fix a memory leak: the output file was not closed on exit.

### Version 0.3

- Add the `-O` option that defines the file where the sky data are
  cached. If the file does not exist, the sky data structures are built
  from scratch and serialized into this new file. If this file exists,
  these data structures are directly read from it, leading to a huge
  speed up of the `htrdr` pre-processing step. Note that if the provided
  file exists but is filled with data that do not match the submitted
  HTGOP, HTCP and HTMie files, an error is detected and the program
  stops.
- Rely on the HTSky library to manage the sky data. This library handles
  the code previously defined into the `htrdr_sky.<c|h>` files. The
  HTCP, HTGOP, HTMie libraries are thus no more dependencies of `htrdr`
  since only the `htrdr_sky` files used them.

### Version 0.2

- Add the `-b` option that controls the BRDF of the ground geometry.
- Make optional the use of a ground geometry (option `-g`).
- Make optional the definition of the optical properties of water
  droplets (option `-m`) when no cloud field is used.

### Version 0.1

- Add the `-V` option that fixes the maximum definition of the octrees
  used to partitioned the radiative properties of the clouds.
- Add a per pixel estimation of the per radiative path computation time.

### Version 0.0.4

- Fix the computation of the surface scattering: there was a bug in how
  Russian roulette was implemented at surface scattering leading to an
  underestimation of the surface reflection.
- Update the thread allocation policy: by default, the number of threads
  is now defined as the maximum between the number of processors
  detected by OpenMP and the number of threads defined by the
  `OMP_NUM_THREADS` environment variable. This variable can be used to
  counteract the number of processors detected by OpenMP that can be
  lower than the real number of processors of the system.

### Version 0.0.3

- Fix compilation on systems with a GNU C Library whose version is less
  than 2.19.
- Fix a possible invalid memory access to cloud data leading to
  segmentation faults.

## Copyright notice

Copyright © 2018-2019, 2022-2024 Centre National de la Recherche Scientifique  
Copyright © 2020-2022 Institut Mines Télécom Albi-Carmaux  
Copyright © 2022-2024 Institut Pierre-Simon Laplace  
Copyright © 2022-2024 Institut de Physique du Globe de Paris  
Copyright © 2018-2024 [|Méso|Star>](http://www.meso-star.com) (contact@meso-star.com)  
Copyright © 2022-2024 Observatoire de Paris  
Copyright © 2022-2024 Université de Reims Champagne-Ardenne  
Copyright © 2022-2024 Université de Versaille Saint-Quentin  
Copyright © 2018-2019, 2022-2024 Université Paul Sabatier

## License

`htrdr` is free software released under the GPL v3+ license: GNU GPL
version 3 or later. You are welcome to redistribute it under certain
conditions; refer to the COPYING file for details.
