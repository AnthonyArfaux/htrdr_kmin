# High-Tune: RenDeRer

This programs is used to test the implementation of radiative transfer
Monte-Carlo algorithms in cloudy atmopsheres.

## How to build

This program is compatible GNU/Linux 64-bits. It relies on the
[CMake](http://www.cmake.org) and the
[RCMake](https://gitlab.com/vaplv/rcmake/) packages to build.  It also depends
on the
[HTCP](https://gitlab.com/meso-star/htcp/),
[HTMIE](https://gitlab.com/meso-star/htmie/),
[RSys](https://gitlab.com/vaplv/rsys/),
[Star-3D](https://gitlab.com/meso-star/star-3d/),
[Star-3DAW](https://gitlab.com/meso-star/star-3daw/),
[Star-SF](https://gitlab.com/meso-star/star-sf/),
[Star-SP](https://gitlab.com/meso-star/stat-sp/) and
[Star-VX](https://gitlab.com/meso-star/star-vx/) libraries and on the
[OpenMP](http://www.openmp.org) 1.2 specification to parallelize its
computations.

First ensure that CMake is installed on your system. Then install the RCMake
package as well as the aforementioned prerequisites. Finally generate the
project from the `cmake/CMakeLists.txt` file by appending to the
`CMAKE_PREFIX_PATH` variable the install directories of its dependencies. The
resulting project can be edited, built, tested and installed as any CMake
project. Refer to the [CMake](https://cmake.org/documentation) for further
informations on CMake.

## License

htrdr is a free software copyright (C) 2018 Université Paul Sabatier
<contact-edstar@laplace.univ-tlse.fr>, |Meso|Star <contact@meso-star.com>. It
is released under the GPL v3+ license: GNU GPL version 3 or later. You are
welcome to redistribute it under certain conditions; refer to the COPYING file
for details.

