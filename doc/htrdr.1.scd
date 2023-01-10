htrdr(1)

; Copyright (C) 2018, 2019, 2020, 2021 |Meso|Star> (contact@meso-star.com)
; Copyright (C) 2018, 2019, 2021 CNRS
; Copyright (C) 2018, 2019 Université Paul Sabatier
;                          (contact-edstar@laplace.univ-tlse.fr)
; 
; This program is free software: you can redistribute it and/or modify
; it under the terms of the GNU General Public License as published by
; the Free Software Foundation, either version 3 of the License, or
; (at your option) any later version.
; 
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
; GNU General Public License for more details.
; 
; You should have received a copy of the GNU General Public License
; along with this program. If not, see <http://www.gnu.org/licenses/>.

# NAME

htrdr - the Monte-Carlo radiative transfert simulator

# SYNOPSIS

htrdr [--version] [--help] <_mode_> [<_args_>]

# DESCRIPTION

*htrdr* computes images of scenes composed of semi transparent materials and
geometries with spectral varying reflectivities. The way the scenes are
described and the solved Monte-Carlo algorithms are controlled by the _mode_
argument that defines the computation context.

# OPTIONS

*--help*
	Print short help and exit.

*--version*
	Display version information and exit.

# MODES

Each *htrdr* _mode_ is actually a mapping to a seperate command line named
_htrdr-<mode>_ with its own set of arguments. Refer to their corresponding
man page for their complete description.

The available _htrdr-<mode>_ commands are:

*htrdr-atmosphere*(1)
	Radiative transfer computations in a cloudy atmosphere.

*htrdr-combustion*(1)
	Radiative transfer computations in a combustion medium.

*htrdr-planeto*(1)
	Radiative transfer computations in 3D planetary atmosphere.

# COPYRIGHT

Copyright © 2018, 2019, 2020, 2021 |Méso|Star> <contact@meso-star.com>++
Copyright © 2018, 2019, 2021 Centre National de la Recherche Scientifique++
Copyright © 2018, 2019 Université Paul Sabatier <contact-edstar@laplace.univ-tlse.fr>

# LICENSE

*htrdr* is free software released under the GPLv3+ license: GNU GPL version 3
or later <https://gnu.org/licenses/gpl.html>. You are free to change and
redistribute it. There is NO WARRANTY, to the extent permitted by law.

# SEE ALSO

*htrdr-atmosphere*(1),
*htrdr-combustion*(1),
*htrdr-planeto*(1)
