htrdr-materials(5)

; Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
; Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
; Copyright (C) 2022-2023 Institut Pierre-Simon Laplace
; Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
; Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
; Copyright (C) 2022-2023 Observatoire de Paris
; Copyright (C) 2022-2023 Université de Reims Champagne-Ardenne
; Copyright (C) 2022-2023 Université de Versaille Saint-Quentin
; Copyright (C) 2018-2019, 2022-2023 Université Paul Sabatier
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

htrdr-materials - list of materials used by the geometries in htrdr(1)

# DESCRIPTION

A *htrdr-materials* file lists the materials that can be used by geometries
provided through a *htrdr-obj*(5) file to the *htrdr*(1) program.  Each line
of the file gives the name of the material. For opaque materials the material
name is followed by the path toward the *mrumtl*(5) file storing the spectral
properties of its associated Bidirectional Reflectance Distribution Function.
Furthermore, the temperature of the material must be defined too.

The material name can be composed of any characters excepted for spaces and
tabulations. The path toward the *mrumtl*(5) file must be a valid path
relative to the file system.

Characters behind the hash mark (#) are considered as comments and are thus
ignored. Empty lines, i.e. lines without any characters or composed of spaces
and tabulations, are simply skipped.

# GRAMMAR

```
<htrdr-materials>     ::= <material>
                        [ <material> ... ]
<material>            ::= <name> <properties>
<name>                ::= STRING
<properties>          ::= none | <mrumtl-path> <temperature>
<mrumtl-path>         ::= PATH
<temperature>         ::= REAL # In Kelvin
```

# EXAMPLE

The following file lists 3 materials. The first one named _grass_ has its
spectral BRDF defines in the _A001.mrumtl_ file. The second one is named
_sand_ and has its spectral properties saved in the _B002.mrumtl_ file. Both
materials have a temperature of 300 Kelvin. The last material is a
semi-transparent material named _air_ with no additionnal properties defined
in this file.

```
grass   /opt/materials/A001.mrumtl 300
sand    /opt/materials/B002.mrumtl 300
air     none
```

# SEE ALSO

*htrdr*(1), *htrdr-obj*(5), *mrumtl*(5)
