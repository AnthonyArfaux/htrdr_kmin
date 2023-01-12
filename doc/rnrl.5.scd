rnrl(5)

; Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
; Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
; Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
; Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
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

rnsr - Rad-Net Radiance List file format

# DESCRIPTION

*rnrl* is a binary file format for storing a list of spectrally varying
radiances. The radiances (in W/m²/sr/m) are listed by wavelength (in nm)
sorted in ascending order.

A *rnrl* file is actually a Star-Buffer file (see *sbuf*(5)). It starts with a
header of 4 64-bit integers describing the layout of the data. The first integer
is a power of two (usually 4096) which defines the size of the memory page in
bytes on which the list of radiances by wavelength aligns (_pagesize_). The
second integer is the _size_ of the array, that is, the number of wavelengths
for which a radiance is defined. Finally, the remaining 2 integers store the
memory size (16 bytes) and memory alignment (16 bytes) of a radiance per
wavelength.

The fill bytes follow the file header to align the radiances by wavelength to
_pagesize_.

Each item in the list is composed of 2 double-precision floating-point values:
the wavelength in nanometers and its corresponding radiance in W/m²/sr/m.

The end of the file is eventually padded with dummy bytes to ensure that the
overall file size is a multiple of _pagesize_.

# BINARY FILE FORMAT

Data are encoded with respect to the little endian bytes ordering, i.e. least
significant bytes are stored first.

```
<rnsr>          ::= <pagesize> <size> 16 16
                    <padding>
                    <radiances>
                    <padding>

<pagesize>      ::= UINT64
<size>          ::= UINT64 # Number of items stored

---

<radiances>     ::= <property> ... ...
<property>      ::= <wavelength> <radiance>
<wavelength>    ::= DOUBLE # in nanometer
<radiance>      ::= DOUBLE # in W/m²/sr/m

---

<padding>       ::= [ BYTE ... ] # Ensure alignement
```

# SEE ALSO

*sbuf*(5)
