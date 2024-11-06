#!/bin/sh

# Copyright (C) 2018-2019, 2022-2024 Centre National de la Recherche Scientifique
# Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
# Copyright (C) 2022-2024 Institut Pierre-Simon Laplace
# Copyright (C) 2022-2024 Institut de Physique du Globe de Paris
# Copyright (C) 2018-2024 |Méso|Star> (contact@meso-star.com)
# Copyright (C) 2022-2024 Observatoire de Paris
# Copyright (C) 2022-2024 Université de Reims Champagne-Ardenne
# Copyright (C) 2022-2024 Université de Versaille Saint-Quentin
# Copyright (C) 2018-2019, 2022-2024 Université Paul Sabatier
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

set -e

mode=$1
prefix=$2
shift 2

mkdir -p "${prefix}"

for i in "$@"; do
  dst="${prefix}/${i##*/}"

  if cmp -s "${i}" "${dst}"; then
    printf 'Up to date %s\n' "${dst}"
  else
    printf 'Installing %s\n' "${dst}"
    cp "${i}" "${prefix}"
    chmod "${mode}" "${prefix}/$(basename "${i}")"
  fi
done
