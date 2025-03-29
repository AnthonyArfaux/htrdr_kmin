# Copyright (C) 2018-2019, 2022-2025 Centre National de la Recherche Scientifique
# Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
# Copyright (C) 2022-2025 Institut Pierre-Simon Laplace
# Copyright (C) 2022-2025 Institut de Physique du Globe de Paris
# Copyright (C) 2018-2025 |Méso|Star> (contact@meso-star.com)
# Copyright (C) 2022-2025 Observatoire de Paris
# Copyright (C) 2022-2025 Université de Reims Champagne-Ardenne
# Copyright (C) 2022-2025 Université de Versaille Saint-Quentin
# Copyright (C) 2018-2019, 2022-2025 Université Paul Sabatier
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

.POSIX:
.SUFFIXES: # Clean up default inference rules

include config.mk

default install uninstall lint clean:
	@$(MAKE) -fMakefile.core $@
	@if [ "$(ATMOSPHERE)" == "ENABLE" ]; then $(MAKE) -fMakefile.atmosphere $@; fi
	@if [ "$(COMBUSTION)" == "ENABLE" ]; then $(MAKE) -fMakefile.combustion $@; fi
	@if [ "$(PLANETS)" == ENABLE ]; then $(MAKE) -fMakefile.planets $@; fi

test:
	@if [ "$(COMBUSTION)" == "ENABLE" ]; then $(MAKE) -fMakefile.combustion $@; fi
	@if [ "$(PLANETS)" == ENABLE ]; then $(MAKE) -fMakefile.planets $@; fi

default: htrdr
install: install_common
uninstall: uninstall_common
lint: lint_common
clean: clean_common

htrdr: src/htrdr.in
	sed -e "s/@MAJOR@/$(VERSION_MAJOR)/" \
	    -e "s/@MINOR@/$(VERSION_MINOR)/" \
	    -e "s/@PATCH@/$(VERSION_PATCH)/" \
	    src/htrdr.in > $@

install_common: htrdr
	@$(SHELL) install.sh 755 "$(DESTDIR)$(BINPREFIX)" htrdr
	@$(SHELL) install.sh 644 "$(DESTDIR)$(DOCPREFIX)/htrdr" COPYING README.md
	@$(SHELL) install.sh 644 "$(DESTDIR)$(MANPREFIX)/man1" doc/htrdr.1
	@$(SHELL) install.sh 644 "$(DESTDIR)$(MANPREFIX)/man5" doc/htrdr-image.5
	@$(SHELL) install.sh 644 "$(DESTDIR)$(MANPREFIX)/man5" doc/htrdr-materials.5
	@$(SHELL) install.sh 644 "$(DESTDIR)$(MANPREFIX)/man5" doc/htrdr-obj.5
	@$(SHELL) install.sh 644 "$(DESTDIR)$(MANPREFIX)/man5" doc/rnrl.5

uninstall_common:
	rm -f "$(DESTDIR)$(BINPREFIX)/htrdr"
	rm -f "$(DESTDIR)$(DOCPREFIX)/htrdr/COPYING"
	rm -f "$(DESTDIR)$(DOCPREFIX)/htrdr/README.md"
	rm -f "$(DESTDIR)$(MANPREFIX)/man1/htrdr.1"
	rm -f "$(DESTDIR)$(MANPREFIX)/man5/htrdr-image.5"
	rm -f "$(DESTDIR)$(MANPREFIX)/man5/htrdr-materials.5"
	rm -f "$(DESTDIR)$(MANPREFIX)/man5/htrdr-obj.5"
	rm -f "$(DESTDIR)$(MANPREFIX)/man5/rnrl.5"

clean_common:
	rm -f htrdr

lint_common: htrdr
	shellcheck -o all install.sh
	shellcheck -o all htrdr
	mandoc -Tlint -Wall doc/htrdr.1 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/htrdr-image.5 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/htrdr-materials.5 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/htrdr-obj.5 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/rnrl.5 || [ $$? -le 1 ]
