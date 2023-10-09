# Copyright (C) 2018-2019, 2022-2023 Centre National de la Recherche Scientifique
# Copyright (C) 2020-2022 Institut Mines Télécom Albi-Carmaux
# Copyright (C) 2022-2023 Institut Pierre-Simon Laplace
# Copyright (C) 2022-2023 Institut de Physique du Globe de Paris
# Copyright (C) 2018-2023 |Méso|Star> (contact@meso-star.com)
# Copyright (C) 2022-2023 Observatoire de Paris
# Copyright (C) 2022-2023 Université de Reims Champagne-Ardenne
# Copyright (C) 2022-2023 Université de Versaille Saint-Quentin
# Copyright (C) 2018-2019, 2022-2023 Université Paul Sabatier
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

PKG_CONFIG_LOCAL = PKG_CONFIG_PATH="./:$${PKG_CONFIG_PATH}" $(PKG_CONFIG)

# Define macros when ATMOSPHERE is set to ENABLE
ATMOSPHERE_CFLAGS_ENABLE = $$($(PKG_CONFIG_LOCAL) --static --cflags htrdr-atmosphere)
ATMOSPHERE_LIBS_ENABLE = $$($(PKG_CONFIG_LOCAL) --static --libs htrdr-atmosphere)
ATMOSPHERE_BUILD_LIB_ENABLE = build_atmosphere
ATMOSPHERE_BUILD_CMD_ENABLE = build_htrdr_atmosphere
ATMOSPHERE_LIBNAME_ENABLE = libhtrdr-atmosphere.a

# Define macros when COMBUSTION is set to ENABLE
COMBUSTION_CFLAGS_ENABLE = $$($(PKG_CONFIG_LOCAL) --static --cflags htrdr-combustion)
COMBUSTION_LIBS_ENABLE = $$($(PKG_CONFIG_LOCAL) --static --libs htrdr-combustion)
COMBUSTION_BUILD_LIB_ENABLE = build_combustion
COMBUSTION_BUILD_CMD_ENABLE = build_htrdr_combustion
COMBUSTION_LIBNAME_ENABLE = libhtrdr-combustion.a

# Define macros when PLANETO is set to ENABLE
PLANETO_CFLAGS_ENABLE = $$($(PKG_CONFIG_LOCAL) --static --cflags htrdr-planeto)
PLANETO_LIBS_ENABLE = $$($(PKG_CONFIG_LOCAL) --static --libs htrdr-planeto)
PLANETO_BUILD_LIB_ENABLE = build_planeto
PLANETO_BUILD_CMD_ENABLE = build_htrdr_planeto
PLANETO_LIBNAME_ENABLE = libhtrdr-planeto.a

# Default target
all:\
 build_htrdr\
 build_htrdr_atmosphere\
 build_htrdr_combustion\
 build_htrdr_planeto\
 man

# Check commands dependencies
.config_commands: config.mk
	@if ! $(PKG_CONFIG) --atleast-version $(RSYS_VERSION) rsys; then \
	  echo "rsys $(RSYS_VERSION) not found" >&2; exit 1; fi
	@echo "config done" > $@


# Inference rules for command build
.SUFFIXES: .c .d .o
.c.d:
	@$(CC) $(CFLAGS_EXE) $(HTRDR_DPDC_CFLAGS) -Isrc -MM -MT "$(@:.d=.o) $@" \
	$< -MF $@

.c.o:
	$(CC) $(CFLAGS_EXE) $(HTRDR_DPDC_CFLAGS) -Isrc -c $< -o $@

################################################################################
# Build the htrdr command
################################################################################
HTRDR_SRC = src/commands/htrdr_cmd.c
HTRDR_OBJ = $(HTRDR_SRC:.c=.o)
HTRDR_DEP = $(HTRDR_SRC:.c=.d)

HTRDR_DPDC_CFLAGS =\
 $(RSYS_CFLAGS)\
 $$($(PKG_CONFIG_LOCAL) $(PCFLAGS) --cflags htrdr-core)\
 $(ATMOSPHERE_CFLAGS_$(ATMOSPHERE))\
 $(COMBUSTION_CFLAGS_$(COMBUSTION))\
 $(PLANETO_CFLAGS_$(PLANETO))

HTRDR_DPDC_LIBS =\
 $(RSYS_LIBS)\
 $$($(PKG_CONFIG_LOCAL) $(PCFLAGS) --libs htrdr-core)\
 $(ATMOSPHERE_LIBS_$(ATMOSPHERE))\
 $(COMBUSTION_LIBS_$(COMBUSTION))\
 $(PLANETO_LIBS_$(PLANETO))

HTRDR_DPDC_BUILD =\
 build_core\
 $(ATMOSPHERE_BUILD_LIB_$(ATMOSPHERE))\
 $(COMBUSTION_BUILD_LIB_$(COMBUSTION))\
 $(PLANETO_BUILD_LIB_$(PLANETO))

HTRDR_DPDC_PREREQ =\
 $(CORE_LIBNAME)\
 $(ATMOSPHERE_LIBNAME_$(ATMOSPHERE))\
 $(COMBUSTION_LIBNAME_$(COMBUSTION))\
 $(PLANETO_LIBNAME_$(PLANETO))

build_htrdr: .config_commands $(HTRDR_DPDC_BUILD) $(HTRDR_DEP)
	@$(MAKE) -fMakefile -f $(HTRDR_DEP) htrdr

htrdr: config.mk $(HTRDR_OBJ) $(HTRDR_DPDC_PREREQ)
	$(CC) $(CFLAGS_EXE) $(HTRDR_DPDC_CFLAGS) -o $@ \
	$(HTRDR_OBJ) $(LDFLAGS_EXE) $(HTRDR_DPDC_LIBS)

$(HTRDR_OBJ) $(HTRDR_DEP): config.mk

clean_htrdr:
	rm -f $(HTRDR_OBJ) htrdr .config_commands

distclean_htrdr: clean_htrdr
	rm -f $(HTRDR_DEP)

################################################################################
# Build the htrdr-atmosphere command
################################################################################
HTRDR_ATMOSPHERE_SRC = src/commands/htrdr_atmosphere_cmd.c
HTRDR_ATMOSPHERE_OBJ = $(HTRDR_ATMOSPHERE_SRC:.c=.o)
HTRDR_ATMOSPHERE_DEP = $(HTRDR_ATMOSPHERE_SRC:.c=.d)

HTRDR_ATMOSPHERE_DPDC_LIBS = $(ATMOSPHERE_LIBS_$(ATMOSPHERE))
HTRDR_ATMOSPHERE_DPDC_BUILD = build_core $(ATMOSPHERE_BUILD_LIB_$(ATMOSPHERE))
HTRDR_ATMOSPHERE_DPDC_PREREQ = $(CORE_LIBNAME) $(ATMOSPHERE_LIBNAME_$(ATMOSPHERE))

build_htrdr_atmosphere:\
 .config_commands\
 $(HTRDR_ATMOSPHERE_DPDC_BUILD)\
 $(HTRDR_ATMOSPHERE_DEP)
	@$(MAKE) -fMakefile -f $(HTRDR_ATMOSPHERE_DEP) htrdr-atmosphere

htrdr-atmosphere: config.mk $(HTRDR_ATMOSPHERE_OBJ) $(HTRDR_ATMOSPHERE_DPDC_PREREQ)
	$(CC) $(CFLAGS_EXE) $(HTRDR_DPDC_CFLAGS) -o $@ \
	$(HTRDR_ATMOSPHERE_OBJ) $(LDFLAGS_EXE) $(HTRDR_ATMOSPHERE_DPDC_LIBS)

$(HTRDR_ATMOSPHERE_OBJ) $(HTRDR_ATMOSPHERE_DEP): config.mk

clean_htrdr-atmosphere:
	rm -f $(HTRDR_ATMOSPHERE_OBJ) htrdr-atmosphere

distclean_htrdr-atmosphere: clean_htrdr-atmosphere
	rm -f $(HTRDR_ATMOSPHERE_DEP) .config_commands

################################################################################
# Build the htrdr-combustion command
################################################################################
HTRDR_COMBUSTION_SRC = src/commands/htrdr_combustion_cmd.c
HTRDR_COMBUSTION_OBJ = $(HTRDR_COMBUSTION_SRC:.c=.o)
HTRDR_COMBUSTION_DEP = $(HTRDR_COMBUSTION_SRC:.c=.d)

HTRDR_COMBUSTION_DPDC_LIBS = $(COMBUSTION_LIBS_$(COMBUSTION))
HTRDR_COMBUSTION_DPDC_BUILD = build_core $(COMBUSTION_BUILD_LIB_$(COMBUSTION))
HTRDR_COMBUSTION_DPDC_PREREQ = $(CORE_LIBNAME) $(COMBUSTION_LIBNAME_$(COMBUSTION))

build_htrdr_combustion:\
 .config_commands\
 $(HTRDR_COMBUSTION_DPDC_BUILD)\
 $(HTRDR_COMBUSTION_DEP)
	@$(MAKE) -fMakefile -f $(HTRDR_COMBUSTION_DEP) htrdr-combustion

htrdr-combustion: config.mk $(HTRDR_COMBUSTION_OBJ) $(HTRDR_COMBUSTION_DPDC_PREREQ)
	$(CC) $(CFLAGS_EXE) $(HTRDR_DPDC_CFLAGS) -o $@ \
	$(HTRDR_COMBUSTION_OBJ) $(LDFLAGS_EXE) $(HTRDR_COMBUSTION_DPDC_LIBS)

$(HTRDR_COMBUSTION_OBJ) $(HTRDR_COMBUSTION_DEP): config.mk

clean_htrdr-combustion:
	rm -f $(HTRDR_COMBUSTION_OBJ) htrdr-combustion

distclean_htrdr-combustion: clean_htrdr-combustion
	rm -f $(HTRDR_COMBUSTION_DEP) .config_commands

################################################################################
# Build the htrdr-planeto command
################################################################################
HTRDR_PLANETO_SRC = src/commands/htrdr_planeto_cmd.c
HTRDR_PLANETO_OBJ = $(HTRDR_PLANETO_SRC:.c=.o)
HTRDR_PLANETO_DEP = $(HTRDR_PLANETO_SRC:.c=.d)

HTRDR_PLANETO_DPDC_LIBS = $(PLANETO_LIBS_$(PLANETO))
HTRDR_PLANETO_DPDC_BUILD = build_core $(PLANETO_BUILD_LIB_$(PLANETO))
HTRDR_PLANETO_DPDC_PREREQ = $(CORE_LIBNAME) $(PLANETO_LIBNAME_$(PLANETO))

build_htrdr_planeto:\
 .config_commands\
 $(HTRDR_PLANETO_DPDC_BUILD)\
 $(HTRDR_PLANETO_DEP)
	@$(MAKE) -fMakefile -f $(HTRDR_PLANETO_DEP) htrdr-planeto

htrdr-planeto: config.mk $(HTRDR_PLANETO_OBJ) $(HTRDR_PLANETO_DPDC_PREREQ)
	$(CC) $(CFLAGS_EXE) $(HTRDR_DPDC_CFLAGS) -o $@ \
	$(HTRDR_PLANETO_OBJ) $(LDFLAGS_EXE) $(HTRDR_PLANETO_DPDC_LIBS)

$(HTRDR_PLANETO_OBJ) $(HTRDR_PLANETO_DEP): config.mk

clean_htrdr-planeto:
	rm -f $(HTRDR_PLANETO_OBJ) htrdr-planeto

distclean_htrdr-planeto: clean_htrdr-planeto
	rm -f $(HTRDR_PLANETO_DEP) .config_commands

################################################################################
# Building the core
################################################################################
CORE_LIBNAME_STATIC = libhtrdr-core.a
CORE_LIBNAME_SHARED = libhtrdr-core.so
CORE_LIBNAME = $(CORE_LIBNAME_$(LIB_TYPE))

CORE_CFLAGS_SHARED = $(CFLAGS_SO)
CORE_CFLAGS_STATIC = $(CFLAGS_EXE)
CORE_CFLAGS = $(CORE_CFLAGS_$(LIB_TYPE))

CORE_SRC =\
 src/core/htrdr.c\
 src/core/htrdr_args.c\
 src/core/htrdr_buffer.c\
 src/core/htrdr_draw_map.c\
 src/core/htrdr_geometry.c\
 src/core/htrdr_log.c\
 src/core/htrdr_materials.c\
 src/core/htrdr_ran_wlen_cie_xyz.c\
 src/core/htrdr_ran_wlen_discrete.c\
 src/core/htrdr_ran_wlen_planck.c\
 src/core/htrdr_rectangle.c\
 src/core/htrdr_slab.c\
 src/core/htrdr_spectral.c
CORE_OBJ = $(CORE_SRC:.c=.o)
CORE_DEP = $(CORE_SRC:.c=.d)

build_core: .config_core htrdr-core.pc $(CORE_DEP)
	@$(MAKE) -fMakefile $$(for i in $(CORE_DEP); do echo -f $${i}; done) \
	$$(if [ -n "$(CORE_LIBNAME)" ]; then \
	     echo "$(CORE_LIBNAME)"; \
	   else \
	     echo "$(CORE_LIBNAME_SHARED)"; \
	   fi)

$(CORE_DEP) $(CORE_OBJ): config.mk src/core/htrdr_args.h src/core/htrdr_version.h

$(CORE_LIBNAME_SHARED): $(CORE_OBJ)
	$(CC) $(CORE_CFLAGS) $(CORE_DPDC_CFLAGS) -o $@ $(CORE_OBJ) $(LDFLAGS_SO) $(CORE_DPDC_LIBS)

$(CORE_LIBNAME_STATIC): $(CORE_OBJ)
	$(AR) -rc $@ $?
	$(RANLIB) $@

.config_core: config.mk
	@if ! $(PKG_CONFIG) --atleast-version $(AW_VERSION) aw; then \
	  echo "aw $(AW_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(MPI_VERSION) $(MPI_PC); then \
	  echo "$(MPI_PC) $(MPI_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(MRUMTL_VERSION) mrumtl; then \
	  echo "mrumtl $(MRUMTL_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(RSYS_VERSION) rsys; then \
	  echo "rsys $(RSYS_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(S3D_VERSION) s3d; then \
	  echo "s3d $(S3D_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(SCAM_VERSION) scam; then \
	  echo "scam $(SCAM_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(SSF_VERSION) ssf; then \
	  echo "ssf $(SSF_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(SSP_VERSION) star-sp; then \
	  echo "star-sp $(SSP_VERSION) not found" >&2; exit 1; fi
	@echo "config done" > $@

src/core/htrdr_args.h: config.mk src/core/htrdr_args.h.in
	sed -e 's/@HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MIN@/$(HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MIN)/g' \
	    -e 's/@HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MAX@/$(HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MAX)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_POS@/$(HTRDR_ARGS_DEFAULT_CAMERA_POS)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_TGT@/$(HTRDR_ARGS_DEFAULT_CAMERA_TGT)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_UP@/$(HTRDR_ARGS_DEFAULT_CAMERA_UP)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_ORTHOGRAPHIC_HEIGHT@/$(HTRDR_ARGS_DEFAULT_CAMERA_ORTHOGRAPHIC_HEIGHT)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOV@/$(HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOV)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_LENS_RADIUS@/$(HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_LENS_RADIUS)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOCAL_DST@/$(HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOCAL_DST)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_RECTANGLE_POS@/$(HTRDR_ARGS_DEFAULT_RECTANGLE_POS)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_RECTANGLE_TGT@/$(HTRDR_ARGS_DEFAULT_RECTANGLE_TGT)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_RECTANGLE_UP@/$(HTRDR_ARGS_DEFAULT_RECTANGLE_UP)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_RECTANGLE_SZ@/$(HTRDR_ARGS_DEFAULT_RECTANGLE_SZ)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_IMG_WIDTH@/$(HTRDR_ARGS_DEFAULT_IMG_WIDTH)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_IMG_HEIGHT@/$(HTRDR_ARGS_DEFAULT_IMG_HEIGHT)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_IMG_SPP@/$(HTRDR_ARGS_DEFAULT_IMG_SPP)/g' \
	    $@.in > $@

src/core/htrdr_version.h: config.mk src/core/htrdr_version.h.in
	sed -e 's/@VERSION_MAJOR@/$(VERSION_MAJOR)/g' \
	    -e 's/@VERSION_MINOR@/$(VERSION_MINOR)/g' \
	    -e 's/@VERSION_PATCH@/$(VERSION_PATCH)/g' \
	    $@.in > $@

$(CORE_DEP):
	@$(CC) $(CORE_CFLAGS) $(CORE_DPDC_CFLAGS) -Isrc -MM -MT "$(@:.d=.o) $@" $(@:.d=.c) -MF $@

$(CORE_OBJ):
	$(CC) $(CORE_CFLAGS) $(CORE_DPDC_CFLAGS) -Isrc -DHTRDR_CORE_SHARED_BUILD -c $(@:.o=.c) -o $@

htrdr-core.pc: config.mk htrdr-core.pc.in
	sed -e 's/@VERSION@/$(VERSION)/g' \
	    -e 's/@AW_VERSION@/$(AW_VERSION)/g' \
	    -e 's/@MPI_PC@/$(MPI_PC)/g' \
	    -e 's/@MPI_VERSION@/$(MPI_VERSION)/g' \
	    -e 's/@MRUMTL_VERSION@/$(MRUMTL_VERSION)/g' \
	    -e 's/@RSYS_VERSION@/$(RSYS_VERSION)/g' \
	    -e 's/@S3D_VERSION@/$(S3D_VERSION)/g' \
	    -e 's/@SCAM_VERSION@/$(SCAM_VERSION)/g' \
	    -e 's/@SSF_VERSION@/$(SSF_VERSION)/g' \
	    -e 's/@SSP_VERSION@/$(SSP_VERSION)/g' \
	    $@.in > $@

clean_core:
	rm -f $(CORE_OBJ) $(CORE_LIBNAME)
	rm -f .config_core src/core/htrdr_args.h src/core/htrdr_version.h
	rm -f htrdr-core.pc

distclean_core: clean_core
	rm -f $(CORE_DEP)

################################################################################
# Building the atmosphere library
################################################################################
ATMOSPHERE_LIBNAME = libhtrdr-atmosphere.a

ATMOSPHERE_SRC =\
 src/atmosphere/htrdr_atmosphere_args.c\
 src/atmosphere/htrdr_atmosphere.c\
 src/atmosphere/htrdr_atmosphere_compute_radiance_lw.c\
 src/atmosphere/htrdr_atmosphere_compute_radiance_sw.c\
 src/atmosphere/htrdr_atmosphere_draw_map.c\
 src/atmosphere/htrdr_atmosphere_ground.c\
 src/atmosphere/htrdr_atmosphere_main.c\
 src/atmosphere/htrdr_atmosphere_sun.c
ATMOSPHERE_OBJ = $(ATMOSPHERE_SRC:.c=.o)
ATMOSPHERE_DEP = $(ATMOSPHERE_SRC:.c=.d)

build_atmosphere: build_core .config_atmosphere htrdr-atmosphere.pc $(ATMOSPHERE_DEP)
	@$(MAKE) -fMakefile $$(for i in $(ATMOSPHERE_DEP); do echo -f $${i}; done) \
	$(ATMOSPHERE_LIBNAME)

$(ATMOSPHERE_DEP) $(ATMOSPHERE_OBJ): config.mk src/atmosphere/htrdr_atmosphere_args.h

$(ATMOSPHERE_LIBNAME): $(ATMOSPHERE_OBJ)
	$(AR) -rc $@ $?
	$(RANLIB) $@

.config_atmosphere: config.mk
	@if ! $(PKG_CONFIG) --atleast-version $(HTSKY_VERSION) htsky; then \
	  echo "htsky $(HTSKY_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(RSYS_VERSION) rsys; then \
	  echo "rsys $(RSYS_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(S3D_VERSION) s3d; then \
	  echo "s3d $(S3D_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(SCAM_VERSION) scam; then \
	  echo "scam $(SCAM_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(SSF_VERSION) ssf; then \
	  echo "ssf $(SSF_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(SSP_VERSION) star-sp; then \
	  echo "star-sp $(SSP_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(SVX_VERSION) svx; then \
	  echo "svx $(SVX_VERSION) not found" >&2; exit 1; fi
	@echo "config done" > $@

src/atmosphere/htrdr_atmosphere_args.h: config.mk src/atmosphere/htrdr_atmosphere_args.h.in
	sed -e 's/@HTRDR_ATMOSPHERE_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD@/$(HTRDR_ATMOSPHERE_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD)/g' \
	    -e 's/@HTRDR_ATMOSPHERE_ARGS_DEFAULT_SKY_MTL_NAME@/$(HTRDR_ATMOSPHERE_ARGS_DEFAULT_SKY_MTL_NAME)/g' \
	    $@.in > $@

$(ATMOSPHERE_DEP):
	@$(CC) $(CFLAGS_EXE) $(ATMOSPHERE_DPDC_CFLAGS) -Isrc -MM -MT "$(@:.d=.o) $@" $(@:.d=.c) -MF $@

$(ATMOSPHERE_OBJ):
	$(CC) $(CFLAGS_EXE) $(ATMOSPHERE_DPDC_CFLAGS) -Isrc -DHTRDR_STATIC -c $(@:.o=.c) -o $@

htrdr-atmosphere.pc: config.mk htrdr-atmosphere.pc.in
	sed -e 's/@VERSION@/$(VERSION)/g' \
	    -e 's/@HTSKY_VERSION@/$(HTSKY_VERSION)/g' \
	    -e 's/@RSYS_VERSION@/$(RSYS_VERSION)/g' \
	    -e 's/@S3D_VERSION@/$(S3D_VERSION)/g' \
	    -e 's/@SCAM_VERSION@/$(SCAM_VERSION)/g' \
	    -e 's/@SSF_VERSION@/$(SSF_VERSION)/g' \
	    -e 's/@SSP_VERSION@/$(SSP_VERSION)/g' \
	    -e 's/@SVX_VERSION@/$(SVX_VERSION)/g' \
	    $@.in > $@

clean_atmosphere:
	rm -f $(ATMOSPHERE_OBJ) $(ATMOSPHERE_LIBNAME)
	rm -f .config_atmosphere src/atmosphere/htrdr_atmosphere_args.h
	rm -f htrdr-atmosphere.pc

distclean_atmosphere: clean_atmosphere
	rm -f $(ATMOSPHERE_DEP)

################################################################################
# Building the combustion library
################################################################################
COMBUSTION_LIBNAME = libhtrdr-combustion.a

COMBUSTION_SRC =\
 src/combustion/htrdr_combustion.c\
 src/combustion/htrdr_combustion_args.c\
 src/combustion/htrdr_combustion_draw_map.c\
 src/combustion/htrdr_combustion_compute_radiance_sw.c\
 src/combustion/htrdr_combustion_geometry_ray_filter.c\
 src/combustion/htrdr_combustion_laser.c\
 src/combustion/htrdr_combustion_main.c\
 src/combustion/htrdr_combustion_phase_func.c
COMBUSTION_OBJ = $(COMBUSTION_SRC:.c=.o)
COMBUSTION_DEP = $(COMBUSTION_SRC:.c=.d)

build_combustion: build_core .config_combustion htrdr-combustion.pc $(COMBUSTION_DEP)
	@$(MAKE) -fMakefile $$(for i in $(COMBUSTION_DEP); do echo -f $${i}; done) \
	$(COMBUSTION_LIBNAME)

$(COMBUSTION_DEP) $(COMBUSTION_OBJ): config.mk src/combustion/htrdr_combustion_args.h

$(COMBUSTION_LIBNAME): $(COMBUSTION_OBJ)
	$(AR) -rc $@ $?
	$(RANLIB) $@

.config_combustion: config.mk
	@if ! $(PKG_CONFIG) --atleast-version $(ATRSTM_VERSION) atrstm; then \
	  echo "atrstm $(ATRSTM_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(RSYS_VERSION) rsys; then \
	  echo "rsys $(RSYS_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(S3D_VERSION) s3d; then \
	  echo "s3d $(S3D_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(SCAM_VERSION) scam; then \
	  echo "scam $(SCAM_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(SSF_VERSION) ssf; then \
	  echo "ssf $(SSF_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(SSP_VERSION) star-sp; then \
	  echo "star-sp $(SSP_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(SVX_VERSION) svx; then \
	  echo "svx $(SVX_VERSION) not found" >&2; exit 1; fi
	@echo "config done" > $@

src/combustion/htrdr_combustion_args.h: config.mk src/combustion/htrdr_combustion_args.h.in
	sed -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_LASER_FLUX_DENSITY@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_LASER_FLUX_DENSITY)/g' \
	    -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_FRACTAL_DIMENSION@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_FRACTAL_DIMENSION)/g' \
	    -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_FRACTAL_PREFACTOR@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_FRACTAL_PREFACTOR)/g' \
	    -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD)/g' \
	    -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_GRID_DEFINITION_HINT@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_GRID_DEFINITION_HINT)/g' \
	    -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_WAVELENGTH@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_WAVELENGTH)/g' \
	    $@.in > $@

$(COMBUSTION_DEP):
	@$(CC) $(CFLAGS_EXE) $(COMBUSTION_DPDC_CFLAGS) -Isrc -MM -MT "$(@:.d=.o) $@" $(@:.d=.c) -MF $@

$(COMBUSTION_OBJ):
	$(CC) $(CFLAGS_EXE) $(COMBUSTION_DPDC_CFLAGS) -Isrc -DHTRDR_STATIC -c $(@:.o=.c) -o $@


htrdr-combustion.pc: config.mk htrdr-combustion.pc.in
	sed -e 's/@VERSION@/$(VERSION)/g' \
	    -e 's/@ATRSTM_VERSION@/$(ATRSTM_VERSION)/g' \
	    -e 's/@RSYS_VERSION@/$(RSYS_VERSION)/g' \
	    -e 's/@S3D_VERSION@/$(S3D_VERSION)/g' \
	    -e 's/@SCAM_VERSION@/$(SCAM_VERSION)/g' \
	    -e 's/@SSF_VERSION@/$(SSF_VERSION)/g' \
	    -e 's/@SSP_VERSION@/$(SSP_VERSION)/g' \
	    -e 's/@SVX_VERSION@/$(SVX_VERSION)/g' \
	    $@.in > $@

clean_combustion:
	rm -f $(COMBUSTION_OBJ) $(COMBUSTION_LIBNAME)
	rm -f .config_combustion src/combustion/htrdr_combustion_args.h
	rm -f htrdr-combustion.pc

distclean_combustion: clean_combustion
	rm -f $(COMBUSTION_DEP)

################################################################################
# Building the planeto library
################################################################################
PLANETO_LIBNAME = libhtrdr-planeto.a

PLANETO_SRC =\
 src/planeto/htrdr_planeto.c\
 src/planeto/htrdr_planeto_args.c\
 src/planeto/htrdr_planeto_compute_radiance.c\
 src/planeto/htrdr_planeto_draw_map.c\
 src/planeto/htrdr_planeto_main.c\
 src/planeto/htrdr_planeto_source.c
PLANETO_OBJ = $(PLANETO_SRC:.c=.o)
PLANETO_DEP = $(PLANETO_SRC:.c=.d)

build_planeto: build_core .config_planeto htrdr-planeto.pc $(PLANETO_DEP)
	@$(MAKE) -fMakefile $$(for i in $(PLANETO_DEP); do echo -f $${i}; done) \
	$(PLANETO_LIBNAME)

$(PLANETO_DEP) $(PLANETO_OBJ): config.mk src/planeto/htrdr_planeto_args.h

$(PLANETO_LIBNAME): $(PLANETO_OBJ)
	$(AR) -rc $@ $?
	$(RANLIB) $@

.config_planeto: config.mk
	@if ! $(PKG_CONFIG) --atleast-version $(RNATM_VERSION) rnatm; then \
	  echo "rnatm $(RNATM_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(RNGRD_VERSION) rngrd; then \
	  echo "rngrd $(RNGRD_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(RSYS_VERSION) rsys; then \
	  echo "rsys $(RSYS_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(S3D_VERSION) s3d; then \
	  echo "s3d $(S3D_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(SBUF_VERSION) sbuf; then \
	  echo "sbuf $(SBUF_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(SCAM_VERSION) scam; then \
	  echo "scam $(SCAM_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(SSF_VERSION) ssf; then \
	  echo "ssf $(SSF_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(SSP_VERSION) star-sp; then \
	  echo "star-sp $(SSP_VERSION) not found" >&2; exit 1; fi
	@if ! $(PKG_CONFIG) --atleast-version $(SVX_VERSION) svx; then \
	  echo "svx $(SVX_VERSION) not found" >&2; exit 1; fi
	@echo "config done" > $@

src/planeto/htrdr_planeto_args.h: config.mk src/planeto/htrdr_planeto_args.h.in
	sed -e 's/@HTRDR_PLANETO_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD@/$(HTRDR_PLANETO_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD)/g' \
	    -e 's/@HTRDR_PLANETO_ARGS_DEFAULT_GRID_DEFINITION_HINT@/$(HTRDR_PLANETO_ARGS_DEFAULT_GRID_DEFINITION_HINT)/g' \
	    $@.in > $@

$(PLANETO_DEP):
	@$(CC) $(CFLAGS_EXE) $(PLANETO_DPDC_CFLAGS) -Isrc -MM -MT "$(@:.d=.o) $@" $(@:.d=.c) -MF $@

$(PLANETO_OBJ):
	$(CC) $(CFLAGS_EXE) $(PLANETO_DPDC_CFLAGS) -Isrc -DHTRDR_STATIC -c $(@:.o=.c) -o $@

htrdr-planeto.pc: config.mk htrdr-planeto.pc.in
	sed -e 's/@VERSION@/$(VERSION)/g' \
	    -e 's/@RNATM_VERSION@/$(RNATM_VERSION)/g' \
	    -e 's/@RNGRD_VERSION@/$(RNGRD_VERSION)/g' \
	    -e 's/@RSYS_VERSION@/$(RSYS_VERSION)/g' \
	    -e 's/@S3D_VERSION@/$(S3D_VERSION)/g' \
	    -e 's/@SBUF_VERSION@/$(SBUF_VERSION)/g' \
	    -e 's/@SCAM_VERSION@/$(SCAM_VERSION)/g' \
	    -e 's/@SSF_VERSION@/$(SSF_VERSION)/g' \
	    -e 's/@SSP_VERSION@/$(SSP_VERSION)/g' \
	    -e 's/@SVX_VERSION@/$(SVX_VERSION)/g' \
	    $@.in > $@

clean_planeto:
	rm -f $(PLANETO_OBJ) $(PLANETO_LIBNAME)
	rm -f .config_planeto src/planeto/htrdr_planeto_args.h
	rm -f htrdr-planeto.pc

distclean_planeto: clean_planeto
	rm -f $(PLANETO_DEP)

################################################################################
# Man pages
################################################################################
man: doc/htrdr-atmosphere.1 doc/htrdr-combustion.1 doc/htrdr-planeto.1

doc/htrdr-atmosphere.1: doc/htrdr-atmosphere.1.in
	sed -e 's/@HTRDR_ATMOSPHERE_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD@/$(HTRDR_ATMOSPHERE_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD)/g' \
	    -e 's/@HTRDR_ATMOSPHERE_ARGS_DEFAULT_SKY_MTL_NAME@/$(HTRDR_ATMOSPHERE_ARGS_DEFAULT_SKY_MTL_NAME)/g' \
	    -e 's/@HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MIN@/$(HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MIN)/g' \
	    -e 's/@HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MAX@/$(HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MAX)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOCAL_DST@/$(HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOCAL_DST)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOV@/$(HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOV)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_LENS_RADIUS@/$(HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_LENS_RADIUS)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_POS@/$(HTRDR_ARGS_DEFAULT_CAMERA_POS)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_TGT@/$(HTRDR_ARGS_DEFAULT_CAMERA_TGT)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_UP@/$(HTRDR_ARGS_DEFAULT_CAMERA_UP)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_IMG_WIDTH@/$(HTRDR_ARGS_DEFAULT_IMG_WIDTH)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_IMG_HEIGHT@/$(HTRDR_ARGS_DEFAULT_IMG_HEIGHT)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_IMG_SPP@/$(HTRDR_ARGS_DEFAULT_IMG_SPP)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_RECTANGLE_POS@/$(HTRDR_ARGS_DEFAULT_RECTANGLE_POS)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_RECTANGLE_TGT@/$(HTRDR_ARGS_DEFAULT_RECTANGLE_TGT)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_RECTANGLE_UP@/$(HTRDR_ARGS_DEFAULT_RECTANGLE_UP)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_RECTANGLE_SZ@/$(HTRDR_ARGS_DEFAULT_RECTANGLE_SZ)/g'\
	    $@.in > $@

doc/htrdr-combustion.1: doc/htrdr-combustion.1.in
	sed -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_LASER_FLUX_DENSITY@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_LASER_FLUX_DENSITY)/g' \
	    -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_FRACTAL_DIMENSION@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_FRACTAL_DIMENSION)/g' \
	    -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_FRACTAL_PREFACTOR@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_FRACTAL_PREFACTOR)/g' \
	    -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD)/g' \
	    -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_GRID_DEFINITION_HINT@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_GRID_DEFINITION_HINT)/g' \
	    -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_WAVELENGTH@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_WAVELENGTH)/g' \
	    -e 's/@HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MIN@/$(HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MIN)/g' \
	    -e 's/@HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MAX@/$(HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MAX)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOCAL_DST@/$(HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOCAL_DST)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOV@/$(HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOV)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_LENS_RADIUS@/$(HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_LENS_RADIUS)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_POS@/$(HTRDR_ARGS_DEFAULT_CAMERA_POS)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_TGT@/$(HTRDR_ARGS_DEFAULT_CAMERA_TGT)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_UP@/$(HTRDR_ARGS_DEFAULT_CAMERA_UP)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_IMG_WIDTH@/$(HTRDR_ARGS_DEFAULT_IMG_WIDTH)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_IMG_HEIGHT@/$(HTRDR_ARGS_DEFAULT_IMG_HEIGHT)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_IMG_SPP@/$(HTRDR_ARGS_DEFAULT_IMG_SPP)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_RECTANGLE_POS@/$(HTRDR_ARGS_DEFAULT_RECTANGLE_POS)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_RECTANGLE_TGT@/$(HTRDR_ARGS_DEFAULT_RECTANGLE_TGT)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_RECTANGLE_UP@/$(HTRDR_ARGS_DEFAULT_RECTANGLE_UP)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_RECTANGLE_SZ@/$(HTRDR_ARGS_DEFAULT_RECTANGLE_SZ)/g'\
	    $@.in > $@

doc/htrdr-planeto.1: doc/htrdr-planeto.1.in
	sed -e 's/@HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MIN@/$(HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MIN)/g' \
	    -e 's/@HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MAX@/$(HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MAX)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOCAL_DST@/$(HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOCAL_DST)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOV@/$(HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOV)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_LENS_RADIUS@/$(HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_LENS_RADIUS)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_POS@/$(HTRDR_ARGS_DEFAULT_CAMERA_POS)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_TGT@/$(HTRDR_ARGS_DEFAULT_CAMERA_TGT)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_CAMERA_UP@/$(HTRDR_ARGS_DEFAULT_CAMERA_UP)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_IMG_WIDTH@/$(HTRDR_ARGS_DEFAULT_IMG_WIDTH)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_IMG_HEIGHT@/$(HTRDR_ARGS_DEFAULT_IMG_HEIGHT)/g' \
	    -e 's/@HTRDR_ARGS_DEFAULT_IMG_SPP@/$(HTRDR_ARGS_DEFAULT_IMG_SPP)/g' \
	    -e 's/@HTRDR_PLANETO_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD@/$(HTRDR_PLANETO_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD)/g' \
	    -e 's/@HTRDR_PLANETO_ARGS_DEFAULT_GRID_DEFINITION_HINT@/$(HTRDR_PLANETO_ARGS_DEFAULT_GRID_DEFINITION_HINT)/g' \
	    $@.in > $@

clean_man:
	rm -f doc/htrdr-atmosphere.1 doc/htrdr-combustion.1 doc/htrdr-planeto.1

################################################################################
# Installation
################################################################################
install: all
	@$(SHELL) make.sh install "$(DESTDIR)$(PREFIX)/bin" htrdr
	@$(SHELL) make.sh install "$(DESTDIR)$(PREFIX)/bin" htrdr-atmosphere
	@$(SHELL) make.sh install "$(DESTDIR)$(PREFIX)/bin" htrdr-combustion
	@$(SHELL) make.sh install "$(DESTDIR)$(PREFIX)/bin" htrdr-planeto
	@if [ "$(LIB_TYPE)" = "SHARED" ]; then \
	 $(SHELL) make.sh install "$(DESTDIR)$(PREFIX)/lib" $(CORE_LIBNAME_SHARED); fi
	@$(SHELL) make.sh install "$(DESTDIR)$(PREFIX)/share/doc/htrdr" COPYING README.md
	@$(SHELL) make.sh install "$(DESTDIR)$(PREFIX)/share/man/man1" doc/htrdr.1
	@if [ "$(ATMOSPHERE)" = "ENABLE" ]; then \
	 $(SHELL) make.sh install "$(DESTDIR)$(PREFIX)/share/man/man1" doc/htrdr-atmosphere.1; fi
	@if [ "$(COMBUSTION)" = "ENABLE" ]; then \
	 $(SHELL) make.sh install "$(DESTDIR)$(PREFIX)/share/man/man1" doc/htrdr-combustion.1; fi
	@if [ "$(PLANETO)" = "ENABLE" ]; then \
	 $(SHELL) make.sh install "$(DESTDIR)$(PREFIX)/share/man/man1" doc/htrdr-planeto.1; fi
	@$(SHELL) make.sh install "$(DESTDIR)$(PREFIX)/share/man/man5" doc/htrdr-image.5
	@$(SHELL) make.sh install "$(DESTDIR)$(PREFIX)/share/man/man5" doc/htrdr-materials.5
	@$(SHELL) make.sh install "$(DESTDIR)$(PREFIX)/share/man/man5" doc/htrdr-obj.5
	@$(SHELL) make.sh install "$(DESTDIR)$(PREFIX)/share/man/man5" doc/rnrl.5

uninstall:
	rm -f "$(DESTDIR)$(PREFIX)/bin/htrdr"
	rm -f "$(DESTDIR)$(PREFIX)/bin/htrdr-atmosphere"
	rm -f "$(DESTDIR)$(PREFIX)/bin/htrdr-combustion"
	rm -f "$(DESTDIR)$(PREFIX)/bin/htrdr-planeto"
	rm -f "$(DESTDIR)$(PREFIX)/lib/$(CORE_LIBNAME_SHARED)"
	rm -f "$(DESTDIR)$(PREFIX)/share/doc/htrdr/COPYING"
	rm -f "$(DESTDIR)$(PREFIX)/share/doc/htrdr/README.md"
	rm -f "$(DESTDIR)$(PREFIX)/share/man/man1/htrdr.1"
	rm -f "$(DESTDIR)$(PREFIX)/share/man/man1/htrdr-atmosphere.1"
	rm -f "$(DESTDIR)$(PREFIX)/share/man/man1/htrdr-combustion.1"
	rm -f "$(DESTDIR)$(PREFIX)/share/man/man1/htrdr-planeto.1"
	rm -f "$(DESTDIR)$(PREFIX)/share/man/man5/htrdr-image.5"
	rm -f "$(DESTDIR)$(PREFIX)/share/man/man5/htrdr-materials.5"
	rm -f "$(DESTDIR)$(PREFIX)/share/man/man5/htrdr-obj.5"
	rm -f "$(DESTDIR)$(PREFIX)/share/man/man5/rnrl.5"

################################################################################
# Miscellaneous targets
################################################################################
clean:\
 clean_htrdr\
 clean_htrdr-atmosphere\
 clean_htrdr-combustion\
 clean_htrdr-planeto\
 clean_atmosphere\
 clean_combustion\
 clean_planeto\
 clean_core\
 clean_man

distclean:\
 distclean_htrdr\
 distclean_htrdr-atmosphere\
 distclean_htrdr-combustion\
 distclean_htrdr-planeto\
 distclean_atmosphere\
 distclean_combustion\
 distclean_planeto\
 distclean_core\
 clean_man

lint: doc/htrdr-atmosphere.1 doc/htrdr-combustion.1 doc/htrdr-planeto.1
	shellcheck -o all make.sh
	mandoc -Tlint -Wall doc/htrdr.1 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/htrdr-atmosphere.1 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/htrdr-combustion.1 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/htrdr-planeto.1 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/htrdr-image.5 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/htrdr-materials.5 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/htrdr-obj.5 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/rnrl.5 || [ $$? -le 1 ]
