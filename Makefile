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

ATMOSPHERE_LIBNAME = libhtrdr-atmosphere.a
COMBUSTION_LIBNAME = libhtrdr-combustion.a
PLANETS_LIBNAME = libhtrdr-planets.a

PKG_CONFIG_LOCAL = PKG_CONFIG_PATH="./:$${PKG_CONFIG_PATH}" $(PKG_CONFIG)

# Define macros when ATMOSPHERE is set to ENABLE
ATMOSPHERE_CFLAGS_ENABLE = $$($(PKG_CONFIG_LOCAL) --static --cflags htrdr-atmosphere)
ATMOSPHERE_LIBS_ENABLE = $$($(PKG_CONFIG_LOCAL) --static --libs htrdr-atmosphere)
ATMOSPHERE_BUILD_LIB_ENABLE = build_atmosphere
ATMOSPHERE_BUILD_CMD_ENABLE = build_htrdr_atmosphere
ATMOSPHERE_LIBNAME_ENABLE = $(ATMOSPHERE_LIBNAME)

# Define macros when COMBUSTION is set to ENABLE
COMBUSTION_CFLAGS_ENABLE = $$($(PKG_CONFIG_LOCAL) --static --cflags htrdr-combustion)
COMBUSTION_LIBS_ENABLE = $$($(PKG_CONFIG_LOCAL) --static --libs htrdr-combustion)
COMBUSTION_BUILD_LIB_ENABLE = build_combustion
COMBUSTION_BUILD_CMD_ENABLE = build_htrdr_combustion
COMBUSTION_LIBNAME_ENABLE = $(COMBUSTION_LIBNAME)

# Define macros when PLANETS is set to ENABLE
PLANETS_CFLAGS_ENABLE = $$($(PKG_CONFIG_LOCAL) --static --cflags htrdr-planets)
PLANETS_LIBS_ENABLE = $$($(PKG_CONFIG_LOCAL) --static --libs htrdr-planets)
PLANETS_BUILD_LIB_ENABLE = build_planets
PLANETS_BUILD_CMD_ENABLE = build_htrdr_planets
PLANETS_LIBNAME_ENABLE = $(PLANETS_LIBNAME)

# Default target
all:\
 build_htrdr\
 build_htrdr_atmosphere\
 build_htrdr_combustion\
 build_htrdr_planets\
 man

# Check commands dependencies
.config_commands: config.mk
	$(PKG_CONFIG) --atleast-version $(RSYS_VERSION) rsys
	echo 'config done' > $@

# Inference rules for command build
.SUFFIXES: .c .d .o
.c.d:
	@$(CC) $(CFLAGS_EXE) $(HTRDR_DPDC_CFLAGS) -MM -MT "$(@:.d=.o) $@" \
	$< -MF $@

.c.o:
	$(CC) $(CFLAGS_EXE) $(HTRDR_DPDC_CFLAGS) -c $< -o $@

################################################################################
# Build the htrdr command
################################################################################
HTRDR_SRC = src/commands/htrdr_cmd.c
HTRDR_OBJ = $(HTRDR_SRC:.c=.o)
HTRDR_DEP = $(HTRDR_SRC:.c=.d)

HTRDR_DPDC_CFLAGS =\
 $$($(PKG_CONFIG_LOCAL) $(PCFLAGS) --cflags htrdr-core)\
 $(ATMOSPHERE_CFLAGS_$(ATMOSPHERE))\
 $(COMBUSTION_CFLAGS_$(COMBUSTION))\
 $(PLANETS_CFLAGS_$(PLANETS))\
 $(RSYS_CFLAGS)

HTRDR_DPDC_LIBS =\
 $$($(PKG_CONFIG_LOCAL) $(PCFLAGS) --libs htrdr-core)\
 $(ATMOSPHERE_LIBS_$(ATMOSPHERE))\
 $(COMBUSTION_LIBS_$(COMBUSTION))\
 $(PLANETS_LIBS_$(PLANETS))\
 $(RSYS_LIBS)

HTRDR_DPDC_BUILD =\
 build_core\
 $(ATMOSPHERE_BUILD_LIB_$(ATMOSPHERE))\
 $(COMBUSTION_BUILD_LIB_$(COMBUSTION))\
 $(PLANETS_BUILD_LIB_$(PLANETS))

HTRDR_DPDC_PREREQ =\
 $(CORE_LIBNAME)\
 $(ATMOSPHERE_LIBNAME_$(ATMOSPHERE))\
 $(COMBUSTION_LIBNAME_$(COMBUSTION))\
 $(PLANETS_LIBNAME_$(PLANETS))

build_htrdr: .config_commands $(HTRDR_DPDC_BUILD) $(HTRDR_DEP)
	@$(MAKE) -fMakefile -f $(HTRDR_DEP) htrdr

htrdr: config.mk $(HTRDR_OBJ) $(HTRDR_DPDC_PREREQ)
	$(CC) $(CFLAGS_EXE) $(HTRDR_DPDC_CFLAGS) -o $@ \
	$(HTRDR_OBJ) $(LDFLAGS_EXE) $(HTRDR_DPDC_LIBS)

$(HTRDR_OBJ) $(HTRDR_DEP): config.mk

clean_htrdr:
	rm -f htrdr $(HTRDR_OBJ) $(HTRDR_DEP)
	rm -f .config_commands

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
	rm -f htrdr-atmosphere $(HTRDR_ATMOSPHERE_OBJ) $(HTRDR_ATMOSPHERE_DEP)
	rm -f .config_commands

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
	rm -f htrdr-combustion $(HTRDR_COMBUSTION_OBJ) $(HTRDR_COMBUSTION_DEP)
	rm -f .config_commands

################################################################################
# Build the htrdr-planets command
################################################################################
HTRDR_PLANETS_SRC = src/commands/htrdr_planets_cmd.c
HTRDR_PLANETS_OBJ = $(HTRDR_PLANETS_SRC:.c=.o)
HTRDR_PLANETS_DEP = $(HTRDR_PLANETS_SRC:.c=.d)

HTRDR_PLANETS_DPDC_LIBS = $(PLANETS_LIBS_$(PLANETS))
HTRDR_PLANETS_DPDC_BUILD = build_core $(PLANETS_BUILD_LIB_$(PLANETS))
HTRDR_PLANETS_DPDC_PREREQ = $(CORE_LIBNAME) $(PLANETS_LIBNAME_$(PLANETS))

build_htrdr_planets:\
 .config_commands\
 $(HTRDR_PLANETS_DPDC_BUILD)\
 $(HTRDR_PLANETS_DEP)
	@$(MAKE) -fMakefile -f $(HTRDR_PLANETS_DEP) htrdr-planets

htrdr-planets: config.mk $(HTRDR_PLANETS_OBJ) $(HTRDR_PLANETS_DPDC_PREREQ)
	$(CC) $(CFLAGS_EXE) $(HTRDR_DPDC_CFLAGS) -o $@ \
	$(HTRDR_PLANETS_OBJ) $(LDFLAGS_EXE) $(HTRDR_PLANETS_DPDC_LIBS)

$(HTRDR_PLANETS_OBJ) $(HTRDR_PLANETS_DEP): config.mk

clean_htrdr-planets:
	rm -f htrdr-planets $(HTRDR_PLANETS_OBJ) $(HTRDR_PLANETS_DEP)
	rm -f .config_commands

################################################################################
# Building the core
################################################################################
CORE_LIBNAME_STATIC = libhtrdr-core.a
CORE_LIBNAME_SHARED = libhtrdr-core.so
CORE_LIBNAME = $(CORE_LIBNAME_$(LIB_TYPE))

CORE_SRC =\
 src/core/htrdr.c\
 src/core/htrdr_args.c\
 src/core/htrdr_buffer.c\
 src/core/htrdr_draw_map.c\
 src/core/htrdr_geometry.c\
 src/core/htrdr_log.c\
 src/core/htrdr_materials.c\
 src/core/htrdr_proc_work.c\
 src/core/htrdr_ran_wlen_cie_xyz.c\
 src/core/htrdr_ran_wlen_discrete.c\
 src/core/htrdr_ran_wlen_planck.c\
 src/core/htrdr_rectangle.c\
 src/core/htrdr_slab.c\
 src/core/htrdr_solve_buffer.c\
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
	$(CC) $(CFLAGS_SO) $(CORE_DPDC_CFLAGS) -o $@ $(CORE_OBJ) $(LDFLAGS_SO) $(CORE_DPDC_LIBS)

$(CORE_LIBNAME_STATIC): libhtrdr-core.o
	$(AR) -rc $@ $?
	$(RANLIB) $@

libhtrdr-core.o: $(CORE_OBJ)
	$(LD) -r $(CORE_OBJ) -o $@
	$(OBJCOPY) $(OCPFLAGS) $@

.config_core: config.mk
	$(PKG_CONFIG) --atleast-version $(AW_VERSION) aw
	$(PKG_CONFIG) --atleast-version $(MPI_VERSION) $(MPI_PC)
	$(PKG_CONFIG) --atleast-version $(MRUMTL_VERSION) mrumtl
	$(PKG_CONFIG) --atleast-version $(RSYS_VERSION) rsys
	$(PKG_CONFIG) --atleast-version $(S3D_VERSION) s3d
	$(PKG_CONFIG) --atleast-version $(SCAM_VERSION) scam
	$(PKG_CONFIG) --atleast-version $(SSF_VERSION) ssf
	$(PKG_CONFIG) --atleast-version $(SSP_VERSION) star-sp
	echo 'config done' > $@

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
	@$(CC) $(CFLAGS_SO) $(CORE_DPDC_CFLAGS) -Isrc -MM -MT "$(@:.d=.o) $@" $(@:.d=.c) -MF $@

$(CORE_OBJ):
	$(CC) $(CFLAGS_SO) $(CORE_DPDC_CFLAGS) -Isrc -DHTRDR_SHARED_BUILD -c $(@:.o=.c) -o $@

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
	rm -f $(CORE_LIBNAME) $(CORE_OBJ) $(CORE_DEP)
	rm -f libhtrdr-core.o .config_core htrdr-core.pc
	rm -f src/core/htrdr_args.h src/core/htrdr_version.h

################################################################################
# Building the atmosphere library
################################################################################
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

$(ATMOSPHERE_LIBNAME): libhtrdr-atmosphere.o
	$(AR) -rc $@ $?
	$(RANLIB) $@

libhtrdr-atmosphere.o: $(ATMOSPHERE_OBJ)
	$(LD) -r $(ATMOSPHERE_OBJ) -o $@
	$(OBJCOPY) $(OCPFLAGS) $@

.config_atmosphere: config.mk
	$(PKG_CONFIG) --atleast-version $(HTSKY_VERSION) htsky
	$(PKG_CONFIG) --atleast-version $(RSYS_VERSION) rsys
	$(PKG_CONFIG) --atleast-version $(S3D_VERSION) s3d
	$(PKG_CONFIG) --atleast-version $(SCAM_VERSION) scam
	$(PKG_CONFIG) --atleast-version $(SSF_VERSION) ssf
	$(PKG_CONFIG) --atleast-version $(SSP_VERSION) star-sp
	$(PKG_CONFIG) --atleast-version $(SVX_VERSION) svx
	echo 'config done' > $@

src/atmosphere/htrdr_atmosphere_args.h: config.mk src/atmosphere/htrdr_atmosphere_args.h.in
	sed -e 's/@HTRDR_ATMOSPHERE_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD@/$(HTRDR_ATMOSPHERE_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD)/g' \
	    -e 's/@HTRDR_ATMOSPHERE_ARGS_DEFAULT_SKY_MTL_NAME@/$(HTRDR_ATMOSPHERE_ARGS_DEFAULT_SKY_MTL_NAME)/g' \
	    $@.in > $@

$(ATMOSPHERE_DEP):
	@$(CC) $(CFLAGS_SO) $(ATMOSPHERE_DPDC_CFLAGS) -Isrc -MM -MT "$(@:.d=.o) $@" $(@:.d=.c) -MF $@

$(ATMOSPHERE_OBJ):
	$(CC) $(CFLAGS_SO) $(ATMOSPHERE_DPDC_CFLAGS) -Isrc -DHTRDR_SHARED_BUILD -c $(@:.o=.c) -o $@

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
	rm -f $(ATMOSPHERE_LIBNAME) $(ATMOSPHERE_OBJ) $(ATMOSPHERE_DEP)
	rm -f .config_atmosphere libhtrdr-atmosphere.o htrdr-atmosphere.pc
	rm -f src/atmosphere/htrdr_atmosphere_args.h

################################################################################
# Building the combustion library
################################################################################
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

$(COMBUSTION_LIBNAME): libhtrdr-combustion.o
	$(AR) -rc $@ $?
	$(RANLIB) $@

libhtrdr-combustion.o: $(COMBUSTION_OBJ)
	$(LD) -r $(COMBUSTION_OBJ) -o $@
	$(OBJCOPY) $(OCPFLAGS) $@

.config_combustion: config.mk
	$(PKG_CONFIG) --atleast-version $(ATRSTM_VERSION) atrstm
	$(PKG_CONFIG) --atleast-version $(RSYS_VERSION) rsys
	$(PKG_CONFIG) --atleast-version $(S3D_VERSION) s3d
	$(PKG_CONFIG) --atleast-version $(SCAM_VERSION) scam
	$(PKG_CONFIG) --atleast-version $(SSF_VERSION) ssf
	$(PKG_CONFIG) --atleast-version $(SSP_VERSION) star-sp
	$(PKG_CONFIG) --atleast-version $(SVX_VERSION) svx
	echo 'config done' > $@

src/combustion/htrdr_combustion_args.h: config.mk src/combustion/htrdr_combustion_args.h.in
	sed -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_LASER_FLUX_DENSITY@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_LASER_FLUX_DENSITY)/g' \
	    -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_FRACTAL_DIMENSION@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_FRACTAL_DIMENSION)/g' \
	    -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_FRACTAL_PREFACTOR@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_FRACTAL_PREFACTOR)/g' \
	    -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD)/g' \
	    -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_GRID_DEFINITION_HINT@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_GRID_DEFINITION_HINT)/g' \
	    -e 's/@HTRDR_COMBUSTION_ARGS_DEFAULT_WAVELENGTH@/$(HTRDR_COMBUSTION_ARGS_DEFAULT_WAVELENGTH)/g' \
	    $@.in > $@

$(COMBUSTION_DEP):
	@$(CC) $(CFLAGS_SO) $(COMBUSTION_DPDC_CFLAGS) -Isrc -MM -MT "$(@:.d=.o) $@" $(@:.d=.c) -MF $@

$(COMBUSTION_OBJ):
	$(CC) $(CFLAGS_SO) $(COMBUSTION_DPDC_CFLAGS) -Isrc -DHTRDR_SHARED_BUILD -c $(@:.o=.c) -o $@

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
	rm -f $(COMBUSTION_LIBNAME) $(COMBUSTION_OBJ) $(COMBUSTION_DEP)
	rm -f .config_combustion libhtrdr-combustion.o htrdr-combustion.pc
	rm -f src/combustion/htrdr_combustion_args.h

################################################################################
# Building the planets library
################################################################################
PLANETS_LIBNAME = libhtrdr-planets.a

PLANETS_SRC =\
 src/planets/htrdr_planets.c\
 src/planets/htrdr_planets_args.c\
 src/planets/htrdr_planets_compute_radiance.c\
 src/planets/htrdr_planets_draw_map.c\
 src/planets/htrdr_planets_main.c\
 src/planets/htrdr_planets_source.c
PLANETS_OBJ = $(PLANETS_SRC:.c=.o)
PLANETS_DEP = $(PLANETS_SRC:.c=.d)

build_planets: build_core .config_planets htrdr-planets.pc $(PLANETS_DEP)
	@$(MAKE) -fMakefile $$(for i in $(PLANETS_DEP); do echo -f $${i}; done) \
	$(PLANETS_LIBNAME)

$(PLANETS_DEP) $(PLANETS_OBJ): config.mk src/planets/htrdr_planets_args.h

$(PLANETS_LIBNAME): libhtrdr-planets.o
	$(AR) -rc $@ $?
	$(RANLIB) $@

libhtrdr-planets.o: $(PLANETS_OBJ)
	$(LD) -r $(PLANETS_OBJ) -o $@
	$(OBJCOPY) $(OCPFLAGS) $@

.config_planets: config.mk
	$(PKG_CONFIG) --atleast-version $(RNATM_VERSION) rnatm
	$(PKG_CONFIG) --atleast-version $(RNGRD_VERSION) rngrd
	$(PKG_CONFIG) --atleast-version $(RSYS_VERSION) rsys
	$(PKG_CONFIG) --atleast-version $(S3D_VERSION) s3d
	$(PKG_CONFIG) --atleast-version $(SBUF_VERSION) sbuf
	$(PKG_CONFIG) --atleast-version $(SCAM_VERSION) scam
	$(PKG_CONFIG) --atleast-version $(SSF_VERSION) ssf
	$(PKG_CONFIG) --atleast-version $(SSP_VERSION) star-sp
	$(PKG_CONFIG) --atleast-version $(SVX_VERSION) svx
	echo 'config done' > $@

src/planets/htrdr_planets_args.h: config.mk src/planets/htrdr_planets_args.h.in
	sed -e 's/@HTRDR_PLANETS_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD@/$(HTRDR_PLANETS_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD)/g' \
	    -e 's/@HTRDR_PLANETS_ARGS_DEFAULT_GRID_DEFINITION_HINT@/$(HTRDR_PLANETS_ARGS_DEFAULT_GRID_DEFINITION_HINT)/g' \
	    -e 's/@HTRDR_PLANETS_ARGS_DEFAULT_VOLRAD_BUDGET_SPT@/$(HTRDR_PLANETS_ARGS_DEFAULT_VOLRAD_BUDGET_SPT)/g' \
	    $@.in > $@

$(PLANETS_DEP):
	@$(CC) $(CFLAGS_SO) $(PLANETS_DPDC_CFLAGS) -Isrc -MM -MT "$(@:.d=.o) $@" $(@:.d=.c) -MF $@

$(PLANETS_OBJ):
	$(CC) $(CFLAGS_SO) $(PLANETS_DPDC_CFLAGS) -Isrc -DHTRDR_SHARED_BUILD -c $(@:.o=.c) -o $@

htrdr-planets.pc: config.mk htrdr-planets.pc.in
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

clean_planets:
	rm -f $(PLANETS_LIBNAME) $(PLANETS_OBJ) $(PLANETS_DEP)
	rm -f .config_planets libhtrdr-planets.o htrdr-planets.pc
	rm -f src/planets/htrdr_planets_args.h

################################################################################
# Man pages
################################################################################
man: doc/htrdr-atmosphere.1 doc/htrdr-combustion.1 doc/htrdr-planets.1

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

doc/htrdr-planets.1: doc/htrdr-planets.1.in
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
	    -e 's/@HTRDR_PLANETS_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD@/$(HTRDR_PLANETS_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD)/g' \
	    -e 's/@HTRDR_PLANETS_ARGS_DEFAULT_GRID_DEFINITION_HINT@/$(HTRDR_PLANETS_ARGS_DEFAULT_GRID_DEFINITION_HINT)/g' \
	    -e 's/@HTRDR_PLANETS_ARGS_DEFAULT_VOLRAD_BUDGET_SPT@/$(HTRDR_PLANETS_ARGS_DEFAULT_VOLRAD_BUDGET_SPT)/g' \
	    $@.in > $@

clean_man:
	rm -f doc/htrdr-atmosphere.1 doc/htrdr-combustion.1 doc/htrdr-planets.1

################################################################################
# Installation
################################################################################
install: all
	@$(SHELL) install.sh 755 "$(DESTDIR)$(BINPREFIX)" htrdr
	@$(SHELL) install.sh 755 "$(DESTDIR)$(BINPREFIX)" htrdr-atmosphere
	@$(SHELL) install.sh 755 "$(DESTDIR)$(BINPREFIX)" htrdr-combustion
	@$(SHELL) install.sh 755 "$(DESTDIR)$(BINPREFIX)" htrdr-planets
	@if [ "$(LIB_TYPE)" = "SHARED" ]; then \
	 $(SHELL) install.sh 755 "$(DESTDIR)$(LIBPREFIX)" $(CORE_LIBNAME_SHARED); fi
	@$(SHELL) install.sh 644 "$(DESTDIR)$(DOCPREFIX)/htrdr" COPYING README.md
	@$(SHELL) install.sh 644 "$(DESTDIR)$(MANPREFIX)/man1" doc/htrdr.1
	@if [ "$(ATMOSPHERE)" = "ENABLE" ]; then \
	 $(SHELL) install.sh 644 "$(DESTDIR)$(MANPREFIX)/man1" doc/htrdr-atmosphere.1; fi
	@if [ "$(COMBUSTION)" = "ENABLE" ]; then \
	 $(SHELL) install.sh 644 "$(DESTDIR)$(MANPREFIX)/man1" doc/htrdr-combustion.1; fi
	@if [ "$(PLANETS)" = "ENABLE" ]; then \
	 $(SHELL) install.sh 644 "$(DESTDIR)$(MANPREFIX)/man1" doc/htrdr-planets.1; fi
	@$(SHELL) install.sh 644 "$(DESTDIR)$(MANPREFIX)/man5" doc/htrdr-image.5
	@$(SHELL) install.sh 644 "$(DESTDIR)$(MANPREFIX)/man5" doc/htrdr-materials.5
	@$(SHELL) install.sh 644 "$(DESTDIR)$(MANPREFIX)/man5" doc/htrdr-obj.5
	@$(SHELL) install.sh 644 "$(DESTDIR)$(MANPREFIX)/man5" doc/rnrl.5

uninstall:
	rm -f "$(DESTDIR)$(BINPREFIX)/htrdr"
	rm -f "$(DESTDIR)$(BINPREFIX)/htrdr-atmosphere"
	rm -f "$(DESTDIR)$(BINPREFIX)/htrdr-combustion"
	rm -f "$(DESTDIR)$(BINPREFIX)/htrdr-planets"
	rm -f "$(DESTDIR)$(LIBPREFIX)/$(CORE_LIBNAME_SHARED)"
	rm -f "$(DESTDIR)$(DOCPREFIX)/htrdr/COPYING"
	rm -f "$(DESTDIR)$(DOCPREFIX)/htrdr/README.md"
	rm -f "$(DESTDIR)$(MANPREFIX)/man1/htrdr.1"
	rm -f "$(DESTDIR)$(MANPREFIX)/man1/htrdr-atmosphere.1"
	rm -f "$(DESTDIR)$(MANPREFIX)/man1/htrdr-combustion.1"
	rm -f "$(DESTDIR)$(MANPREFIX)/man1/htrdr-planets.1"
	rm -f "$(DESTDIR)$(MANPREFIX)/man5/htrdr-image.5"
	rm -f "$(DESTDIR)$(MANPREFIX)/man5/htrdr-materials.5"
	rm -f "$(DESTDIR)$(MANPREFIX)/man5/htrdr-obj.5"
	rm -f "$(DESTDIR)$(MANPREFIX)/man5/rnrl.5"

################################################################################
# Miscellaneous targets
################################################################################
clean:\
 clean_htrdr\
 clean_htrdr-atmosphere\
 clean_htrdr-combustion\
 clean_htrdr-planets\
 clean_atmosphere\
 clean_combustion\
 clean_planets\
 clean_core\
 clean_man

lint: doc/htrdr-atmosphere.1 doc/htrdr-combustion.1 doc/htrdr-planets.1
	shellcheck -o all install.sh
	mandoc -Tlint -Wall doc/htrdr.1 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/htrdr-atmosphere.1 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/htrdr-combustion.1 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/htrdr-planets.1 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/htrdr-image.5 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/htrdr-materials.5 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/htrdr-obj.5 || [ $$? -le 1 ]
	mandoc -Tlint -Wall doc/rnrl.5 || [ $$? -le 1 ]
