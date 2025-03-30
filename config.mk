VERSION_MAJOR = 0
VERSION_MINOR = 11
VERSION_PATCH = 0

VERSION = $(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)

PREFIX = /usr/local
BINPREFIX = $(PREFIX)/bin
DOCPREFIX = $(PREFIX)/share/doc
INCPREFIX = $(PREFIX)/include
LIBPREFIX = $(PREFIX)/lib
MANPREFIX = $(PREFIX)/share/man

# Define the features supported, i.e. the htrdr commands to be built.
# Any value other than ENABLE disables the corresponding functionality.
# So, simply comment on a feature to deactivate it.
ATMOSPHERE = ENABLE
COMBUSTION = ENABLE
PLANETS = ENABLE

LIB_TYPE = SHARED
#LIB_TYPE = STATIC

BUILD_TYPE = RELEASE
#BUILD_TYPE = DEBUG

# MPI pkg-config file
MPI_PC = ompi

################################################################################
# Default argument values
################################################################################
# Core
HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MIN = 0.0
HTRDR_ARGS_CAMERA_PERSPECTIVE_FOV_EXCLUSIVE_MAX = 180.0
HTRDR_ARGS_DEFAULT_CAMERA_POS = 0,0,0
HTRDR_ARGS_DEFAULT_CAMERA_TGT = 0,1,0
HTRDR_ARGS_DEFAULT_CAMERA_UP = 0,0,1
HTRDR_ARGS_DEFAULT_CAMERA_ORTHOGRAPHIC_HEIGHT = 1
HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOV = 70
HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_LENS_RADIUS = 0
HTRDR_ARGS_DEFAULT_CAMERA_PERSPECTIVE_FOCAL_DST = 1
HTRDR_ARGS_DEFAULT_RECTANGLE_POS = 0,0,0
HTRDR_ARGS_DEFAULT_RECTANGLE_TGT = 0,0,1
HTRDR_ARGS_DEFAULT_RECTANGLE_UP = 0,1,0
HTRDR_ARGS_DEFAULT_RECTANGLE_SZ = 1,1
HTRDR_ARGS_DEFAULT_IMG_WIDTH = 320
HTRDR_ARGS_DEFAULT_IMG_HEIGHT = 240
HTRDR_ARGS_DEFAULT_IMG_SPP = 1

# Atmosphere
HTRDR_ATMOSPHERE_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD = 1
HTRDR_ATMOSPHERE_ARGS_DEFAULT_SKY_MTL_NAME = "air"

# Combustion
HTRDR_COMBUSTION_ARGS_DEFAULT_LASER_FLUX_DENSITY = 1
HTRDR_COMBUSTION_ARGS_DEFAULT_FRACTAL_DIMENSION = 1.80
HTRDR_COMBUSTION_ARGS_DEFAULT_FRACTAL_PREFACTOR = 1.30
HTRDR_COMBUSTION_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD = 1.0
HTRDR_COMBUSTION_ARGS_DEFAULT_GRID_DEFINITION_HINT  = 256
HTRDR_COMBUSTION_ARGS_DEFAULT_WAVELENGTH = 532

# Planets
HTRDR_PLANETS_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD = 1
HTRDR_PLANETS_ARGS_DEFAULT_GRID_DEFINITION_HINT = 512
HTRDR_PLANETS_ARGS_DEFAULT_VOLRAD_BUDGET_SPT = 1000

################################################################################
# Tools
################################################################################
AR = ar
CC = cc
LD = ld
OBJCOPY = objcopy
PKG_CONFIG = pkg-config
PKG_CONFIG_LOCAL = PKG_CONFIG_PATH="./:$${PKG_CONFIG_PATH}" $(PKG_CONFIG)
RANLIB = ranlib

################################################################################
# Dependencies
################################################################################
PCFLAGS_SHARED =
PCFLAGS_STATIC = --static
PCFLAGS = $(PCFLAGS_$(LIB_TYPE))

AW_VERSION = 2.1
ATRSTM_VERSION = 0.1
HTSKY_VERSION = 0.3
MPI_VERSION = 2
MRUMTL_VERSION = 0.2
RNATM_VERSION = 0.1
RNGRD_VERSION = 0.1
RSYS_VERSION = 0.14
S3D_VERSION = 0.10
SBUF_VERSION = 0.1
SCAM_VERSION = 0.2
SSF_VERSION = 0.9
SMSH_VERSION = 0.1
SSP_VERSION = 0.14
SVX_VERSION = 0.3

# Atmosphere
ATMOSPHERE_INCS = $$($(PKG_CONFIG_LOCAL) $(PCFLAGS) --cflags\
 htrdr-core htsky rsys s3d scam ssf star-sp svx)
ATMOSPHERE_LIBS = $$($(PKG_CONFIG_LOCAL) $(PCFLAGS) --libs\
 htrdr-core htsky rsys s3d scam ssf star-sp svx) -lm

# Combustion
COMBUSTION_INCS = $$($(PKG_CONFIG_LOCAL) $(PCFLAGS) --cflags\
 atrstm htrdr-core rsys s3d scam ssf star-sp svx)
COMBUSTION_LIBS = $$($(PKG_CONFIG_LOCAL) $(PCFLAGS) --libs\
 atrstm htrdr-core rsys s3d scam ssf star-sp svx) -lm

# Core
CORE_INCS = $$($(PKG_CONFIG_LOCAL) $(PCFLAGS) --cflags\
 aw $(MPI_PC) mrumtl rsys s3d scam ssf star-sp) -fopenmp
CORE_LIBS = $$($(PKG_CONFIG_LOCAL) $(PCFLAGS) --libs\
 aw $(MPI_PC) mrumtl rsys s3d scam ssf star-sp) -fopenmp -lm

# Planets
PLANETS_INCS = $$($(PKG_CONFIG_LOCAL) $(PCFLAGS) --cflags\
 $(MPI_PC) htrdr-core rnatm rngrd rsys s3d sbuf scam smsh ssf star-sp svx)
PLANETS_LIBS = $$($(PKG_CONFIG_LOCAL) $(PCFLAGS) --libs\
 $(MPI_PC) htrdr-core rnatm rngrd rsys s3d sbuf scam smsh ssf star-sp svx) -lm

################################################################################
# Compilation options
################################################################################
WFLAGS =\
 -Wall\
 -Wcast-align\
 -Wconversion\
 -Wextra\
 -Wmissing-declarations\
 -Wmissing-prototypes\
 -Wshadow

# Increase the security and robustness of generated binaries
CFLAGS_HARDENED =\
 -D_FORTIFY_SOURCES=2\
 -fcf-protection=full\
 -fstack-clash-protection\
 -fstack-protector-strong

CFLAGS_COMMON =\
 -std=c89\
 -pedantic\
 -fvisibility=hidden\
 -fstrict-aliasing\
 $(CFLAGS_HARDENED)\
 $(WFLAGS)

CFLAGS_DEBUG = -g $(CFLAGS_COMMON)
CFLAGS_RELEASE = -O2 -DNDEBUG $(CFLAGS_COMMON)
CFLAGS = $(CFLAGS_$(BUILD_TYPE))

CFLAGS_SO = $(CFLAGS) -fPIC
CFLAGS_EXE = $(CFLAGS) -fPIE

################################################################################
# Linker options
################################################################################
LDFLAGS_HARDENED = -Wl,-z,relro,-z,now
LDFLAGS_DEBUG = $(LDFLAGS_HARDENED)
LDFLAGS_RELEASE = -s $(LDFLAGS_HARDENED)
LDFLAGS = $(LDFLAGS_$(BUILD_TYPE))

LDFLAGS_SO = $(LDFLAGS) -shared -Wl,--no-undefined
LDFLAGS_EXE = $(LDFLAGS) -pie

OCPFLAGS_DEBUG = --localize-hidden
OCPFLAGS_RELEASE = --localize-hidden --strip-unneeded
OCPFLAGS = $(OCPFLAGS_$(BUILD_TYPE))
