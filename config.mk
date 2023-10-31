VERSION_MAJOR = 0
VERSION_MINOR = 10
VERSION_PATCH = 0

VERSION = $(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)
PREFIX = /usr/local

LIB_TYPE = SHARED
#LIB_TYPE = STATIC

BUILD_TYPE = RELEASE
#BUILD_TYPE = DEBUG

# MPI pkg-config file
MPI_PC = ompi

# Define the features supported, i.e. the htrdr commands to be built.
# Any value other than ENABLE disables the corresponding functionality.
# So, simply comment on a feature to deactivate it.
ATMOSPHERE = ENABLE
COMBUSTION = ENABLE
PLANETO = ENABLE

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

# Planeto
HTRDR_PLANETO_ARGS_DEFAULT_OPTICAL_THICKNESS_THRESHOLD = 1
HTRDR_PLANETO_ARGS_DEFAULT_GRID_DEFINITION_HINT = 512

################################################################################
# Tools
################################################################################
AR = ar
CC = cc
LD = ld
OBJCOPY = objcopy
PKG_CONFIG = pkg-config
RANLIB = ranlib

################################################################################
# Dependencies
################################################################################
PCFLAGS_SHARED =
PCFLAGS_STATIC = --static
PCFLAGS = $(PCFLAGS_$(LIB_TYPE))

AW_VERSION = 2.1
AW_CFLAGS = $$($(PKG_CONFIG) $(PCFLAGS) --cflags aw)
AW_LIBS = $$($(PKG_CONFIG) $(PCFLAGS) --libs aw)

ATRSTM_VERSION = 0.1
ATRSTM_CFLAGS = $$($(PKG_CONFIG) $(PCFLAGS) --cflags atrstm)
ATRSTM_LIBS = $$($(PKG_CONFIG) $(PCFLAGS) --libs atrstm)

HTSKY_VERSION = 0.3
HTSKY_CFLAGS = $$($(PKG_CONFIG) $(PCFLAGS) --cflags htsky)
HTSKY_LIBS = $$($(PKG_CONFIG) $(PCFLAGS) --libs htsky)

MPI_VERSION = 2
MPI_CFLAGS = $$($(PKG_CONFIG) $(PCFLAGS) --cflags $(MPI_PC))
MPI_LIBS = $$($(PKG_CONFIG) $(PCFLAGS) --libs $(MPI_PC))

MRUMTL_VERSION = 0.2
MRUMTL_CFLAGS = $$($(PKG_CONFIG) $(PCFLAGS) --cflags mrumtl)
MRUMTL_LIBS = $$($(PKG_CONFIG) $(PCFLAGS) --libs mrumtl)

RNATM_VERSION = 0.1
RNATM_CFLAGS = $$($(PKG_CONFIG) $(PCFLAGS) --cflags rnatm)
RNATM_LIBS = $$($(PKG_CONFIG) $(PCFLAGS) --libs rnatm)

RNGRD_VERSION = 0.1
RNGRD_CFLAGS = $$($(PKG_CONFIG) $(PCFLAGS) --cflags rngrd)
RNGRD_LIBS = $$($(PKG_CONFIG) $(PCFLAGS) --libs rngrd)

RSYS_VERSION = 0.14
RSYS_CFLAGS = $$($(PKG_CONFIG) $(PCFLAGS) --cflags rsys)
RSYS_LIBS = $$($(PKG_CONFIG) $(PCFLAGS) --libs rsys)

S3D_VERSION = 0.10
S3D_CFLAGS = $$($(PKG_CONFIG) $(PCFLAGS) --cflags s3d)
S3D_LIBS = $$($(PKG_CONFIG) $(PCFLAGS) --libs s3d)

SBUF_VERSION = 0.1
SBUF_CFLAGS = $$($(PKG_CONFIG) $(PCFLAGS) --cflags sbuf)
SBUF_LIBS = $$($(PKG_CONFIG) $(PCFLAGS) --libs sbuf)

SCAM_VERSION = 0.2
SCAM_CFLAGS = $$($(PKG_CONFIG) $(PCFLAGS) --cflags scam)
SCAM_LIBS = $$($(PKG_CONFIG) $(PCFLAGS) --libs scam)

SSF_VERSION = 0.9
SSF_CFLAGS = $$($(PKG_CONFIG) $(PCFLAGS) --cflags ssf)
SSF_LIBS = $$($(PKG_CONFIG) $(PCFLAGS) --libs ssf)

SSP_VERSION = 0.14
SSP_CFLAGS = $$($(PKG_CONFIG) $(PCFLAGS) --cflags star-sp)
SSP_LIBS = $$($(PKG_CONFIG) $(PCFLAGS) --libs star-sp)

SVX_VERSION = 0.3
SVX_CFLAGS = $$($(PKG_CONFIG) $(PCFLAGS) --cflags svx)
SVX_LIBS = $$($(PKG_CONFIG) $(PCFLAGS) --libs svx)

# Atmosphere
ATMOSPHERE_DPDC_CFLAGS =\
 $(HTSKY_CFLAGS)\
 $(RSYS_CFLAGS)\
 $(S3D_CFLAGS)\
 $(SCAM_CFLAGS)\
 $(SSF_CFLAGS)\
 $(SSP_CFLAGS)\
 $(SVX_CFLAGS)
ATMOSPHERE_DPDC_LIBS =\
 $(HTSKY_LIBS)\
 $(RSYS_LIBS)\
 $(S3D_LIBS)\
 $(SCAM_LIBS)\
 $(SSF_LIBS)\
 $(SSP_LIBS)\
 $(SVX_LIBS)\
 -lm

# Combustion
COMBUSTION_DPDC_CFLAGS =\
 $(ATRSTM_CFLAGS)\
 $(RSYS_CFLAGS)\
 $(S3D_CFLAGS)\
 $(SCAM_CFLAGS)\
 $(SSF_CFLAGS)\
 $(SSP_CFLAGS)\
 $(SVX_CFLAGS)
COMBUSTION_DPDC_LIBS =\
 $(ATRSTM_LIBS)\
 $(RSYS_LIBS)\
 $(S3D_LIBS)\
 $(SCAM_LIBS)\
 $(SSF_LIBS)\
 $(SSP_LIBS)\
 $(SVX_LIBS)\
 -lm

# Core
CORE_DPDC_CFLAGS =\
 $(AW_CFLAGS)\
 $(MPI_CFLAGS)\
 $(MRUMTL_CFLAGS)\
 $(RSYS_CFLAGS)\
 $(S3D_CFLAGS)\
 $(SCAM_CFLAGS)\
 $(SSF_CFLAGS)\
 $(SSP_CFLAGS)\
 -fopenmp
CORE_DPDC_LIBS =\
 $(AW_LIBS)\
 $(MPI_LIBS)\
 $(MRUMTL_LIBS)\
 $(RSYS_LIBS)\
 $(S3D_LIBS)\
 $(SCAM_LIBS)\
 $(SSF_LIBS)\
 $(SSP_LIBS)\
 -fopenmp\
 -lm

# Planeto
PLANETO_DPDC_CFLAGS=\
 $(RNATM_CFLAGS)\
 $(RNGRD_CFLAGS)\
 $(RSYS_CFLAGS)\
 $(S3D_CFLAGS)\
 $(SBUF_CFLAGS)\
 $(SCAM_CFLAGS)\
 $(SSF_CFLAGS)\
 $(SSP_CFLAGS)\
 $(SVX_CFLAGS)
PLANETO_DPDC_LIBS=\
 $(RNATM_LIBS)\
 $(RNGRD_LIBS)\
 $(RSYS_LIBS)\
 $(S3D_LIBS)\
 $(SBUF_LIBS)\
 $(SCAM_LIBS)\
 $(SSF_LIBS)\
 $(SSP_LIBS)\
 $(SVX_LIBS)\
-lm

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
