#!/bin/bash

# INFO: Installs alberta3 (for dune-only), gmsh, dune and amdis2 modules based on a toolchain and dune release
# AUTHOR: M. Mokbel (with help by L. Wittwer and Simon Praeterius)

### installing necessary software

sudo apt-get install \
libalberta-dev \
libmetis-dev \
libopenmpi-dev \
libparmetis-dev \
libsuitesparse-dev



DUNE_RELEASE="2.9"

## PUT YOUR DIRECTORY HERE ##

DUNE_CORE_MODULES="dune-common dune-geometry dune-grid dune-istl dune-localfunctions"
DUNE_STAGING_MODULES="dune-functions dune-typetree dune-uggrid"
DUNE_EXTENSIONS_MODULES="dune-alugrid dune-foamgrid dune-vtk dune-grid-glue"
DUNE_BACKENDS="ISTL"

#########################################################################
# YOU SHOULD NOT HAVE TO CHANGE ANYTHING BELOW
#########################################################################


# Paths, variables etc.


GMSH_INSTALL_DIR="gmsh"
ALBERTA_INSTALL_DIR="alberta3"
DUNE_INSTALL_DIR="dune"
AMDIS_INSTALL_DIR="dune"

mkdir -p ${GMSH_INSTALL_DIR}
mkdir -p ${ALBERTA_INSTALL_DIR}
mkdir -p ${DUNE_INSTALL_DIR}


#####################################################
# Installing gmsh
#####################################################
GMSH_GIT="https://gitlab.onelab.info/gmsh/gmsh.git"
GMSH_BUILD_DIR=${GMSH_INSTALL_DIR}/gmsh/build

cd ${GMSH_INSTALL_DIR}
git clone ${GMSH_GIT}
cd gmsh

mkdir build
cd build
cmake -DENABLE_BUILD_DYNAMIC=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
make -j4
sudo make install

cd ../../../

#####################################################
# Installing Alberta3, release/3.0.3 (for dune-only)
#####################################################
# NOTE: Could be optimised to have only one installation for each toolchain (and not dune release)

ALBERTA_GIT="https://gitlab.com/alberta-fem/alberta3.git"
ALBERTA_BUILD_DIR=${ALBERTA_INSTALL_DIR}/build

cd ${ALBERTA_INSTALL_DIR}

# clonig git repo and checkout release 3.0.3
git clone ${ALBERTA_GIT} .

git checkout releases/3.0.3

mkdir build

# config and install
./generate-alberta-automakefiles.sh
autoreconf --force --install

./configure --prefix=${ALBERTA_BUILD_DIR} --disable-fem-toolbox

make -j4
make install

cd ..

# PKG_CONFIG_PATH is not correctly detected for the dune installation, lets set it
export PKG_CONFIG_PATH="${PKG_CONFIG_PATH}:${ALBERTA_BUILD_DIR}/lib/pkgconfig"


#############################################
# Create opts-Files
#############################################
mkdir ${DUNE_INSTALL_DIR}/opts
OPTS_FILES_DIR=${DUNE_INSTALL_DIR}/opts
OPTS_FILES=""

## change compiler options to match your compiler (find your compiler with "whereis clang") ##
for be in $DUNE_BACKENDS; do
	RELEASE_OPTS_TEXT=$(cat <<-EOM
    CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=RelWithDebInfo
                 -DCMAKE_CXX_COMPILER=/usr/bin/clang++
                 -DCMAKE_C_COMPILER=/usr/bin/clang
                 -DCMAKE_CXX_FLAGS=-w
                 -DDUNE_ENABLE_PYTHONBINDINGS:BOOL=0
                 -DDUNE_ENABLE_PYTHONBINDINGS=OFF"
    MAKE_FLAGS="-j2"
    BUILDDIR="build-clang"
EOM
)
	OPTS_FILE_NAME=clang.opts
	echo "${RELEASE_OPTS_TEXT}" > ${OPTS_FILES_DIR}/${OPTS_FILE_NAME}
	OPTS_FILES="${OPTS_FILES} ${OPTS_FILE_NAME}"
done


#############################################
# Installing DUNE and AMDiS2
#############################################

# downloading dune modules
cd ${DUNE_INSTALL_DIR}

for module in $DUNE_CORE_MODULES; do
	git clone https://gitlab.dune-project.org/core/${module}.git
done

for module in $DUNE_STAGING_MODULES; do
	git clone https://gitlab.dune-project.org/staging/${module}.git
done

for module in $DUNE_EXTENSIONS_MODULES; do
	git clone https://gitlab.dune-project.org/extensions/${module}.git
done

git clone https://gitlab.mn.tu-dresden.de/iwr/dune-curvedgrid.git
git clone https://gitlab.mn.tu-dresden.de/iwr/dune-curvedgeometry.git

PATH=dune-common/bin:$PATH
export DUNE_ROOT=.

# downloading amdis + amdis-extensions
cd ../${AMDIS_INSTALL_DIR}

git clone https://gitlab.com/amdis/amdis.git
git clone https://gitlab.mn.tu-dresden.de/amdis/amdis-extensions.git
git clone https://gitlab.com/aland-lab/amdis-surfacebasis.git

# set the DUNE_CONTROL_PATH to include the amdis dir, exporting it s.t. dunecontrol finds it
export DUNE_CONTROL_PATH="./amdis:."


# checking out the correct release
#dune-common/bin/dunecontrol exec git checkout releases/2.9
: <<'END'
cd ../dune/dune-common
git checkout 1c46e0ddd00f9dedb66c163fd4918b0dbebcae87
cd ../dune-geometry
git checkout 8ea54dbb585d17bad235917c4bd546eb6317721a
cd ../dune-grid
git checkout 38117baece23d2f83f6bf6569802a40a06bbf476
cd ../dune-istl
git checkout 38a50dc30f476aca6b1b4ced8524366549ba5c86
cd ../dune-localfunctions
git checkout 08e08658bdc982917418f11799a1d00b07ff1aaa


cd ../dune-functions
git checkout 904ec38e924c415e6f231bdfd60769d9306d806e
cd ../dune-typetree
git checkout 8423973ef7cecaa474e8f12b8b0c4df7fba7a888
cd ../dune-uggrid
git checkout 0d43ffd8d782185f8a28074a4a8c01fe60fa8380

cd ../dune-alugrid
git checkout 6b0d597391d4690e98a1b05f919cf46173e73b88
cd ../dune-vtk
git checkout e6152c937303b0ff385de2b3f9e285c3ca7bb168
cd ../dune-foamgrid
git checkout 25d1bf102273bd92f21c88c992522c5de738d93e
cd ../dune-grid-glue
git checkout 68dbc53c55e8a993578d4e449f58542e4951203e

cd ../dune-curvedgrid
git checkout b6b7d42ebc344fae5b4a48d7aa82b53c2fd8342f
cd ../dune-curvedgeometry
git checkout 6f3d2b0a654921b4764dcd48aa74908667001d84

cd ../amdis
git checkout c9efc4ee88cdf1640fad13e00b9157e0d7678aa5
cd ../amdis-extensions
git checkout 3a7af0200f868f880c18e6e8c25a7c58cf5aca3a
cd ../amdis-surfacebasis
git checkout 237fc7b94baaa0c3f01e41114980c0d61648ba5c
cd ../..
git clone https://gitlab.com/aland-lab/marcelmokbel/alemodelamdis2.git
cd alemodelamdis2
git checkout 08fcaf9bd39d7bc679aed567be4c4c18ec9bb69b
END
# replace file in dune-grid-glue
cd ..
cp alemodelamdis2/codim0extractor.hh dune/dune-grid-glue/dune/grid-glue/extractors/codim0extractor.hh

# building the modules for all the opts files
cd dune
for opt in $OPTS_FILES; do
	dune-common/bin/dunecontrol --opts=../${OPTS_FILES_DIR}/${opt} all
done

cd ../alemodelamdis2
git checkout origin/hivProject
for opt in $OPTS_FILES; do
  export DUNE_CONTROL_PATH=../dune
	../dune/dune-common/bin/dunecontrol --opts=../${OPTS_FILES_DIR}/${opt} --current all
done

echo "Installation done!"