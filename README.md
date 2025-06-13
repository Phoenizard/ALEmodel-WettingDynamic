This project is mainly cloned from ALE Lab Project:

https://gitlab.com/aland-lab/marcelmokbel/alemodelamdis2

The key modification is changing the initial setting and the condition in NavierStokes 
CahnHiliard Equation.

---

This code is designed to perform ALE (Arbitrary Lagrangian Eulerian) simulations with Dune-AMDiS.
It also contains a model for wetting of biomembranes, using a phase field approach for two phase flow, which is used in the papers

**[1] A simualtion method for the wetting dynamics of liquid droplets on deformable membranes** by M. Mokbel, D. Mokbel, S. Liese, C. Weber, and S. Aland, submitted to _SIAM Journal on Scientific Computing_

and 

**[2] ESCRT condensed phases facilitate formation of multivesicular bodies by endosomal membrane 
bending and scission**
 by Y. Wang, S. Li, M. Mokbel, A. May, Z. Liang, Y. Zeng, W. Wang, H- Zhang, F. Yu, K. Sporbeck, L. Jiang, S. Aland, J. Agudo-Canalejo, R. Knorr, and X. Fang, submitted to _Nature_

Operating Systems and Versions
========================
This code has to be installed on a Linux System from 2020 or newer (e.g. Ubuntu20.04). For Windows, it is highly recommended to use the WSL
(Windows Subsystem Linux). Mac systems are currently not supported.

You can install the code and all dependencies using the provided install script (see the section
_Automatic installation_ below), of you can install it manually using the description in the section _Preparing
the Sources_
. The installation should be done within a maximum of one working day (including resolving any issues during installation).

The code has been tested on Ubuntu systems 20.04, 22.04. 

Automatic installation
=========================

Download the file `installScript.sh`, copy it into your local folder. Check your clang version and location
(typically something like `/usr/bin/clang++` - it can be found by typing `clang --version` and `whereis clang`
in a terminal). Compare
to `lines 105 and 106` of `installAScript.sh` and change paths to your clang version.

Then, in your local folder, run
```
sh installScript.sh
```
This will automatically install all dune and AMDiS modules and it will also download and install the code for
this project to the folder `alemodelamdis2`.  
So, after successful installation, you can just type
```
cd alemodelamdis2
build-clang/src/ALEModelAmdis2 init/stokes.dat.2d
```

to run your code. And if you want to re-compile your code, just type (replacing `path/to/dune` and `/path/to/amdis` with your local
directories)
```
export DUNE_CONTROL_PATH=/path/to/amdis:/path/to/dune  
/path/to/dune/dune-common/bin/dunecontrol --current all 
```
Preparing the Sources
=========================

If you want to install manually, You will need to install the DUNE and AMDIS libraries. Information about AMDiS (and how to install DUNE and AMDiS) can be found on  
`https://amdis.readthedocs.io/en/latest/`

If you do not use the automatic install script, you'll need the
following programs installed on your system:

cmake >= 3.1, Alberta, metis, openmpi, parmetis, suitesparse. In Linux, run the following, to get them

````
sudo apt-get install \
libalberta-dev \
libmetis-dev \
libopenmpi-dev \
libparmetis-dev \
libsuitesparse-dev
````  
You'll also need the following dune modules:

- amdis,
- dune-common,
- dune-istl,
- dune-geometry,
- dune-uggrid,
- dune-grid,
- dune-localfunctions,
- dune-typetree,
- dune-functions,
- dune-alugrid,
- dune-foamgrid,
- dune-vtk,
- dune-curvedgeometry,
- dune-curvedgrid,
- dune-grid-glue
- amdis-extensions
- amdis-surfacebasis

Therefore, do the following:
````
mkdir DUNE_DIR  
cd DUNE_DIR 
git clone https://gitlab.com/amdis/amdis.git
git clone https://gitlab.dune-project.org/core/dune-common.git  
git clone https://gitlab.dune-project.org/core/dune-istl.git  
git clone https://gitlab.dune-project.org/core/dune-geometry.git  
git clone https://gitlab.dune-project.org/staging/dune-uggrid.git  
git clone https://gitlab.dune-project.org/core/dune-grid.git  
git clone https://gitlab.dune-project.org/core/dune-localfunctions.git  
git clone https://gitlab.dune-project.org/staging/dune-typetree.git    
git clone https://gitlab.dune-project.org/staging/dune-functions.git  
git clone https://gitlab.dune-project.org/extensions/dune-alugrid.git  
git clone https://gitlab.dune-project.org/extensions/dune-foamgrid.git  
git clone https://gitlab.dune-project.org/extensions/dune-grid-glue.git  
git clone https://gitlab.dune-project.org/extensions/dune-vtk.git   
git clone https://gitlab.mn.tu-dresden.de/iwr/dune-curvedgeometry.git  
git clone https://gitlab.mn.tu-dresden.de/iwr/dune-curvedgrid.git  
git clone https://gitlab.mn.tu-dresden.de/amdis/amdis-extensions.git  
git clone https://gitlab.com/aland-lab/amdis-surfacebasis.git

dune-common/bin/dunecontrol all
````

Getting started
---------------

If these preliminaries are met, you should clone this repository to your local folder and compile the code
```
git clone https://gitlab.com/aland-lab/marcelmokbel/alemodelamdis2.git
cd alemodelamdis2  
(git checkout _myBranch_)
export DUNE_CONTROL_PATH=path/to/amdis:path/to/dune 
path/to/dune/dune-common/bin/dunecontrol --current all  
```
Note that the git checkout _myBranch_ is optional, if you want to checkout a specific branch. You have to replace 
_myBranch_ with your branch name (e.g. "hivProject" or "activeDropletsVsVesicles")

To run your code, you can now execute
```
build-cmake/src/ALEModelAmdis2 init/stokes.dat.2d
```

Note that if dune is not installed properly you will either
have to add the directory where the dunecontrol script resides (probably
./dune-common/bin) to your path or specify the relative path of the script.

Most probably you'll have to provide additional information to dunecontrol
(e.g. compilers, configure options) and/or make options.

The most convenient way is to use options files in this case. The files
define four variables:

CMAKE_FLAGS      flags passed to cmake (during configure)

An example options file might look like this:
```
CMAKE_FLAGS="-DCMAKE_CXX_COMPILER=clang++ 
-DCMAKE_C_COMPILER=clang
-DDUNE_ENABLE_PYTHONBINDINGS:BOOL=0"
MAKE_FLAGS="-j2"
BUILDDIR="build-clang"
```
If you save this information into example.opts you can pass the opts file to
dunecontrol via the --opts option, e.g.
```
dunecontrol --opts=example.opts all
```
Running the code
---------

To run your code, execute
```
build-cmake/src/ALEModelAmdis2 init/stokes.dat.2d
```

_The following only holds for the branch vesicleBudding and is an explanation of how to use the parameters for the paper [2] mentioned above!_

The file `init/stokes.dat.2d` contains the parameters used in the simulation. You can change these to your needs.
An additional parameter file `init/smallVesicles.2d` is provided for convenience, using the parameters for smaller vesicles (MVB simulations).
The parameter file(s) contain(s) several parameters that are relevant to perform simulations for GUVs (Giant Unilamellar Vesicles, typically 10µm of size) and MVBs (Multi-Vesicular Bodies, typically 1µm of size), especially
the bending stiffness "Kb", area dilation modulus "Ka", lineTension, and the membrane and surface tensions. Here
`parameters->sigma` refers to the condensate surface tension and the parameters
`parameters->sigma0` and `parameters->sigma1` refer to the ambient membrane tension and the condensate membrane tension,
respectively. Note that the parameters are nondimensionalized by length scale `L` and time scale `T`.

A large number of other parameters can be changed within the parameter file. You can try them out.

An example is provided already when using the `init/stokes.dat.2d` parameter file. There, a GUV simulation is performed
for the first 10 time steps. This should take less than 5 minutes to run on a 'normal' computer. A longer simulation can be tested by increasing the parameter `adapt->end time`. 

The results
of the simualtion are written into a folder, that is created separately, which is specified by
`output directory`. In this case it is `../ResultsALEModelAmdis2/test`. There, you will find two files called
`surface.pvd` and `bulk.pvd`. These can be used to visualize the results, e.g. with ParaView (www.paraview.org). In
ParaView, load the two files and choose the resulting field that you want to see, e.g. `phi` in the bulk and
`phiS` on the surface to see the values of the phase field. These indicate, where the condensate (`phi=1`) and the ambient
`phi = 0` are located.

You should see 10 time steps there. The simulation is axisymmetric, hence you will only see half of the cross section of the membrane and droplet. The membrane is initially prolate shaped, with the droplet being a half sphere attatched on the rigth of the membrane. During those first 10 time steps, the membrane should be deformed slightly by the droplet, especially near the three phase contact region. 

More info on the dunecontrol installation script
---------

See

     dunecontrol --help

for further options.


The full build system is described in the dune-common/doc/buildsystem (Git version) or under share/doc/dune-common/buildsystem if you installed DUNE!


