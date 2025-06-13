#include <filesystem>
#include <amdis/AMDiS.hpp>
#include <amdis/AdaptInstationary.hpp>
#include <amdis/ProblemStat.hpp>
#include <amdis/StandardProblemIteration.hpp>
#include <amdis/ProblemStatTraits.hpp>
#include <amdis/ElementVector.hpp>
#include <amdis/gridfunctions/ElementGridFunction.hpp>
#include <dune/grid/uggrid.hh>
#include <dune/alugrid/grid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <amdis/io/VTKWriter.hpp>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/vtk/datacollectors/discontinuousdatacollector.hh>
#include <dune/curvedgrid/grid.hh>
#include <dune/curvedgrid/gridfunctions/discretegridviewfunction.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/alugrid/common/backuprestore.hh>

#include <amdis/surfacebasis/surfacelagrangebasis.hh>
#include <amdis/surfacebasis/SurfaceToBulkInterpolator.hpp>
#include <amdis/surfacebasis/SurfaceDerivative.hpp>
#include <amdis/surfacebasis/InterpolationDataTransferSurface.hpp>
#include <amdis/ALEModelAmdis2/myIndices.hpp>
#include <amdis/ALEModelAmdis2/HeaderIncludes.hpp>

// the base problems for solving the system and computing the grid velocity
#include "Laplace.hpp"
#include "NavierStokesCahnHilliard.hpp"

// utility classes for model extensions
#include <amdis/ALEModelAmdis2/modelExtensions/ExtensionInterface.hpp>
#include <amdis/ALEModelAmdis2/modelExtensions/ExtensionManager.hpp>

// other classes
#include "ShellGeometry.hpp"
#include "src/TimeHandler.hpp"
#include <amdis/surfacegridfunctions/SurfaceDerivativeGridFunction.hpp>
#include <amdis/surfacegridfunctions/SurfaceDiscreteFunction.hpp>

//modelExtensions
#include <amdis/ALEModelAmdis2/modelExtensions/Activity.hpp>
#include <amdis/ALEModelAmdis2/modelExtensions/Binodal.hpp>

using namespace AMDiS;
using namespace Dune::Functions::BasisFactory;
using namespace Dune::Indices;
using namespace MyIndices;

/**
 * @brief The main method where the whole simulation happens
 * @param argc
 * @param argv
 * @return 0
 */
int main(int argc, char** argv)
{
    Environment env(argc, argv);

    // time variables
    double time = 0.0;
    double endTime = 1000;
    Parameters::get("adapt->start time", time);
    Parameters::get("adapt->end time", endTime);
    std::vector<double> timesteps;

    std::string path;
    Parameters::get("output directory", path);

    // copy source files to output folder for making the simulation reproducable
    copyFilesToOutputFolder();

    // read mesh file and store it into newMesh.msh, which is used in the code and updated during remeshing
    copyMeshFileToOutputFolder(time);

    // the grid type (ALUGrid)
    auto const dow = Grid::dimensionworld;

    // using a file reader to read the mesh
    MeshCreator<Grid> meshCreator("stokesMesh");

    // store the mesh into a host grid
    auto gridPtr = meshCreator.create();

    // perform global refinements of the mesh
    auto refine = Parameters::get<int>("stokesMesh->global refinements").value_or(0);
    gridPtr->globalRefine(refine);

    // define Partition-Indices of the mesh according to Physical Surface (2D) definition (e.g. from gmsh)
    std::vector<int> partitions = meshCreator.elementIds();
    std::vector<int> boundaryIds = meshCreator.boundaryIds();

    for (int i = 0; i < partitions.size(); ++i) partitions[i] = (partitions[i] > 0);

    // create a time handler to handle data in time
    GlobalBasis auxB(gridPtr->leafGridView(), power<2>(lagrange<1>(),flatInterleaved()));
    auto aux = std::make_shared<DOFVector<decltype(auxB)>>(auxB);
    TimeHandler<AdaptiveGrid<Grid>> timeHandler(time);

    // define facets vector for surface basis (facets = partitions on elements touching the interface, -1 else)
    double shellSize0 = 0.0;
    std::vector<int> facets = computeFacets(*gridPtr,partitions);

    // compute minimum and average element volumes (needed for refinement at the interface)
    timeHandler.estimatedElementSizes(gridPtr->leafGridView(),partitions,time);
    /// loop over remeshing steps (after remeshing, simulation restarts from here with the new mesh)
    while (time < endTime) {
        // if time > 0: update the grid to the new one and also partitions, facets, boundaryIds (, vMin)
        timeHandler.updateGridAfterRemeshing(gridPtr, partitions, facets, boundaryIds);

        // make a DOFVector filled with coordinates as GridFunction for the CurvedGrid
        GlobalBasis basis1(gridPtr->leafGridView(), power<2>(lagrange<1>(), flatInterleaved()));
        auto coordinates = std::make_shared<DOFVector<decltype(basis1)>>(basis1);
        valueOf(*coordinates) << X();

        // the moving mesh...changing coordsDF updates the mesh!
        Dune::CurvedGrid grid{*gridPtr, valueOf(*coordinates)};

        // for Cahn-Hilliard basis functions, determine the polynomial degree (P1 or P2) from parameter file
        auto deg = Parameters::get<int>("polynomial degree ch").value_or(2);
        auto activeWetting = Parameters::get<int>("extensions->activity").value_or(0);

        std::vector<int> partitionsEmpty(partitions.size());
        if (activeWetting) partitionsEmpty = partitions;
        // define pre basis factory for Navier-Stokes problem, including all solution components
        auto pbf = composite(power<2>(lagrange<2>(), flatInterleaved()),//velocity
                             surface(lagrange<1>(), facets),//pressure
                             surface(lagrange(deg), facets),//phi
                             surface(lagrange(deg), facets),//mu
                             power<2>(Dune::Surfacebasis::surfaceLagrange<1>(partitions), flatInterleaved()),//newCoords
                             power<2>(Dune::Surfacebasis::surfaceLagrange<1>(partitions, true),
                                      flatInterleaved()),//kappaVec
                             Dune::Surfacebasis::surfaceLagrange<1>(partitions, true),//lambda1 (principal stretches)
                             Dune::Surfacebasis::surfaceLagrange<1>(partitions, true),//lambda2
                             power<2>(Dune::Surfacebasis::surfaceLagrange<1>(partitions, true),
                                      flatInterleaved()),//force
                             power<2>(Dune::Surfacebasis::surfaceLagrange<1>(partitions, true), flatInterleaved()),//vS
                             Dune::Surfacebasis::surfaceLagrange<1>(partitions, true),//phiS
                             Dune::Surfacebasis::surfaceLagrange<1>(partitionsEmpty, true),//phiGamma
                             Dune::Surfacebasis::surfaceLagrange<1>(partitionsEmpty, true),//muGamma
                             flatLexicographic());

        // define problem for Navier-Stokes Cahn-Hilliard with basis on subdomains
        NavierStokesCahnHilliard nschProb("stokes", grid, pbf, boundaryIds, partitions);

        // define problem for Laplace with basis on subdomains, using the partitions-Vector
        Laplace<decltype(grid)> laplaceProb("laplace", nschProb.grid(), partitions, boundaryIds);

        AdaptInfo adaptInfo("adapt");
        adaptInfo.setTime(time);

        // fill in extensions of the model
        ExtensionManager extensionManager(nschProb);

        // set the problems into the timeHandler
        timeHandler.setProblems(nschProb, laplaceProb, extensionManager, adaptInfo);

        // initialize data
        timeHandler.initializeData();

        // the grid view and level 0 gridView of the nschProb
        auto const &gridView = nschProb.gridView();

        // some solution variables needed later
        auto v = nschProb.solution(_v);
        auto u = laplaceProb.solution();

        // kappaVec of the old time step but on the updated grid!
        DOFVector kappaVecOld(gridView, power<2>(Dune::Surfacebasis::surfaceLagrange<1>(partitions), flatInterleaved()),
                              tag::no_datatransfer{});

        // y value at start (for stretching force parameter lambda2, axisymmetric case)
        DOFVector y0(gridView, lagrange<1>());
        valueOf(y0) << X(1); // set y0 to y at start;

        // after remeshing, restore values from last time step before remeshing
        timeHandler.restoreDataAfterRemeshing(facets);

        // set intitial values for distances and lambda2 (for stretching force)
        timeHandler.updateForStretching(y0);

        // extract surface mesh from bulk mesh and print grid info
        auto const& surfacePreBasis = nschProb.globalBasis()->preBasis().subPreBasis(_phiS);
        auto const& surfaceGrid = surfacePreBasis.surfaceGrid();

        // print grid information
        timeHandler.printGridInfo(grid,surfaceGrid);

        // create shell geometry object
        // for interpolation between bulk and surface and sorting surface vertices for remeshing
        ShellGeometry shellGeometry{surfacePreBasis,gridView};

        // add internal boundary condition for laplace problem on the shell (incl. permeability)
        auto perm = Parameters::get<double>("parameters->permeability").value_or(0);
        auto normalMovement = Parameters::get<int>("move grid points with normal velocity").value_or(0);
        auto const& vGrid = nschProb.vGrid();

        auto v_corrected = makeGridFunction(
                invokeAtQP([perm, &vGrid, normalMovement](auto const& n,
                                                          auto const& v,
                                                          auto const& deltaP,
                                                          auto const& deltaPhi)
                                                          {
                    return (normalMovement ? v.dot(n)*n : v) - perm*(double)deltaP*std::abs((double)deltaPhi)*n + vGrid;
                    }, valueOf(*nschProb.normalDOF()), v, valueOf(*nschProb.deltaP()), valueOf(*nschProb.deltaPhi())),gridView);

        laplaceProb.problem().addConstraint(
                MeshMovementBC{*laplaceProb.problem().globalBasis(), partitions, v_corrected});

        // if no laplace solved -> no velocity in x direction for line mesh...todo: delete later
        auto ignoreLaplace = Parameters::get<int>("ignore laplace problem").value_or(0);
        if (ignoreLaplace) {
            nschProb.problem().addConstraint(
                    MeshMovementBC{*nschProb.problem().globalBasis(), _v, _v,
                                   partitions,
                                   makeGridFunction(invokeAtQP([]() { return FieldVector<double,2>{0.0,0.0}; }), gridView)
                    });
        }

        // initialize data of the nschProb
        timeHandler.initializeMainProb(shellSize0);

        // print surface grid info after initial refinement
        timeHandler.printSurfaceGridInfo(surfaceGrid);

        // compute shell and droplet volumes
        timeHandler.computeInitialVolumes();

        // compute the jump of pressure and phi (for needed permeability only)
        timeHandler.updateForPermeability();

        /// write initial files
        Dune::VTKWriter bulkWriter(gridView, Dune::VTK::nonconforming);
        writeFiles(nschProb,bulkWriter,kappaVecOld,y0,timesteps,adaptInfo,path,extensionManager,1);

        timeHandler.countAfterRemeshing(0);

        /// now the actual time step loop starts (everything before was preparation)
        DOFVector coords(gridView,power<2>(lagrange(1),flatInterleaved()));
        while (adaptInfo.time() < adaptInfo.endTime()) {
            Dune::Timer timestepTime;

            //depending on parameter in initfile, change timestep size in the first few time steps and/or after remeshing
            timeHandler.adaptTimestepIfNecessary();

            // update the time by adding the time step size
            timeHandler.updateTime();

            // print shell area, velocity and volume
            timeHandler.printShellInfo();

            // iteration of main problem and laplace problem with all necessary data updates
            timeHandler.iteration(kappaVecOld,y0);

            // move the grid with solution of Laplace as displacement
            // actual grid update: the values of coords are set so the coodinates DOFVector
            timeHandler.moveGrid(coords,coordinates);

            // update (surface) basis after grid movement, to get correct surfaceGrid
            timeHandler.updateBasis();
            shellGeometry.updateGridFunctions();

            // print surface grid information
            timeHandler.printSurfaceGridInfo(surfaceGrid);

            /// write the files
            auto every = Parameters::get<int>("every").value_or(1);
            if (timeHandler.count() % every == 0)
                writeFiles(nschProb,bulkWriter,kappaVecOld,y0,timesteps,adaptInfo,path,extensionManager);

            // end of timestep reached, print info
            timeHandler.printEndOfTimestep(timestepTime);

            /// if necessary and desired, correct surface interface grid points if they have moved due to spurious currents...
            auto doCorrect = Parameters::get<int>("correct point distances").value_or(0);
            if (doCorrect)
               correctSurfaceInterfacePointCoords(nschProb,shellGeometry,coords,coordinates,u,adaptInfo.timestep());

            /// if necessary, perform remeshing and leave the current loop to restart the simulation with
            /// the backup files from the remeshing function
            if (remeshing(nschProb, gridPtr, shellGeometry, partitions,
                          facets, boundaryIds, adaptInfo, path, timeHandler.vMin())) {
                timeHandler.noRemeshingHappenedBefore(0);
                break;
            }
        }
    }
}

// include ccache with:  -DCMAKE_CXX_COMPILER_LAUNCHER=ccache