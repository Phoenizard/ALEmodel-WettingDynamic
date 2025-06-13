#pragma once

#include "amdis/ALEModelAmdis2/localoperators/SecondOrderGradTestGradTrialAxi.hpp"
#include "amdis/ALEModelAmdis2/myIndices.hpp"

using namespace AMDiS;
using namespace Dune::Functions::BasisFactory;
using namespace Dune::Indices;
using namespace MyIndices;

using Grid = Dune::ALUGrid<2,2,Dune::simplex,Dune::ALUGridRefinementType::conforming>;

/** \brief interpolation of a grid function onto a Container (e.g. DiscreteGridViewFunction)
 * @param container - the gridFunction to be interpolated on
 * @param gridFct - the source gridFunction
 */
template <class Container, class GridFunction>
void interpolate(Container& container, GridFunction const& gridFct)
{
    // create a temporary container
    auto tmp = container.coefficients();

    // interpolate the gridFunction
    auto&& gf = makeGridFunction(gridFct, container.basis().gridView());
    Dune::Functions::interpolate(container.basis(), tmp, gf);

    // store the result of the interpolation in the container
    container.coefficients() = std::move(tmp);
}

/** \brief compute DOF corrdinates as a DiscreteGridViewFunction
 *
 * @param gridView - the gridView of the bulk grid
 * @return
 */
template <class GridView>
auto vertexCoordinates(GridView const& gridView)
{
    // create a special (discrete) gridFunction with internal storage of coefficients
    Dune::DiscreteGridViewFunction coordinates{gridView, 1};

    // interpolate the coordinate function into this gridFunction container
    Dune::Functions::interpolate(coordinates.basis(), coordinates.coefficients(), [](auto const& x) { return x; });

    return coordinates;
}

/** \brief interpolate normal as ElementVector to the surface grid
 *
 * @param normalBulk - a vector of element vectors on the bulk grid, containing the components of the normal on \Gamma
 * @param normal0 - an ElementVector on the surface containing the x-component of the normal, to be computed
 * @param normal1 - an ElementVector on the surface containing the y-component of the normal, to be computed
 * @param pb - the surface basis, containing the surface grid and index maps to interpolate between surface and bulk
 */
template <class EVV, class EVS, class SPB>
void interpolateEVToSurface(EVV const& normalBulk, EVS& normal0, EVS& normal1, SPB const& pb) {
    auto const& map = pb.surfaceToFluidMap();
    auto surfaceGridView = pb.surfaceGridView();
    auto gridView = normalBulk[0].gridView();

    auto const& isS = surfaceGridView.indexSet();
    auto const& isB = gridView.indexSet();

    // the map maps between bulk and surface elements
    for (const auto& sbmap : map) {
        auto const& es = surfaceGridView.grid().entity(sbmap.first);
        auto const& eb = gridView.grid().entity(sbmap.second);

        normal0.data()[isS.index(es)] = normalBulk[0].data()[isB.index(eb)];
        normal1.data()[isS.index(es)] = normalBulk[1].data()[isB.index(eb)];
    }
}

/** compute current kappaVec with finite elements
 *
 * @param kappaVecOld - the curvature vector of the current surface, to be computed here
 * @param surfaceGrid - the current surface grid
 * @param surfacePreBasis - the surface basis instance, containing index maps for interpolation from bulk to surface
 * @param normal - an elementVector containing the element-wise normal vectors for elements with a facet on \Gamma
 */
template <class DV, class SG, class SPB, class EV>
void computeKappaVecOld(DV& kappaVecOld,
                        std::shared_ptr<SG> const& surfaceGrid,
                        SPB const &surfacePreBasis,
                        EV const& normal) {
    auto axi = Parameters::get<int>("axisymmetric").value_or(0);
    auto preBasis = composite(power<2>(lagrange<1>(),flatInterleaved()), //coords
                              power<2>(lagrange<1>(),flatInterleaved()),flatLexicographic());//kappaVecOld

    ProblemStat prob("kappaVecOld", *surfaceGrid, preBasis);
    auto k0 = Parameters::get<double>("spontaneous curvature").value_or(0);
    prob.initialize(INIT_ALL);

    // x =
    auto opX = makeOperator(tag::testvec_trialvec{},1.0);
    prob.addMatrixOperator(opX,_0,_0);

    // = oldCoords
    auto opV = makeOperator(tag::testvec{},X(),4);
    prob.addVectorOperator(opV,_0);

    // L =
    auto opL = makeOperator(tag::testvec_trialvec{},(X(1)*axi + 1.0-axi));
    prob.addMatrixOperator(opL,_1,_1);

    // subtract spontaneous curvature if k0 nonzero
    // normal needs to be transferred to a surface ElementVector first
    ElementVector normalXS(*prob.grid(),0.0);
    ElementVector normalYS(*prob.grid(),0.0);

    if (k0) {
       interpolateEVToSurface(normal,normalXS,normalYS,surfacePreBasis);

        // then subtract spontaneous curvature term
       auto opSpontaneous1 = makeOperator(tag::testvec{},
                                           invokeAtQP([k0, axi](auto const& nx, auto const& ny, auto const& y) {
                                                          auto n = FieldVector<double,2>{nx,ny};
                                                          return -k0 * (y*axi + 1.0-axi) * n;
                                                      },
                                                      valueOf(normalXS),valueOf(normalYS), X(1)), 2);
    }
    // <grad x, grad v>
    if (axi) {
        for (int i = 0; i < surfaceGrid->dimensionworld; i++) {
            auto opLB = makeOperator(tag::gradtest_gradtrial_axi{}, X(1));
            prob.addMatrixOperator(opLB, makeTreePath(_1, i), makeTreePath(_0, i));
        }
    } else {
        for (int i = 0; i < surfaceGrid->dimensionworld; i++) {
            auto opLB = makeOperator(tag::gradtest_gradtrial{}, 1.0);
            prob.addMatrixOperator(opLB, makeTreePath(_1, i), makeTreePath(_0, i));
        }
    }

    if (axi) {
        auto opAxi4 = makeOperator(tag::test_trial{},1.0/(X(1)),4);
        prob.addMatrixOperator(opAxi4,makeTreePath(_1,1),makeTreePath(_0,1));


        auto predicate = [](auto const& x){ return x[1] < 1.e-9; }; // define boundary
        int use;
        Parameters::get("axisymmetric dirichlet kappaVecOld",use);
        if (use) prob.addDirichletBC(predicate, makeTreePath(_1,1), makeTreePath(_1,1), 0.0);
    }

    int line = Parameters::get<int>("line mesh").value_or(0);
    double ymin = 100000, ymax = -100000;
    if (line) {
        //find corner coordinates of the grid
        for (auto const &el: elements(prob.gridView())) {
                for (int i = 0; i < el.geometry().corners(); i++) {
                    ymin = std::min(ymin, el.geometry().corner(i)[1]);
                    ymax = std::max(ymax, el.geometry().corner(i)[1]);
                }
            }
        auto predicate = [ymin,ymax](auto const& x){ return x[1] < ymin + 1.e-5 || x[1] > ymax - 1.e-5; }; // define boundary
        auto zeroTopBottom = Parameters::get<int>("line mesh no movement on top and bottom").value_or(0);
        if (zeroTopBottom) prob.addDirichletBC(predicate, makeTreePath(_1,0), makeTreePath(_1,0), 0.0);
        prob.addDirichletBC(predicate, makeTreePath(_1,1), makeTreePath(_1,1), 0.0);
    }

    AdaptInfo adaptInfo("adapt");
    prob.assemble(adaptInfo);
    prob.solve(adaptInfo);

    DOFVector kappaVecSurf(prob.gridView(), power<SG::dimensionworld>(lagrange<1>(),flatInterleaved()));
    valueOf(kappaVecSurf) << prob.solution(_1);
    interpolateFromSurfaceGrid(kappaVecSurf,kappaVecOld,surfacePreBasis.surfaceToFluidMap(),1);

    prob.writeFiles(adaptInfo);
}

/** compute the phase field values on the surfaceGrid (only implemented for one-sided phase field)
 *
 * @param phiDGVF - a DiscreteGridViewFunction containing the phase field
 * @param nschProb - the main problem
 */
template <class DGVFS, class Prob>
void phiSurface(DGVFS& phiDGVF, Prob& nschProb) {
    auto const& spb = nschProb.globalBasis()->preBasis().subPreBasis(_phiGamma);
    auto const& partitions = nschProb.partitions();

    Dune::DiscreteGridViewFunction<decltype(nschProb.gridView()), 1> phiBulkDGVF{nschProb.gridView(),1};
    DOFVector phiSurf(nschProb.gridView(), lagrange<1>());
    auto outside = Parameters::get<int>("phase field outside").value_or(0);
    auto inside = Parameters::get<int>("phase field inside").value_or(0);

    auto phi = nschProb.solution(_phi);
    if (outside) interpolateOneSided(phiSurf,phi,0,partitions,_0);
    if (inside) interpolateOneSided(phiSurf,phi,1,partitions,_0);
    Dune::Functions::interpolate(phiBulkDGVF.basis(), phiBulkDGVF.coefficients(), valueOf(phiSurf));
    auto map = spb.surfaceToFluidMapIdx();
    for (auto it : map) {
        phiDGVF.coefficients()[it.first] = phiBulkDGVF.coefficients()[it.second];
    }
}

static double shellSize00 = 0.0;
/** compute the facets vector based on a given partitions vector, the result is the facet number whithin the element for
 * the facet belonging to \Gamma
 * @param grid - the bulk grid
 * @param partitions - a vector containing elementIds (1 inside the shell, 0 outside)
 * @param doComputeArea - flag to decide, whether the area of the shell needs to be computed
 * @param size0 - the area of the shell (initial state)
 * @return
 */
template <class Grid>
std::vector<int> computeFacets(Grid const& grid,
                               std::vector<int> const& partitions,
                               int doComputeArea = 0,
                               double& size0 = shellSize00) {
    int closedShell = Parameters::get<int>("closed shell").value_or(1); // 1 if closed shell, 0 else
    std::vector<int> facets(partitions.size(),-1);
    double shellSize0 = 0.0;
    if (!closedShell) {
        for(const auto& e : elements(grid.levelGridView(0))) {
            auto const& indexSet = grid.levelGridView(0).indexSet();
            for (const auto& is : intersections(grid.levelGridView(0),e)) {
                if (!is.neighbor())
                    continue;

                auto p = partitions[indexSet.index(e)];
                auto q = partitions[indexSet.index(is.outside())];
                if (p == 0 && q == 1) {
                    facets[indexSet.index(e)] = 1.0;
                    facets[indexSet.index(is.outside())] = 0.0;
                }
            }
        }
    } else {
        facets = partitions;
        if (doComputeArea) {
            // compute initial shell area
            for (const auto &e: elements(grid.levelGridView(0))) {
                auto const &indexSet = grid.levelGridView(0).indexSet();
                for (const auto &is: intersections(grid.levelGridView(0), e)) {
                    if (!is.neighbor())
                        continue;

                    auto p = partitions[indexSet.index(e)];
                    auto q = partitions[indexSet.index(is.outside())];
                    if (p == 0 && q == 1) {
                        size0 += is.geometry().volume();
                    }
                }
            }
            std::cout << "initial shell area = " << size0 << "\n";
        }
    }
    return facets;
}

/** \brief this function copies the current source files to the desired output folder
 * this way, the simulations can later be reproduced easier
 */

void copyFilesToOutputFolder() {
    // create directory for all the output files
    std::string path;
    Parameters::get("output directory", path);
    filesystem::create_directories(path);
    filesystem::create_directories(path + "/backup");
    filesystem::create_directories(path + "/src");

    // copy parameter file to output folder
    std::string parameterFileName, parameterOutputFileName;
    Parameters::get("parameter file name", parameterFileName);
    Parameters::get("parameter output file name", parameterOutputFileName);
    std::filesystem::copy_file(parameterFileName,
                               parameterOutputFileName,
                               std::filesystem::copy_options::overwrite_existing);

    std::filesystem::path src = "./src";
    for (const auto& entry : std::filesystem::recursive_directory_iterator("./src")) {
        if (!(entry.is_directory() || entry.path().string().find("/CMakeFiles/") != std::string::npos)) {
            std::filesystem::path dest = path + "/src";
            std::filesystem::copy(entry.path(),
                                  dest / entry.path().lexically_relative(src),
                                  std::filesystem::copy_options::overwrite_existing);
        }
    }
    std::filesystem::copy("./src",
                          path + "/src",
                          std::filesystem::copy_options::overwrite_existing);
    std::cout << "Copying done!" << std::endl;
    std::filesystem::copy("./amdis",
                          path + "/amdis",
                          std::filesystem::copy_options::recursive |
                          std::filesystem::copy_options::overwrite_existing);
    std::cout << "Copying done!" << std::endl;
}

/** this function copes the mesh file to the output folder, to make it usable later and make reproduction
 * of the simulation easier
 * @param time - the current time of the simulation
 */
void copyMeshFileToOutputFolder(double time) {
    if (time < 1.e-10) {
        std::string inputMesh, outputMesh;
        Parameters::get("stokesMesh->macro file name", inputMesh);
        Parameters::get("newMesh->macro file name", outputMesh);
        std::ifstream src(inputMesh, std::ios::binary);
        std::ofstream dst(outputMesh, std::ios::binary);
        dst << src.rdbuf();
    }
}