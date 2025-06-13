#pragma once

using namespace AMDiS;
using namespace Dune::Functions::BasisFactory;
using namespace Dune::Indices;

/// compute the minimum angle on the mesh
template<typename GridView>
double minAngle(GridView const& gridView) {
    //compute minimum angle between two grid points on the grid
    double minangle = 1000;
    for (auto const &el: elements(gridView)) {
        for (int i = 0; i < el.geometry().corners(); i++) {
            auto v1 = el.geometry().corner(i);
            auto v2 = el.geometry().corner((i + 1) % el.geometry().corners());
            auto v3 = el.geometry().corner((i + 2) % el.geometry().corners());
            auto e1 = v2 - v1;
            auto e2 = v3 - v1;
            double angle = acos(dot(e1, e2) / (e1.two_norm() * e2.two_norm()));
            minangle = std::min(minangle, angle);
        }
    }
    std::cout << "Minimum angle on the mesh = " << minangle << "\n";
    return minangle;
}

template <class DV>
void setLambda2EqLambda1Axi(DV const& lambda1, DV& lambda2, std::vector<int> const& partitions) {
    auto const &indexSet = lambda1.basis().gridView().indexSet();
    auto const &indexSet0 = lambda1.basis().gridView().grid().levelGridView(0).indexSet();
    auto localView = lambda2.basis().localView();

    for (auto const &e: elements(lambda1.basis().gridView())) {
        for (const auto &is: intersections(lambda1.basis().gridView(), e)) {
            if (!is.neighbor())
                continue;

            auto father = e;
            while (father.hasFather())
                father = father.father();

            auto fatherOut = is.outside();
            while (fatherOut.hasFather())
                fatherOut = fatherOut.father();

            auto p = partitions[indexSet0.index(father)];
            auto q = partitions[indexSet0.index(fatherOut)];
            if (p != q) {
                localView.bind(e);
                auto const &node = localView.tree();

                for (std::size_t j = 0; j < e.geometry().corners(); ++j) {
                    auto x = e.geometry().corner(j);
                    if (x[1] < 1.e-10) {
                        auto val = lambda1.coefficients().get(localView.index(node.localIndex(j)));
                        lambda2.coefficients().set(localView.index(node.localIndex(j)), val);
                    }
                }
            }
        }
    }
    lambda2.coefficients().finish();
}

template <class GF>
double findVal(GF const& l1S, FieldVector<double,2> const& x) {
    Dune::DiscreteGridViewFunction coords{l1S.basis().gridView(), 1};
    Dune::Functions::interpolate(coords.basis(), coords.coefficients(),
                                 [](auto const &x) { return x; });
    for (int i = 0; i < coords.coefficients().size(); ++i) {
        if ((x - coords.coefficients()[i]).two_norm() < 1.e-7)
            return l1S.coefficients()[i];
    }

    // if point did not exist on old grid, find face, where the point is located in and interpolate values
    auto localView = l1S.basis().localView();

    for (auto const &e: elements(l1S.basis().gridView())) {
        localView.bind(e);
        auto node = localView.tree().child(0);
        auto left = e.geometry().corner(0);
        auto right = e.geometry().corner(1);
        auto face = right - left;
        auto xFace = x - left;
        auto t1 = xFace[0] / face[0];
        auto t2 = xFace[1] / face[1];
        if (std::abs(t1 - t2) < 1.e-7) { // x lies on this face
            auto bulk_idxLeft = localView.index(node.localIndex(0))[0];
            auto bulk_idxRight = localView.index(node.localIndex(1))[0];
            return ((1.0 - t1) * l1S.coefficients()[bulk_idxLeft] + t1 * l1S.coefficients()[bulk_idxRight]);
        }
    }
    // nothing found, which means that a point from the old grid has been removed but a new one has been created by gmsh
    // hence: set value to value of nearest point on old grid
    double dist = 100000;
    int idx = -1;
    for (int i = 0; i < coords.coefficients().size(); ++i) {
        if ((x - coords.coefficients()[i]).two_norm() < dist) {
            dist = (x - coords.coefficients()[i]).two_norm();
            idx = i;
        }
    }
    return l1S.coefficients()[idx];
}

template <class P, class GF, class Grid, class EV>
void interpolateToNewSurface(P const& nschProb,
                             GF& functionNew,
                             Grid const& gridPtrNew,
                             std::vector<int>& partitionsOld,
                             std::vector<int>& partitions,
                             EV const& gfOld,
                             double val = 0,
                             int axi = 0,
                             int isLambda2 = 0) {
    DOFVector lambda1S(nschProb.gridView(),
                       power<1>(Dune::Surfacebasis::surfaceLagrange<1>(partitionsOld), flatInterleaved()));
    DOFVector lambda1SNew(gridPtrNew.leafGridView(),
                          power<1>(Dune::Surfacebasis::surfaceLagrange<1>(partitions), flatInterleaved()));

    // interpolate to surfaceLagrange DOFVector
    valueOf(lambda1S).interpolate(valueOf(gfOld), tag::average{});

    // interpolate then to the surface grid to do the interpolation to new grid there
    auto const &spb = lambda1S.basis().preBasis().subPreBasis();
    auto l1S = surfaceGridFct(lambda1S, spb, _0);

    auto const &spbNew = lambda1SNew.basis().preBasis().subPreBasis();
    auto l1SNew = surfaceGridFct(lambda1SNew, spbNew, _0);

    Dune::DiscreteGridViewFunction coords{l1SNew.basis().gridView(), 1};
    Dune::Functions::interpolate(coords.basis(), coords.coefficients(),
                                 [](auto const &x) { return x; });

    // the actual interpolation step
    for (int i = 0; i < coords.coefficients().size(); ++i) {
        auto x = coords.coefficients()[i];

        auto val = findVal(l1S, x);
        l1SNew.coefficients()[i] = val;

        if (axi && isLambda2) {
            auto val2 = findVal(l1S, x);
            l1SNew.coefficients()[i] = val2;
        }
    }

    // now bring the data back top the bulk, first to a DiscreteGridViewFunction, then to the DOFVector
    using GV = decltype(gridPtrNew.leafGridView());
    Dune::DiscreteGridViewFunction<GV,1> aux1{gridPtrNew.leafGridView(), 1};
    for (auto it: spbNew.surfaceToFluidMapIdx()) {
        aux1.coefficients()[it.second] = l1SNew.coefficients()[it.first];
    }
    valueOf(lambda1SNew) << makeGridFunction(aux1,gridPtrNew.leafGridView()) - val;

    // then make the DOFVectors readable for outputs, by setting them to 1.0
    // everywhere and interpolate surface values afterward
    surfaceValueOf(functionNew,partitions) << valueOf(lambda1SNew, _0);
    valueOf(functionNew) += val;
}

template <class Problem, class CGrid, class SG>
bool remeshing(Problem& nschProb,
               std::unique_ptr<CGrid>& gridPtr,
               SG& shellGeometry,
               std::vector<int>& partitions,
               std::vector<int>& facets,
               std::vector<int>& boundaryIds,
               AdaptInfo& adaptInfo,
               std::string const& path,
               double const& vMin) {
    auto const& distances = nschProb.distances();
    auto const& lambda2 = nschProb.lambda2();

    // perform remeshing
    using Grid = Dune::ALUGrid<2,2,Dune::simplex,Dune::ALUGridRefinementType::conforming>;

    int remeshType = 0;
    int everyRemesh = 1000;
    double minRemeshAngle;
    Parameters::get("remeshing type", remeshType);
    Parameters::get("remeshing every i-th time step", everyRemesh);
    Parameters::get("remeshing at mimimum angle", minRemeshAngle);

    auto const &gridView = nschProb.gridView();
    if (remeshType == 0 && minAngle(gridView) < minRemeshAngle ||
        (remeshType == 1 && adaptInfo.timestepNumber()%everyRemesh == 0)) {

        // define phi as DiscreteGridViewFunction on the surface
        Dune::DiscreteGridViewFunction<decltype(shellGeometry.gridView()), 1> phiDGVF{shellGeometry.gridView(),1};
        phiSurface(phiDGVF, nschProb);

        // create new mesh out of surface coordinates
        generateMesh("newMesh", path, nschProb, adaptInfo);

        //store the mesh into the host grid
        MeshCreator<Grid> meshCreatorNew("newMesh");
        auto gridPtrNew = meshCreatorNew.create();
        auto partitionsNew = meshCreatorNew.elementIds();
        auto boundaryIdsNew = meshCreatorNew.boundaryIds();
        auto partitionsOld = partitions;
        partitions = std::move(partitionsNew);
        boundaryIds = std::move(boundaryIdsNew);
        facets.clear();
        facets = computeFacets(*gridPtr,partitions);

        // create a gridFunction with coordinates of the new grid and make CurvedGrid object
        // coordinates DOFVector for remeshing
        GlobalBasis basis2(gridPtrNew->leafGridView(), power<2>(lagrange<1>(),flatInterleaved()));
        auto coordinatesNewGrid = std::make_shared<DOFVector<decltype(basis2)>>(basis2);
        valueOf(*coordinatesNewGrid) << X();
        auto gridNew = Dune::CurvedGrid{*gridPtrNew, valueOf(*coordinatesNewGrid)};

        using Gr = decltype(gridNew);
        auto deg = Parameters::get<int>("polynomial degree ch").value_or(2);

        auto pbfIpol = composite(power<2>(lagrange<2>(),flatInterleaved()),//velocity
                                 lagrange<1>(),//lagrange multiplier for v
                                 surface(lagrange(deg),facets),//phi
                                 surface(lagrange(deg),facets),flatLexicographic());//mu

        int ref_int  = Parameters::get<int>("refinement->interface").value_or(10);
        for (int i = 0; i <= ref_int; ++i) {
            // create InterpolateGrids Object for the saddle point problem
            InterpolateGrids<Gr,decltype(pbfIpol)> ipolProb("ipolProb", gridNew, *nschProb.grid()->hostGrid(), pbfIpol, vMin);
            ipolProb.initialize(INIT_ALL);

            // initialize the problem
            ipolProb.initBaseProblem(adaptInfo);
            ipolProb.initTimeInterface();

            ipolProb.beginIteration(adaptInfo);

            // assemble the matrix
            auto markFlag = ipolProb.problem().markElements(adaptInfo);
            ipolProb.problem().buildAfterAdapt(adaptInfo, markFlag, true, true);
            //ipolProb.oneIteration(adaptInfo);

            auto &&oldVelocity = nschProb.solution(_v);
            auto &&oldPhi = nschProb.solution(_phi);
            auto &&oldMu = nschProb.solution(_mu);

            //lambda1 on the new grid
            GlobalBasis basis1(ipolProb.problem().gridView(), lagrange<1>());
            auto oldLamb = makeDOFVector(basis1);

            GlobalBasis basisSubdomains(ipolProb.problem().gridView(), power<1>(surface(lagrange(deg),facets),flatInterleaved()));

            auto inside = Parameters::get<int>("phase field inside").value_or(0);
            auto outside = Parameters::get<int>("phase field outside").value_or(0);
            // assemble the RHS of the system for velocity on the new grid
            assembleRHSVec(ipolProb.glue_,*ipolProb.problem().rhsVector(), oldVelocity, _v);
            if (outside) {
                assembleRHSAvg(ipolProb.glue_,*ipolProb.problem().rhsVector(), oldPhi, _phi, partitions, partitionsOld);
                assembleRHSAvg(ipolProb.glue_,*ipolProb.problem().rhsVector(), oldMu, _mu, partitions, partitionsOld);
            }
            if (inside) {
                assembleRHSAvg(ipolProb.glue_,*ipolProb.problem().rhsVector(), oldPhi, _phi, partitions, partitionsOld,1);
                assembleRHSAvg(ipolProb.glue_,*ipolProb.problem().rhsVector(), oldMu, _mu, partitions, partitionsOld,1);
            }

            // solve the problem
            ipolProb.problem().solve(adaptInfo, true, false);
            ipolProb.problem().estimate(adaptInfo);
            ipolProb.closeTimestep(adaptInfo);
            ipolProb.endIteration(adaptInfo);

            // if points on surface have been removed, correct the phase field near the shell accordingly
            if (Parameters::get<int>("remove interface points").value_or(0))
                correctPhaseShell(ipolProb.solution(_phi),phiDGVF,_phi,partitions,partitionsOld);

            if (i <= ref_int) {
                ipolProb.problem().markElements(adaptInfo);
                ipolProb.problem().adaptGrid(adaptInfo);
            }

            ipolProb.writeFiles(adaptInfo);

            // store oldVelocity into grid function
            GlobalBasis basis2D(ipolProb.problem().gridView(), power<2>(lagrange<2>(),flatInterleaved()));
            auto oldV = makeDOFVector(basis2D);
            auto oldPh = makeDOFVector(basisSubdomains);
            auto oldM = makeDOFVector(basisSubdomains);

            valueOf(oldV) << ipolProb.solution(_v);
            valueOf(oldPh) << ipolProb.solution(_phi);
            valueOf(oldM) << ipolProb.solution(_mu);

            oldV.backup(path + "/backup/uhOld" + std::to_string(adaptInfo.time()));
            oldPh.backup(path + "/backup/phiOld" + std::to_string(adaptInfo.time()));
            oldM.backup(path + "/backup/muOld" + std::to_string(adaptInfo.time()));
        }

        // interpolate lambda1 + lambda2 + phiGamma + muGamma to new grid
        {
            auto incompressibleShell = Parameters::get<int>("is shell incompressible").value_or(0);
            auto axi = Parameters::get<int>("axisymmetric").value_or(0);
            auto phiGamma0 = Parameters::get<int>("initial phiGamma").value_or(0.5);
            GlobalBasis basisNew(gridPtrNew->leafGridView(), lagrange<1>());
            auto lambda1New = makeDOFVector(basisNew);
            auto lambda2New = makeDOFVector(basisNew);
            auto phiGammaNew = makeDOFVector(basisNew);
            auto muGammaNew = makeDOFVector(basisNew);

            auto phiGamma = nschProb.solution(_phiGamma);
            auto muGamma = nschProb.solution(_muGamma);

            interpolateToNewSurface(nschProb,phiGammaNew,*gridPtrNew,partitionsOld,partitions,phiGamma,phiGamma0,axi,0);
            interpolateToNewSurface(nschProb,muGammaNew,*gridPtrNew,partitionsOld,partitions,muGamma,0.0,axi,0);
            if (!incompressibleShell) {
                interpolateToNewSurface(nschProb,lambda1New,*gridPtrNew,partitionsOld,partitions,distances,1.0,axi,0);
                interpolateToNewSurface(nschProb,lambda2New,*gridPtrNew,partitionsOld,partitions,*lambda2,1.0,axi,1);

                // correct lambda2 on the symmetry axis
                if (axi) {
                    setLambda2EqLambda1Axi(lambda1New, lambda2New, partitions);
                }
            } else {
                valueOf(lambda1New) << 1.0;
                valueOf(lambda2New) << 1.0;
            }
            lambda1New.backup(path + "/backup/lambda1Old" + std::to_string(adaptInfo.time()));
            lambda2New.backup(path + "/backup/lambda2Old" + std::to_string(adaptInfo.time()));
            phiGammaNew.backup(path + "/backup/phiGammaOld" + std::to_string(adaptInfo.time()));
            muGammaNew.backup(path + "/backup/muGammaOld" + std::to_string(adaptInfo.time()));
        }

        //update grid to new one
        gridPtr.swap(gridPtrNew);
        std::string filename(path + "/backup/gridRefined" + std::to_string(adaptInfo.time()));
        std::ofstream file( filename.c_str());
        gridPtr->hostGrid()->backup(file,ALUGrid::MacroFileHeader::Format::ascii);

        // leave the current loop and start over with the new one on the new grid
        return true;
    }
    return false;
}