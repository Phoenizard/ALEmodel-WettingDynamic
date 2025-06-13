/**
 * Implementation of the Grdi Interpolation problem for remeshing
 * Therefore, we use the AMDiS::BaseProblem
 * You'll need the amdis-extensions module for that
 **/

#pragma once

#include <amdis/AMDiS.hpp>
#include <amdis/ProblemStat.hpp>
#include <amdis/ProblemStatTraits.hpp>
#include <amdis/LocalOperators.hpp>
#include <amdis/Integrate.hpp>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

using namespace AMDiS;
using namespace Dune::Functions::BasisFactory;
using namespace Dune::Indices;

// Evaluate DiscreteFunction in global coordinates according to
// https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
template <class DF, class Domain>
double evalAtPoint(DF const& f, Domain const& x, std::vector<int> const& partitions, int side=0)
{
    auto const& gv = f.basis().gridView();
    auto const& gv0 = f.basis().gridView().grid().levelGridView(0);

    for (auto const& e : elements(gv)) { //search for all elements of base grid
        auto father = e;
        while (father.hasFather())
            father = father.father();

        int p = partitions[gv0.indexSet().index(father)];
        if (p != side) //only look at elements of the desired subdomain
            continue;

        //corners: (only for triangle grids)
        auto v0 = e.geometry().corner(0);
        auto v1 = e.geometry().corner(1);
        auto v2 = e.geometry().corner(2);
        auto area = 0.5 *(-v1[1]*v2[0] + v0[1]*(-v1[0] + v2[0]) + v0[0]*(v1[1] - v2[1]) + v1[0]*v2[1]);

        auto s = 1/(2*area)*(v0[1]*v2[0] - v0[0]*v2[1] + (v2[1] - v0[1])*x[0] + (v0[0] - v2[0])*x[1]);
        auto t = 1/(2*area)*(v0[0]*v1[1] - v0[1]*v1[0] + (v0[1] - v1[1])*x[0] + (v1[0] - v0[0])*x[1]);

        if (s > -1.e-10 && t > -1.e-10 && 1-s-t > 0) { //this is the desired element
            auto geometry = e.geometry();
            auto localFct = localFunction(f);
            localFct.bind(e);
            return localFct(geometry.local(x));
        }
    }
    return 0.0;
}

// assemble the rhs vector with the term (v_old,psi) where v_old lives on the old grid
// and psi is the test function on the new grid (glue is GridGlue object that returns
// intersection iterator, with inside: element on old grid, outside: element on new grid
template <class Glue, class RHS, class GF, class TP>
void assembleRHSVec(Glue& glue, RHS& rhsVector, GF& rhsFct, TP treePath) {
    msg("ipolProb::assembleRHSVec()");

    // loop over all intersections
    auto localView = rhsVector.basis().localView();
    int const dimension = GF::EntitySet::GlobalCoordinate::dimension;
    auto lf = localFunction(rhsFct);
    rhsVector.init(rhsVector.basis(),false);
    for (const auto& intersection : intersections(glue)) {
        // bind local rhsFct to the element on the old grid
        lf.bind(intersection.inside());

        // bind localView to element on the new grid
        auto&& el = intersection.outside();
        localView.bind(el);

        // create node, local finite element and intersectionVector
        auto const &node = localView.tree().child(treePath);
        auto const &lfe = node.child(0).finiteElement();
        std::vector<double> intersectionVector(localView.size());

        // loop over quadrature points of the intersection
        int order = 2; //*(dimension*lfe.localBasis().order()-1);
        const auto &quadRule = Dune::QuadratureRules<double, dimension>::rule(intersection.type(), order);
        for (const auto &qp: quadRule) {
            // The multiplicative factor in the integral transformation formula
            double ds = intersection.geometry().integrationElement(qp.position()) * qp.weight();

            // Position of the current quadrature point in the current element (new grid)
            auto &&localOut = intersection.geometryInOutside().global(qp.position());
            // Position of the current quadrature point in the current element (old grid)
            auto &&localIn = intersection.geometryInInside().global(qp.position());

            // Values of the local basis and the RHS in the quadrature Points
            const auto rhsValues = lf(localIn);
            std::vector<FieldVector<double, 1> > lagrangeValues;
            lfe.localBasis().evaluateFunction(localOut, lagrangeValues);

            for (std::size_t i = 0; i < lagrangeValues.size(); ++i) {
                for (std::size_t k = 0; k < node.degree(); ++k) {
                    auto localIndex = node.child(k).localIndex(i);
                    auto value = ds * rhsValues[k] * lagrangeValues[i];
                    intersectionVector[localIndex] += value;
                }
            }
        }
        //fill rhsVector with calculated Values stored in intersectionVector
        rhsVector.scatter(localView, intersectionVector, Assigner::plus_assign{});
        lf.unbind();
    }
    rhsVector.finish();
}

// assemble the rhs vector with the term (phi_old,psi) where phi_old lives on the old grid
// and psi is the test function on the new grid (glue is GridGlue object that returns
// intersection iterator, with inside: element on old grid, outside: element on new grid
// partitions and partitions old arte vectors to determin from which subdomain an element is taken
// only evaluate the values on the subdomain with index side
template <class Glue, class Coeff, class GF, class TP>
void assembleRHSAvg(Glue& glue, Coeff&& rhsVector, GF const& rhsFct, TP treePath, std::vector<int> const& partitions, std::vector<int> const& partitionsOld,
                    int side = 0) {
    msg("ipolProb::assembleRHSVec()");

    // loop over all intersections
    auto localView = rhsVector.basis().localView();
    int const dimension = GF::EntitySet::GlobalCoordinate::dimension;
    auto lf = localFunction(rhsFct);
    rhsVector.init(rhsVector.basis(),false);
    for (const auto& intersection : intersections(glue)) {
        auto fatherIn = intersection.inside();
        while (fatherIn.hasFather())
            fatherIn = fatherIn.father();

        auto fatherOut = intersection.outside();
        while (fatherOut.hasFather())
            fatherOut = fatherOut.father();

        auto p = partitions[glue.template gridView<1>().grid().levelGridView(0).indexSet().index(fatherOut)];
        auto q = partitionsOld[glue.template gridView<0>().grid().levelGridView(0).indexSet().index(fatherIn)];

        if (p != side || q != side)
            continue;

        // bind local rhsFct to the element on the old grid
        lf.bind(intersection.inside());

        // bind localView to element on the new grid
        auto&& el = intersection.outside();
        localView.bind(el);

        // create node, local finite element and intersectionVector
        auto const &node = localView.tree().child(treePath);
        auto const &lfe = node.finiteElement();
        std::vector<double> intersectionVector(localView.size());

        // loop over quadrature points of the intersection
        int order = 2; //*(dimension*lfe.localBasis().order()-1);
        const auto &quadRule = Dune::QuadratureRules<double, dimension>::rule(intersection.type(), order);
        for (const auto &qp: quadRule) {
            // The multiplicative factor in the integral transformation formula
            double ds = intersection.geometry().integrationElement(qp.position()) * qp.weight();

            // Position of the current quadrature point in the current element (new grid)
            auto &&localOut = intersection.geometryInOutside().global(qp.position());
            // Position of the current quadrature point in the current element (old grid)
            auto &&localIn = intersection.geometryInInside().global(qp.position());

            // Values of the local basis and the RHS in the quadrature Points
            const auto rhsValues = lf(localIn);
            std::vector<FieldVector<double, 1> > lagrangeValues;
            lfe.localBasis().evaluateFunction(localOut, lagrangeValues);

            for (std::size_t i = 0; i < lagrangeValues.size(); ++i) {
                auto localIndex = node.localIndex(i);
                auto value = ds * rhsValues * lagrangeValues[i];
                intersectionVector[localIndex] += value;
            }
        }
        //fill rhsVector with calculated Values stored in intersectionVector
        rhsVector.scatter(localView, intersectionVector, Assigner::plus_assign{});
        lf.unbind();
    }
    rhsVector.finish();
}

template <class GF, class DGVF, class TP>
void correctPhaseShell(GF&& gf, DGVF const& phiDGVF, TP treePath, std::vector<int> const& partitions, std::vector<int> const& partitionsOld) {
        //correct phase field at the shell if points were removed
    auto removedPoints = Parameters::get<double>("distance for removal of interface points").value_or(0.0);
    auto const& indexSet0 = gf.basis().gridView().grid().levelGridView(0).indexSet();
    gf.coefficients().init(gf.basis(),false);
    auto localView = gf.basis().localView();

    if (removedPoints > 1.e-7) {
        Dune::DiscreteGridViewFunction coordsDGVF{phiDGVF.basis().gridView(), 1};
        Dune::Functions::interpolate(coordsDGVF.basis(), coordsDGVF.coefficients(), [](auto const &x) { return x; });

        for (const auto &el: elements(gf.basis().gridView())) {
            localView.bind(el);
            auto father = el;
            while (father.hasFather())
                father = father.father();
            for (const auto &is: intersections(gf.basis().gridView(), el)) {
                if (!is.neighbor())
                    continue;
                auto fatherOut = is.outside();
                while (fatherOut.hasFather())
                    fatherOut = fatherOut.father();

                int p = partitions[indexSet0.index(father)];
                int q = partitions[indexSet0.index(fatherOut)];
                if (p == 0 && q == 1) {
                    int subEntityCodim = 1;
                    int subEntityIndex = is.indexInInside();

                    auto const &node = localView.tree().child(treePath);
                    auto const &lfe = node.finiteElement();
                    auto const &lc = lfe.localCoefficients();
                    auto re = Dune::referenceElement(el);

                    for (std::size_t i = 0, k = 0; i < lc.size() && k < 2; ++i) {
                        auto localKey = lc.localKey(i);
                        if (re.subEntities(subEntityIndex, subEntityCodim, localKey.codim()).contains(
                                localKey.subEntity())) {
                            FieldVector<double, 2> const &coord = is.geometry().corner(i);
                            for (int j = 0; j < phiDGVF.coefficients().size(); j++) {
                                if ((coord - coordsDGVF.coefficients()[j]).two_norm() < 1.e-5) {
                                    gf.coefficients().set(localView.index(node.localIndex(i)), phiDGVF.coefficients()[j]);
                                }
                            }
                        }
                    }
                }
            }
        }
        gf.coefficients().finish();
    }
}


template <class Element>
auto truePredicate = [](const Element&, unsigned int) { return true; };

template <class Grid, class PBF>
InterpolateGrids<Grid,PBF>::InterpolateGrids(std::string const& name, Grid& grid, Grid& gridOld, PBF& pbf, double const& vMin)
        : ProblemIterationInterface()
        , ProblemTimeInterface()
        , ProblemInstatBase(name)
        , problem_(name + "->space", grid, pbf)
        , name_(name)
        , grid_(&grid)
        , gridOld_(&gridOld)
        , glue_(std::make_shared<Extractor>(gridOld.leafGridView(), truePredicate<Element>),
                std::make_shared<Extractor>(grid.leafGridView(), truePredicate<Element>),
                std::make_shared<OverlappingMerge>())
        , elementSizes(grid,-1)
        , elementLevels(grid,-1)
        , vMin_(vMin)
{
    glue_.build();
    Dune::GridGlue::GridGlueVtkWriter::write(glue_, "mergedGrids");
}


// define the boundary conditions of the problem
template <class Grid, class PBF>
void InterpolateGrids<Grid,PBF>::fillBoundaryConditions(AdaptInfo& adaptInfo) {

    //find corner coordinates of the grid
    double xmin=100000,ymin=100000,xmax=-100000,ymax=-100000;
    for (auto const& el : elements(grid_->leafGridView())) {
        for (auto const& is : intersections(grid_->leafGridView(),el)) {
            if (is.neighbor()) //only consider boundary intersections
                continue;
            for (int i = 0; i<is.geometry().corners(); i++) {
                xmin = std::min(xmin,is.geometry().corner(i)[0]);
                ymin = std::min(ymin,is.geometry().corner(i)[1]);
                xmax = std::max(xmax,is.geometry().corner(i)[0]);
                ymax = std::max(ymax,is.geometry().corner(i)[1]);
            }
        }
    }
    FieldVector<double,2> zero(0);

    auto bottom = [ymin](auto const& x){ return x[1] <= ymin + 1e-7; }; // define boundary
    auto left = [xmin](auto const& x){ return x[0] <= xmin + 1e-7; }; // define boundary
    auto right = [xmax](auto const& x){ return x[0] >= xmax - 1e-7; }; // define boundary
    auto top = [ymax](auto const& x){ return x[1] >= ymax - 1e-7; }; // define boundary

    // same BCs as for velocity on old grid
    auto axi = Parameters::get<int>("axisymmetric").value_or(0);
    int line = Parameters::get<int>("line mesh").value_or(0);

    if (axi || line) {
        problem().addDirichletBC(bottom, makeTreePath(_v,1), makeTreePath(_v,1), 0.0);
    } else {
        problem().addDirichletBC(bottom, _v, _v, zero);
    }

    int noMovementOnBoundariesLine = Parameters::get<int>("line mesh no movement on top and bottom").value_or(0);
    if (noMovementOnBoundariesLine && line) {
        problem().addDirichletBC(top, _v, _v, zero);
        problem().addDirichletBC(bottom, _v, _v, zero);
    }

    if(line) {
        problem().addDirichletBC(top, makeTreePath(_v,1), makeTreePath(_v,1), 0.0);
    } else {
        problem().addDirichletBC(top, _v, _v, zero);
    }

    auto zeroLeftRight = Parameters::get<int>("no movement on left and right").value_or(0);
    if (zeroLeftRight) {
        problem().addDirichletBC(left, _v, _v, zero);
        problem().addDirichletBC(right, _v, _v, zero);
    }
    problem().addDirichletBC(left, _lambda1, _lambda1, 0.0);

    problem().addDirichletBC(33, _v, _v, zero);
    problem().addDirichletBC(44, _v, _v, zero);
}

// implement the operators of the equations that you want to solve
template <class Grid, class PBF>
void InterpolateGrids<Grid, PBF>::fillOperators(AdaptInfo& adaptInfo) {

    // <w_i,v_i>
    auto op0 = makeOperator(tag::testvec_trialvec{}, 1.0);
    problem().addMatrixOperator(op0,_v,_v);

    // <d_i(w_i), \lambda>
    auto opP = makeOperator(tag::divtestvec_trial{}, 1.0);
    problem().addMatrixOperator(opP, _v, _lambda1);

    // <q, d_i(u_i)>
    auto opDiv = makeOperator(tag::test_divtrialvec{}, 1.0);
    problem().addMatrixOperator(opDiv, _lambda1, _v);

    // <w_i,v_i>
    auto opPhi0 = makeOperator(tag::test_trial{}, 1.0);
    problem().addMatrixOperator(opPhi0,_phi,_phi);

    // <w_i,v_i>
    auto opMu0 = makeOperator(tag::test_trial{}, 1.0);
    problem().addMatrixOperator(opMu0,_mu,_mu);
}

template <class Grid, class PBF>
void InterpolateGrids<Grid, PBF>::initData(AdaptInfo& adaptInfo) {
    int ref_bulk = Parameters::get<int>("refinement->bulk").value_or(2);
    int ref_int  = Parameters::get<int>("refinement->interface").value_or(10);

    elementSizes.resizeZero();
    elementLevels.resizeZero();
    for (auto const& e : elements(problem_.gridView())) {
        elementSizes.data()[problem_.gridView().indexSet().index(e)] = e.geometry().volume();
        elementLevels.data()[problem_.gridView().indexSet().index(e)] = e.level();
    }
    auto vMin = vMin_;
    auto marker = GridFunctionMarker("b-interface", problem().grid(),
                                     invokeAtQP([ref_int, ref_bulk, vMin](double phi, double size, int level) -> int {
                                         return phi > 0.05 && phi < 0.95 && size > vMin ? ref_int :
                                                (phi > 0.05 && phi < 0.95 && size <= vMin ? std::min(level,ref_int) : ref_bulk);
                                     }, this->getPhase(), valueOf(elementSizes), valueOf(elementLevels)));
    problem().addMarker(Dune::wrap_or_move(std::move(marker)));
}


template <class Grid, class PBF>
void InterpolateGrids<Grid, PBF>::initTimestep(AdaptInfo& adaptInfo) {
    glue_.build();
}

template <class Grid, class PBF>
void InterpolateGrids<Grid, PBF>::closeTimestep(AdaptInfo& adaptInfo) {
    elementSizes.resizeZero();
    elementLevels.resizeZero();
    for (auto const& e : elements(problem_.gridView())) {
        elementSizes.data()[problem_.gridView().indexSet().index(e)] = e.geometry().volume();
        elementLevels.data()[problem_.gridView().indexSet().index(e)] = e.level();
    }
    //writeFiles(adaptInfo);
}
template <class Grid, class PBF>
void InterpolateGrids<Grid, PBF>::solveInitialProblem(AdaptInfo &adaptInfo) {
    transferInitialSolution(adaptInfo);
}