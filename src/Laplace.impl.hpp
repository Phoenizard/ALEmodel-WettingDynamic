/**
 * Implementation of the Laplace problem
 * Therefore, we use the AMDiS::BaseProblem
 * You'll need the amdis-extensions module for that
 **/

#pragma once

#include <amdis/AMDiS.hpp>
#include <amdis/ProblemStat.hpp>
#include <amdis/ProblemStatTraits.hpp>
#include <amdis/LocalOperators.hpp>
#include <amdis/Integrate.hpp>
#include <amdis/localoperators/StokesOperator.hpp>
#include "amdis/ALEModelAmdis2/MeshMovementBC.hpp"
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

using namespace AMDiS;
using namespace Dune::Functions::BasisFactory;
using namespace Dune::Indices;

template <class Grid>
Laplace<Grid>::Laplace(std::string const& name, std::shared_ptr<AdaptiveGrid> grid, std::vector<int> partitions_, std::vector<int> boundaryIDs_)
: Super(name, grid), partitions(partitions_), boundaryIDs(boundaryIDs_)
{
   // set the correct boundary IDs to the mesh
    problem().boundaryManager()->setBoundaryIds(boundaryIDs);
}

//set initial values or reserve memory for variables if existent
template <class Grid>
void Laplace<Grid>::initData(AdaptInfo& adaptInfo) {
}

// do something at the beginning of each time step
template <class Grid>
void Laplace<Grid>::initTimestep(AdaptInfo& adaptInfo)
{
}

// define the boundary conditions of the problem
template <class Grid>
void Laplace<Grid>::fillBoundaryConditions(AdaptInfo& adaptInfo) {
    double xmin=100000,ymin=100000,xmax=-100000,ymax=-100000;
    for (auto const& el : elements(this->gridView())) {
        for (auto const& is : intersections(this->gridView(),el)) {
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

    auto axi = Parameters::get<int>("axisymmetric").value_or(0);
    problem().addDirichletBC(left,0,0,0.0);
    problem().addDirichletBC(right,0,0,0.0);
    problem().addDirichletBC(left,1,1,0.0);
    problem().addDirichletBC(right,1,1,0.0);
    problem().addDirichletBC(top,1,1,0.0);

    if (!axi) problem().addDirichletBC(bottom,0,0,0.0);
    problem().addDirichletBC(bottom,1,1,0.0);

    int line = Parameters::get<int>("line mesh").value_or(0);
    int noMovementOnBoundariesLine = Parameters::get<int>("line mesh no movement on top and bottom").value_or(0);
    if (noMovementOnBoundariesLine && line) {
        problem().addDirichletBC(top,0,0,0.0);
        problem().addDirichletBC(bottom,0,0,0.0);
    }
    // set BCs for hivProject
    auto use_hiv = Parameters::get<int>("use hivProject").value_or(0);
    if (use_hiv) {
        problem().addDirichletBC(33,0,0,0.0);
        problem().addDirichletBC(33,1,1,0.0);
        problem().addDirichletBC(44,0,0,0.0);
        problem().addDirichletBC(44,1,1,0.0);

        problem().addDirichletBC(top,0,0,0.0);
    }
}

// implement the operators of the equations that you want to solve
template <class Grid>
void Laplace<Grid>::fillOperators(AdaptInfo& adaptInfo) {
    for (int i = 0; i < 2; ++i)
       problem().addMatrixOperator(sot(1.0), i, i);
}

