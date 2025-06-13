#pragma once

#include <any>

// operators
#include <amdis/LocalOperators.hpp>
#include <amdis/localoperators/StokesOperator.hpp>
#include <amdis/surfacebasis/SurfaceLocalOperators.hpp>

#include <amdis/GridFunctionOperator.hpp>
#include <amdis/ALEModelAmdis2/localoperators/PointIntegralOperatorAxi.hpp>
#include <amdis/ALEModelAmdis2/localoperators/PointIntegralOperatorAxiKappa.hpp>
#include <amdis/ALEModelAmdis2/localoperators/SurfaceZeroOrderTestVecNormal.hpp>
#include <amdis/ALEModelAmdis2/localoperators/SurfaceZeroOrderTestVecTrialSpecial.hpp>
#include <amdis/ALEModelAmdis2/localoperators/SurfaceFirstOrderTestDivTrialvecAxi.hpp>
#include <amdis/ALEModelAmdis2/localoperators/SurfaceZeroOrderTestTrialAxi2.hpp>
#include <amdis/ALEModelAmdis2/localoperators/SurfaceFShellOperator.hpp>
#include <amdis/ALEModelAmdis2/localoperators/CahnHilliardOperator.hpp>
#include <amdis/ALEModelAmdis2/localoperators/NavierStokesOperator.hpp>
#include <amdis/ALEModelAmdis2/localoperators/Surface_xS_kappaVecOperator.hpp>
#include <amdis/ALEModelAmdis2/localoperators/Surface_lambda1_lambda2Operator.hpp>
#include <amdis/ALEModelAmdis2/localoperators/Surface_couplingOperator.hpp>
#include <amdis/ALEModelAmdis2/localoperators/SurfaceCahnHilliardOperator.hpp>

// interpolators, helper functions etc.
#include <amdis/ALEModelAmdis2/OneSidedInterpolator.hpp>
#include "src/helperFunctions.hh"
#include <amdis/ALEModelAmdis2/CorrectInterfaceMembraneCoords.hpp>
#include "src/jumpAcrossGamma.hpp"
#include "src/updateFunctions.hpp"
#include <amdis/ALEModelAmdis2/hostGridInterpolator.hpp>
#include <amdis/ALEModelAmdis2/writeFiles.hh>
#include "amdis/ALEModelAmdis2/remeshing/InterpolateGrids.hpp"
#include "amdis/ALEModelAmdis2/remeshing/MeshGenerator.h"
#include "src/remeshingLoop.hh"
#include <amdis/ALEModelAmdis2/AxisymmetricDirichletBC.hpp>


