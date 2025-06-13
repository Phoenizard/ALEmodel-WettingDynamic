#ifndef BREASTCANCER_MYINDICES_H
#define BREASTCANCER_MYINDICES_H

// subdomain basis for pressure and phase field spaces
#include <amdis/ALEModelAmdis2/subdomains/subdomainbasis.hh>
#include <amdis/ALEModelAmdis2/subdomains/surfacebasis.hh>

using namespace AMDiS;
using namespace Dune::Functions::BasisFactory;
using namespace Dune::Indices;// the basis used in the main problem

using Grid = Dune::ALUGrid<2,2,Dune::simplex,Dune::ALUGridRefinementType::conforming>;
using AdGr = AdaptiveGrid<Grid>;
using GridView = typename AdGr::LeafGridView;
using MyPBF = decltype(power<2>(lagrange<1>(),flatInterleaved()));
using MyBasis = decltype(GlobalBasis{std::declval<GridView>(),std::declval<MyPBF>()});

using GF = decltype(valueOf(std::declval<DOFVector<MyBasis>&>()));

int deg = Parameters::get<int>("polynomial degree ch").value_or(2);
using PBF = decltype(composite(power<2>(lagrange<2>(),flatInterleaved()),//velocity
                               surface(lagrange<1>(), std::declval<std::vector<int>>()),//pressure
                               surface(lagrange(deg), std::declval<std::vector<int>>()),//phi
                               surface(lagrange(deg), std::declval<std::vector<int>>()),//mu
                               power<2>(Dune::Surfacebasis::surfaceLagrange<1>(std::declval<std::vector<int>>()),flatInterleaved()),//newCoords
                               power<2>(Dune::Surfacebasis::surfaceLagrange<1>(std::declval<std::vector<int>>(),true),flatInterleaved()),//kappaVec
                               Dune::Surfacebasis::surfaceLagrange<1>(std::declval<std::vector<int>>(),true),//lambda1 (principal stretches)
                               Dune::Surfacebasis::surfaceLagrange<1>(std::declval<std::vector<int>>(),true),//lambda2
                               power<2>(Dune::Surfacebasis::surfaceLagrange<1>(std::declval<std::vector<int>>(),true),flatInterleaved()),//force
                               power<2>(Dune::Surfacebasis::surfaceLagrange<1>(std::declval<std::vector<int>>(),true),flatInterleaved()),//vS
                               Dune::Surfacebasis::surfaceLagrange<1>(std::declval<std::vector<int>>(),true),//phiS
                               Dune::Surfacebasis::surfaceLagrange<1>(std::declval<std::vector<int>>(),true),//phiGamma
                               Dune::Surfacebasis::surfaceLagrange<1>(std::declval<std::vector<int>>(),true),//muGamma
                               flatLexicographic()));
using CurvedGrid = Dune::CurvedGrid<AdGr,GF>;

namespace MyIndices {
    static constexpr auto _v = Dune::Indices::_0; //velocity
    static constexpr auto _p = Dune::Indices::_1; //pressure
    static constexpr auto _phi = Dune::Indices::_2; //bulk phase field
    static constexpr auto _mu = Dune::Indices::_3; //bulk chemical potential
    static constexpr auto _xS = Dune::Indices::_4; //surface coordinates
    static constexpr auto _kappaVec = Dune::Indices::_5; //curvature vector on the surface
    static constexpr auto _lambda1 = Dune::Indices::_6; //principal stretch; artificial surfTen in incompressible case
    static constexpr auto _lambda2 = Dune::Indices::_7; //principal stretch
    static constexpr auto _fShell = Dune::Indices::_8; //shell force
    static constexpr auto _vS = Dune::Indices::_9; //restriction of velocity to surface
    static constexpr auto _phiS = Dune::Indices::_10; //restriction of phi to surface

    //modelExtensions
    static constexpr auto _phiGamma = Dune::Indices::_11; //additional phase field on the surface
    static constexpr auto _muGamma = Dune::Indices::_12; //chem. potential of additional phase field
}

#endif //BREASTCANCER_MYINDICES_H
