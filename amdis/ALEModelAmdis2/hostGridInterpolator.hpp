#pragma once

#include <vector>

#include <amdis/common/FakeContainer.hpp>
#include <amdis/common/FieldMatVec.hpp>
#include <amdis/functions/EntitySet.hpp>
#include <amdis/functions/HierarchicNodeToRangeMap.hpp>
#include <amdis/functions/NodeIndices.hpp>
#include <amdis/linearalgebra/Traits.hpp>
#include <amdis/operations/Assigner.hpp>
#include <amdis/typetree/Traversal.hpp>

namespace AMDiS
{
    // interpolate a DiscreteGridViewFunction on a CurvedGrid to a DiscreteFunction on the host grid in order to
    // interpolate the values of the coordinates that are used for the grid movement
    template <class Coeff, class GridFct>
    void interpolateHostGrid(Coeff&& result, GridFct const& gf)
    {
        auto&& basisHost_ = result.basis();
        auto localViewHost = basisHost_.localView();

        result.coefficients().init(basisHost_, false);
        for (const auto& eH : entitySet(basisHost_))
        {
            localViewHost.bind(eH);

            auto const &nodeHost = localViewHost.tree();
            for (int i = 0; i < nodeHost.child(0).size(); i++) {
                for (int j = 0; j < nodeHost.degree(); j++) {
                    result.coefficients().set(localViewHost.index(nodeHost.child(j).localIndex(i)),
                                              gf.coefficients().get(localViewHost.index(nodeHost.child(j).localIndex(i))));
                }
            }
        }
        result.coefficients().finish();
    }
} // end namespace AMDiS
