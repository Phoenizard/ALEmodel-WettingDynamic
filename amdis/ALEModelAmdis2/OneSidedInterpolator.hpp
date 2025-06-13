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
        template <class Coeff, class GridFct, class TP>
        void interpolateOneSided(Coeff& result, GridFct const& gf, int side, std::vector<int> partitions, TP treePath)
        {
            // Obtain a local view of the gridFunction
            auto&& basis_ = result.basis();
            auto lf = localFunction(gf);
            auto localView = basis_.localView();
            auto gv = basis_.gridView();
            auto const& indexSet0 = basis_.gridView().grid().levelGridView(0).indexSet();

            std::vector<typename Coeff::value_type> localCoeff;

            result.init(basis_, false);
            for (const auto& e : entitySet(basis_))
            {
                localView.bind(e);
                lf.bind(e);

                auto inside = e;
                while (inside.hasFather()) {
                    inside = inside.father();
                }
                int p = partitions[indexSet0.index(inside)];

                if (p == side) {
                    //auto &&subTree = Dune::TypeTree::child(localView.tree(), treePath);
                    //Traversal::forEachLeafNode(subTree, [&](auto const &node, auto const &tp) {
                        // compute local interpolation coefficients
                        localView.tree().finiteElement().localInterpolation().interpolate(lf,localCoeff);

                        // assign local coefficients to global vector
                        result.scatter(localView, localView.tree(), localCoeff, Assigner::assign{});
                  //  });
                }
                    lf.unbind();

            }
            result.finish();
        }


} // end namespace AMDiS
