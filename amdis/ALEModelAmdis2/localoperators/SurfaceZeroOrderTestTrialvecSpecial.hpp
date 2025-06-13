#pragma once

#include <type_traits>

#include <amdis/GridFunctionOperator.hpp>
#include <amdis/common/StaticSize.hpp>

namespace AMDiS {
    /**
     * \addtogroup operators
     * @{
     **/

    namespace tag {
        template<class IS = int>
        struct surface_test_trialvec_special {
            int comp;
            std::size_t inside = 0; // integrate over the shell from inside (0) or outside (1)
            IS const &indexSet = {}; //the indexSet of the leafGridView of the grid!
            std::vector<int> const &partitions = {}; //partition numbers of each element
            std::vector<int> const &facets = {}; //facet number of the surface facet (-1 if element i not touching the surface)
        };
    }


    /// zero-order operator \f$ \langle\psi, \mathbf{b}\cdot\Phi\rangle \f$
    template<class IS = int>
    class SurfaceZeroOrderTestTrialvecSpecial {
        IS const &indexSet_ = {};
        std::vector<int> const &partitions_;
        std::vector<int> const &facets_;
        std::size_t inside_;
        int comp_;
    public:
        SurfaceZeroOrderTestTrialvecSpecial(tag::surface_test_trialvec_special<IS> t)
                : indexSet_(t.indexSet), partitions_(t.partitions), facets_(t.facets), inside_(t.inside), comp_(t.comp) {}

        template<class CG, class RN, class CN, class Quad, class LocalFct, class Mat>
        void assemble(CG const &contextGeo, RN const &rowNode, CN const &colNode,
                      Quad const &quad, LocalFct const &localFct, Mat &elementMatrix) const {
            static_assert(static_size_v<typename LocalFct::Range> == 1,
                          "Expression must be of scalar type.");
            static_assert(RN::isLeaf && CN::isPower,
                          "RN must be Leaf-Node and CN must be a Power-Node.");

            assert(colNode.degree() == CG::dow);

            if (rowNode.size() == 0 || colNode.child(0).size() == 0)
                return;

            // compute facet number and partition number from surface basis node (if no surface basis node: return -2)
            auto part = std::max(partition(rowNode, Dune::PriorityTag<2>{}),
                                 partition(colNode.child(0), Dune::PriorityTag<2>{}));
            auto face = std::max(facet(rowNode, Dune::PriorityTag<2>{}),
                                 facet(colNode.child(0), Dune::PriorityTag<2>{}));

            if (face == -2) {
                auto index = idx(indexSet_, contextGeo.element(), Dune::PriorityTag<2>{});
                test_exit(!(partitions_.size() == 0 || facets_.size() == 0 || index == -1),
                          "please provide partitions, facets and indexSet of the leafGridView when test and trial "
                          "functions are not from the surfaceLagrange basis!");
                face = facets_[index];
                part = partitions_[index];
            }

            if (face == -1 || part != inside_) // do nothing if element has no facet on the shell
                return;

            std::size_t rowSize = rowNode.size();
            std::size_t colSize = colNode.child(0).size();

            auto element = contextGeo.element();

            auto intersection = element.template subEntity<1>(face);

            constexpr auto intersectionDim = std::decay_t<decltype(intersection)>::mydimension;
            const auto& quadShell = Dune::QuadratureRules<double,intersectionDim>::rule(intersection.type(), 2);

            // apply the operator only to intersections on the shell
            if (part == inside_) {
                auto geo = intersection.geometry();

                for (auto const& qp : quadShell) {
                    auto localGeo = referenceElement(element).template geometry<1>(face);

                    // Position of the current quadrature point in the reference element
                    auto&& local = localGeo.global(qp.position());

                    // The multiplicative factor in the integral transformation formula
                    const auto factor = geo.integrationElement(qp.position()) * qp.weight();
                    const auto b = localFct(local);

                    auto const &rowShapeValues = rowNode.localBasisValuesAt(local);
                    auto const &colShapeValues = colNode.child(0).localBasisValuesAt(local);

                    //todo: make 3D ready
                    auto x0 = geo.corner(0);
                    auto x1 = geo.corner(1);
                    auto t = x0 - x1;
                    t = t / t.two_norm();
                    auto n = t;
                    n[0] = -t[1];
                    n[1] = t[0];

                    for (std::size_t i = 0; i < rowSize; ++i) {
                        const auto local_i = rowNode.localIndex(i);

                        for (std::size_t j = 0; j < colSize; ++j) {
                            const auto value = b * (factor * rowShapeValues[i] * colShapeValues[j]);

                            for (std::size_t k = 0; k < colNode.degree(); ++k) {
                                const auto local_kj = colNode.child(k).localIndex(j);
                                elementMatrix[local_i][local_kj] += value * n[k] * n[comp_];
                            }
                        }
                    }
                }
            }
        }
    };

    template<class LC, class IS>
    struct GridFunctionOperatorRegistry<tag::surface_test_trialvec_special<IS>, LC> {
        static constexpr int degree = 0;
        using type = SurfaceZeroOrderTestTrialvecSpecial<IS>;
    };

    /** @} **/

} // end namespace AMDiS
