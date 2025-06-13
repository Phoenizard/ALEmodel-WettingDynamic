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
        struct surface_coupling {
            unsigned int axi = 0; // axisymmetric or not?
            int inside = 0; // integrate over the shell from inside?
            int outside = 1; // integrate over the shell from outside?
            IS const &indexSet = {}; //the indexSet of the leafGridView of the grid!
            std::vector<int> const &partitions = {}; //partition numbers of each element
            std::vector<int> const &facets = {}; //facet number of the surface facet (-1 if element i not touching the surface)
        };        
    }


    /// Operator for coupling terms between a surface and a bulk solution (e.g. v in the bulk and v_S on the surface,
    /// where v=v_S on \Gamma should hold.
    template<class IS = int>
    class Surface_coupling {
        unsigned int axi_;
        IS const &indexSet_ = {};
        std::vector<int> const &partitions_;
        std::vector<int> const &facets_;
        int inside_;
        int outside_;
    public:
        Surface_coupling(tag::surface_coupling<IS> t)
                : axi_(t.axi), indexSet_(t.indexSet),
                  partitions_(t.partitions), facets_(t.facets), inside_(t.inside) , outside_(t.outside) {}

        // assemble method for LHS
        template<class CG, class Node, class Quad, class LocalFct, class Mat>
        void assemble(CG const &contextGeo, Node const &tree, Node const /*colNode*/,
                      Quad const &quad, LocalFct const &localFct, Mat &elementMatrix) const {
            static_assert(static_size_v<typename LocalFct::Range> == 1,
                          "Expression must be of scalar type.");

            using VelocityNode = Dune::TypeTree::Child<Node, _v.value>;
            using VelocitySNode = Dune::TypeTree::Child<Node, _vS.value>;
            using PhiNode = Dune::TypeTree::Child<Node, _phi.value>;
            using PhiSNode = Dune::TypeTree::Child<Node, _phiS.value>;


            static_assert(VelocityNode::isPower && VelocitySNode::isPower && PhiNode::isLeaf && PhiSNode::isLeaf,
                          "Nodes must have correct dimensions.");

            assert(!(inside_ && outside_)); // phiS on both sides not implemented, use classical coupling operators instead

            auto const &vNode = tree.child(_v);
            auto const &vSNode = tree.child(_vS);
            auto const &phiNode = tree.child(_phi);
            auto const &phiSNode = tree.child(_phiS);

            if (vNode.child(0).size() == 0 || vSNode.child(0).size() == 0
                || phiNode.size() == 0 || phiSNode.size() == 0)
                return;

            // compute facet number and partition number from surface basis node (if no surface basis node: return -2)
            auto part = std::max(partition(phiSNode, Dune::PriorityTag<2>{}),
                                 partition(vSNode.child(0), Dune::PriorityTag<2>{}));
            auto face = std::max(facet(phiSNode, Dune::PriorityTag<2>{}),
                                 facet(vSNode.child(0), Dune::PriorityTag<2>{}));

            if (face == -2) {
                auto index = idx(indexSet_, contextGeo.element(), Dune::PriorityTag<2>{});
                test_exit(!(partitions_.size() == 0 || facets_.size() == 0 || index == -1),
                          "please provide partitions, facets and indexSet of the leafGridView when test and trial "
                          "functions are not from the surfaceLagrange basis!");
                face = facets_[index];
                part = partitions_[index];
            }

            if (face == -1) // do nothing if element has no facet on the shell
                return;

            std::size_t vSize = vNode.child(0).size();
            std::size_t vSSize = vSNode.child(0).size();
            std::size_t phiSize = phiNode.size();
            std::size_t phiS_Size = phiSNode.size();

            auto element = contextGeo.element();

            auto intersection = element.template subEntity<1>(face);

            constexpr auto intersectionDim = std::decay_t<decltype(intersection)>::mydimension;
            const auto &quadShell = Dune::QuadratureRules<double, intersectionDim>::rule(intersection.type(), 2);

            // apply the operator only to intersections on the shell
            auto geo = intersection.geometry();

            if (part == 1 && inside_ || part == 0 && outside_) {
                for (auto const &qp: quadShell) {
                    auto localGeo = referenceElement(element).template geometry<1>(face);

                    // Position of the current quadrature point in the reference element
                    auto &&local = localGeo.global(qp.position());
                    auto &&global = contextGeo.elementGeometry().global(local);

                    // The multiplicative factor in the integral transformation formula
                    const auto factor = geo.integrationElement(qp.position()) * qp.weight();
                    const auto one = localFct(local); //not needed, since = 1.0

                    auto const &phiValues = phiNode.localBasisValuesAt(local);
                    auto const &phiSValues = phiSNode.localBasisValuesAt(local);
                    auto const &vValues = vNode.child(0).localBasisValuesAt(local);
                    auto const &vSValues = vSNode.child(0).localBasisValuesAt(local);

                    // vS = v on Gamma
                    for (std::size_t i = 0; i < vSSize; ++i) {
                        const auto value = factor * vSValues[i];

                        // vS = ...
                        for (std::size_t k = 0; k < vSNode.degree(); ++k) {
                            const auto local_ki = vSNode.child(k).localIndex(i);
                            elementMatrix[local_ki][local_ki] += value * vSValues[i];
                        }
                        for (std::size_t j = i + 1; j < vSSize; ++j) {
                            for (std::size_t k = 0; k < vSNode.degree(); ++k) {
                                const auto local_ki = vSNode.child(k).localIndex(i);
                                const auto local_kj = vSNode.child(k).localIndex(j);

                                elementMatrix[local_ki][local_kj] += value * vSValues[j];
                                elementMatrix[local_kj][local_ki] += value * vSValues[j];
                            }
                        }

                        // vS = v on Gamma
                        for (std::size_t j = 0; j < vSize; ++j) {
                            for (std::size_t k = 0; k < CG::dow; ++k) {
                                const auto local_kj = vNode.child(k).localIndex(j);
                                const auto local_ki = vSNode.child(k).localIndex(i);
                                elementMatrix[local_ki][local_kj] += - value * vValues[j];
                            }
                        }
                    }

                    // phiS = phi on Gamma (only outside, if phi exists on the inside, add another phiS!
                    for (std::size_t i = 0; i < phiS_Size; ++i) {
                        const auto value = factor * phiSValues[i];

                        // phiS = ...
                        const auto local_i = phiSNode.localIndex(i);
                        elementMatrix[local_i][local_i] += value * phiSValues[i];

                        for (std::size_t j = i + 1; j < phiS_Size; ++j) {
                            const auto local_j = phiSNode.localIndex(j);

                            elementMatrix[local_i][local_j] += value * phiSValues[j];
                            elementMatrix[local_j][local_i] += value * phiSValues[j];
                        }

                        // phiS = phi on Gamma
                        for (std::size_t j = 0; j < phiSize; ++j) {
                            const auto local_j = phiNode.localIndex(j);
                            elementMatrix[local_i][local_j] += - value * phiValues[j];
                        }
                    }
                }
            }
        }
    };


    template<class LC, class IS>
    struct GridFunctionOperatorRegistry<tag::surface_coupling<IS>, LC> {
        static constexpr int degree = 0;
        using type = Surface_coupling<IS>;
    };
    /** @} **/

} // end namespace AMDiS
