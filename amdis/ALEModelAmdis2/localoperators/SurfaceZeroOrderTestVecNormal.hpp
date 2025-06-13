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
        struct surface_testvec_normal {
            std::size_t inside = 0; // integrate over the shell from inside (0) or outside (1)
            IS const &indexSet = {}; //the indexSet of the leafGridView of the grid!
            std::vector<int> const &partitions = {}; //partition numbers of each element
            std::vector<int> const &facets = {}; //facet number of the surface facet (-1 if element i not touching the surface)
        };
    }

    /// zero-order vector-operator \f$ (c\, \psi) \f$
    template<class IS = int>
    class SurfaceZeroOrderTestVecNormal {
        IS const& indexSet_;
        std::vector<int> const &partitions_;
        std::vector<int> const &facets_;
        std::size_t inside_;
    public:
        SurfaceZeroOrderTestVecNormal(tag::surface_testvec_normal<IS> t)
                : indexSet_(t.indexSet), partitions_(t.partitions), facets_(t.facets), inside_(t.inside) {}

        template<class CG, class Node, class Quad, class LocalFct, class Vec>
        void assemble(CG const &contextGeo, Node const &node, Quad const &quad,
                      LocalFct const &localFct, Vec &elementVector) const {
            static_assert(static_size_v<typename LocalFct::Range> == 1,
                          "Expression must be of scalar type.");
            static_assert(Node::isPower,
                          "Operator can be applied to Power-Nodes only");

            assert(node.degree() == CG::dow);
            
            if (node.child(0).size() == 0)
                return;

            // compute facet number and partition number from surface basis node (if no surface basis node: return -2)
            auto part = partition(node.child(0), Dune::PriorityTag<2>{});
            auto face = facet(node.child(0), Dune::PriorityTag<2>{});

            // if no surface basis node, determine part and face from partitions_ and facets_ vectors
            if (face == -2) {
                auto index = idx(indexSet_, contextGeo.element(), Dune::PriorityTag<2>{});
                test_exit(!(partitions_.size() == 0 || facets_.size() == 0 || index == -1),
                          "please provide partitions, facets and indexSet of the leafGridView when test and trial "
                          "functions are not from the surfaceLagrange basis!");
                face = facets_[index];
                part = partitions_[index];
            }

            if (face == -1 || part != inside_) // do nothing if element has no facet on shell
                return;

            std::size_t size = node.child(0).size();

            auto element = contextGeo.element();

            auto intersection = element.template subEntity<1>(face);

            constexpr auto intersectionDim = std::decay_t<decltype(intersection)>::mydimension;
            const auto& quadShell = Dune::QuadratureRules<double,intersectionDim>::rule(intersection.type(), 2);

            // apply the operator only to intersections on the shell
            if (part == inside_) {
                auto geo = intersection.geometry();
                auto x1 = geo.corner(0);
                auto x2 = geo.corner(1);
                auto t = x2 - x1;
                auto n = t;
                n[0] = -t[1];
                n[1] = t[0];
                for (int i = 0; i < element.geometry().corners(); i++) {
                    if ((element.geometry().corner(i) - x1).two_norm() < 1.e-9
                        && (element.geometry().corner(i) - x2).two_norm() < 1.e-9) {
                        auto aux = geo.center() - element.geometry().corner(i);
                        if (n.dot(aux) < 0) n=-n;
                    }
                }
                auto aux = element.geometry().corner(0);

                for (auto const& qp : quadShell) {
                    auto localGeo = referenceElement(element).template geometry<1>(face);

                    // Position of the current quadrature point in the reference element
                    auto &&local = localGeo.global(qp.position());

                    // The multiplicative factor in the integral transformation formula
                    const auto factor = geo.integrationElement(qp.position()) * qp.weight();
                    const auto exprValue = localFct(local);

                    auto const &shapeValues = node.child(0).localBasisValuesAt(local);
                    for (std::size_t i = 0; i < size; ++i) {
                        const auto value = exprValue * (factor * shapeValues[i]);
                        for (std::size_t k = 0; k < CG::dow; ++k) {
                            const auto local_ki = node.child(k).localIndex(i);
                            elementVector[local_ki] += value * Dune::at(n,k);
                        }
                    }
                }
            }
        }
    };

    template<class LC, class IS>
    struct GridFunctionOperatorRegistry<tag::surface_testvec_normal<IS>, LC> {
        static constexpr int degree = 0;
        using type = SurfaceZeroOrderTestVecNormal<IS>;
    };

    /** @} **/

} // end namespace AMDiS
