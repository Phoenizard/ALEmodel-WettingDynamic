#pragma once

#include <type_traits>

#include <amdis/GridFunctionOperator.hpp>
#include <amdis/common/StaticSize.hpp>

namespace AMDiS
{
    /**
     * \addtogroup operators
     * @{
     **/

    namespace tag
    {
        struct point_testvec {};

        template <class IS = int>
        struct surface_point_testvec {
            std::size_t inside = 0; // integrate over the shell from inside (0) or outside (1)
            IS const &indexSet = {}; //the indexSet of the leafGridView of the grid!
            std::vector<int> const&  partitions = {}; //partition numbers of each element
            std::vector<int> const&  facets = {}; //facet number of the surface facet (-1 if element i not touching the surface)
        };
    }


    /// zero-order vector-operator \f$ (\mathbf{b}\cdot\Psi) \f$
    template <class IS = int>
    class SurfacePointIntegralOperatorAxi
    {
        IS const &indexSet_ = {};
        std::vector<int> const& partitions_;
        std::vector<int> const& facets_;
        std::size_t inside_;
    public:
        SurfacePointIntegralOperatorAxi(tag::surface_point_testvec<IS> t)
                : indexSet_(t.indexSet), partitions_(t.partitions), facets_(t.facets), inside_(t.inside)
        {}

        template <class CG, class Node, class Quad, class LocalFct, class Vec>
        void assemble(CG const& contextGeo, Node const& node, Quad const& quad,
                      LocalFct const& localFct, Vec& elementVector) const
        {
            static_assert(static_size_v<typename LocalFct::Range> == 1,
                          "Expression must be of scalar type." );
            static_assert(Node::isPower,
                          "Operator can be applied to Power-Nodes only.");

            assert(node.degree() == CG::dow);

            if (node.child(0).size()==0)
                return;

            // compute facet number and partition number from surface basis node (if no surface basis node: return -2)
            auto part = partition(node.child(0), Dune::PriorityTag<2>{});
            auto face = facet(node.child(0), Dune::PriorityTag<2>{});

            if (face == -2) {
                auto index = idx(indexSet_, contextGeo.element(), Dune::PriorityTag<2>{});
                test_exit(!(partitions_.size() == 0 || facets_.size() == 0 || index == -1),
                          "please provide partitions, facets and indexSet of the leafGridView when test and trial "
                          "functions are not from the surfaceLagrange basis!");
                face = facets_[index];
                part = partitions_[index];
            }

            if (face == -1 || part != inside_) // do nothing if the element is on the wrong partition
                return;

            std::size_t size = node.child(0).size();

            auto element = contextGeo.element();

            auto intersection = element.template subEntity<1>(face);

            constexpr auto intersectionDim = std::decay_t<decltype(intersection)>::mydimension;
            const auto& quadShell = Dune::QuadratureRules<double,intersectionDim>::rule(intersection.type(), 2);

            // apply the operator only to intersections on the shell
            if (part == inside_) {
                auto geo = intersection.geometry();

                //do on facet, which touches the symmetry axis
                if (geo.corner(0)[1] < 1e-10 || geo.corner(1)[1] < 1e-10) {
                    Dune::FieldVector<double,CG::dow> conormal;
                    conormal = geo.corner(0) - geo.corner(1);
                    int i = 0;
                    auto s = geo.corner(0);
                    if (geo.corner(1)[1] < 1e-10) {
                        s = geo.corner(1);
                        i = 1;
                        conormal = -conormal;
                    }
                    conormal = conormal/conormal.two_norm();

                    auto pos = s;
                    // Position of the current quadrature point in the reference element
                    auto &&local = contextGeo.coordinateInElement(pos);
                    // The multiplicative factor in the integral transformation formula
                    const auto exprValue = localFct(local);

                    auto const &shapeValues = node.child(0).localBasisValuesAt(local);

                    const auto value = exprValue;
                    for (std::size_t k = 0; k < 2; ++k) {
                        const auto local_ki = node.child(k).localIndex(i);
                        elementVector[local_ki] += value * conormal[k];
                    }
                }
            }
        }
    };

    template <class LC, class IS>
    struct GridFunctionOperatorRegistry<tag::surface_point_testvec<IS>, LC>
    {
        static constexpr int degree = 0;
        using type = SurfacePointIntegralOperatorAxi<IS>;
    };



    /// zero-order vector-operator \f$ (\mathbf{b}\cdot\Psi) \f$
    class PointIntegralOperatorAxi
    {
    public:
        PointIntegralOperatorAxi(tag::point_testvec)
        {}

        template <class CG, class Node, class Quad, class LocalFct, class Vec>
        void assemble(CG const& contextGeo, Node const& node, Quad const& quad,
                      LocalFct const& localFct, Vec& elementVector) const {
            static_assert(static_size_v<typename LocalFct::Range> == 1,
                          "Expression must be of scalar type.");
            static_assert(Node::isPower,
                          "Operator can be applied to Power-Nodes only.");

            assert(node.degree() == CG::dow);


            auto element = contextGeo.element();
            auto geo = element.geometry();

            //do on facet, which touches the symmetry axis
            if (geo.corner(0)[1] < 1e-10 || geo.corner(1)[1] < 1e-10) {
                Dune::FieldVector<double,CG::dow> conormal;
                conormal = geo.corner(0) - geo.corner(1);
                int i = 0;
                auto s = geo.corner(0);
                if (geo.corner(1)[1] < 1e-10) {
                    s = geo.corner(1);
                    i = 1;
                    conormal = -conormal;
                }
                conormal = conormal/conormal.two_norm();
                auto pos = s;
                // Position of the current quadrature point in the reference element
                auto &&local = geo.local(pos);
                // The multiplicative factor in the integral transformation formula
                const auto exprValue = localFct(local);

                auto const &shapeValues = node.child(0).localBasisValuesAt(local);

                const auto value = exprValue;
                for (std::size_t k = 0; k < 2; ++k) {
                    const auto local_ki = node.child(k).localIndex(i);
                    elementVector[local_ki] += value * conormal[k];
                }
            }
        }
    };

    template <class LC>
    struct GridFunctionOperatorRegistry<tag::point_testvec, LC>
    {
        static constexpr int degree = 0;
        using type = PointIntegralOperatorAxi;
    };

    /** @} **/

} // end namespace AMDiS
