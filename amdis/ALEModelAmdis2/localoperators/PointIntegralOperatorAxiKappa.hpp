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
        template <class IS = int>
        struct surface_point_test_trial {
            std::size_t coord = 0; //coordinate direction, 0: x, 1: y
            std::size_t inside = 0; // integrate over the shell from inside (0) or outside (1)
            IS const &indexSet = {}; //the indexSet of the leafGridView of the grid!
            std::vector<int> const&  partitions = {}; //partition numbers of each element
            std::vector<int> const&  facets = {}; //facet number of the surface facet (-1 if element i not touching the surface)
        };
    }


    /// zero-order vector-operator \f$ (\mathbf{b}\cdot\Psi) \f$
    template <class IS = int>
    class SurfacePointIntegralOperatorAxiKappa
    {
        IS const &indexSet_ = {};
        std::vector<int> const& partitions_;
        std::vector<int> const& facets_;
        std::size_t inside_;
        std::size_t coord_; //coordinate direction, 0: x, 1: y
    public:
        SurfacePointIntegralOperatorAxiKappa(tag::surface_point_test_trial<IS> t)
                : indexSet_(t.indexSet), partitions_(t.partitions), facets_(t.facets), inside_(t.inside), coord_(t.coord)
        {}

        template <class CG, class RN, class CN, class Quad, class LocalFct, class Mat>
        void assemble(CG const& contextGeo, RN const& rowNode, CN const& colNode,
                      Quad const& quad, LocalFct const& localFct, Mat& elementMatrix) const
        {
            static_assert(static_size_v<typename LocalFct::Range> == 1,
                          "Expression must be of scalar type." );
            static_assert(RN::isLeaf,
                          "Operator can be applied to Leaf-Nodes only.");

            if (rowNode.size() == 0 || colNode.size() == 0 )
                return;

            // compute facet number and partition number from surface basis node (if no surface basis node: return -2)
            auto part = partition(rowNode, Dune::PriorityTag<2>{});
            auto face = facet(rowNode, Dune::PriorityTag<2>{});

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

            std::size_t size = rowNode.size();

            auto element = contextGeo.element();

            auto intersection = element.template subEntity<1>(face);

            constexpr auto intersectionDim = std::decay_t<decltype(intersection)>::mydimension;
            const auto& quadShell = Dune::QuadratureRules<double,intersectionDim>::rule(intersection.type(), 2);

            // apply the operator only to intersections on the shell
            if (part == 0) {
                auto geo = intersection.geometry();

                //do on facet, which touches the symmetry axis
                if (geo.corner(0)[1] < 1e-10 || geo.corner(1)[1] < 1e-10) {
                    Dune::FieldVector<double,CG::dow> conormal;
                    conormal = geo.corner(0) - geo.corner(1);
                    int i = 0;
                    auto s = geo.corner(0);
                    auto y1 = geo.corner(1)[1];
                    if (geo.corner(1)[1] < 1e-10) {
                        s = geo.corner(1);
                        i = 1;
                        conormal = -conormal;
                        y1 = geo.corner(0)[1]*0.5;
                    }
                    conormal = conormal/conormal.two_norm();
                    auto ny1 = -conormal[0];

                    auto localGeo = referenceElement(element).template geometry<1>(face);
                    // Position of the current quadrature point in the reference element
                    auto&& local = localGeo.global(geo.local(s));

                    // The multiplicative factor in the integral transformation formula
                    const auto exprValue = localFct(local);
                    auto const &shapeValues = rowNode.localBasisValuesAt(local);

                    const auto value = exprValue;
                    const auto local_ki = colNode.localIndex(i);
                    const auto local_kj = rowNode.localIndex(i);
                    if (inside_) {
                        elementMatrix[local_kj][local_ki] += value * ny1/y1 * conormal[coord_];
                    } else {
                        elementMatrix[local_kj][local_ki] += value * conormal[coord_];
                    }
                }
            }
        }

    template <class CG, class RN,  class Quad, class LocalFct, class Mat>
    void assemble(CG const& contextGeo, RN const& rowNode,
                  Quad const& quad, LocalFct const& localFct, Mat& elementMatrix) const
    {
        static_assert(static_size_v<typename LocalFct::Range> == 1,
                      "Expression must be of scalar type." );
        static_assert(RN::isLeaf,
                      "Operator can be applied to Leaf-Nodes only.");

        if (rowNode.size() == 0 )
            return;

        // compute facet number and partition number from surface basis node (if no surface basis node: return -2)
        auto part = partition(rowNode, Dune::PriorityTag<2>{});
        auto face = facet(rowNode, Dune::PriorityTag<2>{});

        if (face == -2) {
            auto index = idx(indexSet_, contextGeo.element(), Dune::PriorityTag<2>{});
            test_exit(!(partitions_.size() == 0 || facets_.size() == 0 || index == -1),
                      "please provide partitions, facets and indexSet of the leafGridView when test and trial "
                      "functions are not from the surfaceLagrange basis!");
            face = facets_[index];
            part = partitions_[index];
        }

        if (face == -1 || part != 0) // do nothing if the element is on the wrong partition
            return;

        std::size_t size = rowNode.size();

        auto element = contextGeo.element();

        auto intersection = element.template subEntity<1>(face);

        constexpr auto intersectionDim = std::decay_t<decltype(intersection)>::mydimension;
        const auto& quadShell = Dune::QuadratureRules<double,intersectionDim>::rule(intersection.type(), 2);

        // apply the operator only to intersections on the shell
        if (part == 0) {
            auto geo = intersection.geometry();

            //do on facet, which touches the symmetry axis
            if (geo.corner(0)[1] < 1e-10 || geo.corner(1)[1] < 1e-10) {
                Dune::FieldVector<double,CG::dow> conormal;
                conormal = geo.corner(0) - geo.corner(1);
                int i = 0;
                auto s = geo.corner(0);
                auto y1 = geo.corner(1)[1];
                auto x1 = geo.corner(1)[0];
                if (geo.corner(1)[1] < 1e-10) {
                    s = geo.corner(1);
                    i = 1;
                    conormal = -conormal;
                    y1 = geo.corner(0)[1];
                    x1 = geo.corner(0)[0];
                }
                conormal = conormal/conormal.two_norm();
                auto ny1 = (s[0] < 0) ? -conormal[0] : conormal[0];
                auto sign = (x1 < 0 && x1 < s[0] || x1 > 0 && x1 > s[0]) ? -1 : 1;

                auto localGeo = referenceElement(element).template geometry<1>(face);
                // Position of the current quadrature point in the reference element
                auto&& local = localGeo.global(geo.local(s));

                // The multiplicative factor in the integral transformation formula
                const auto exprValue = localFct(local);
                auto const &shapeValues = rowNode.localBasisValuesAt(local);

                const auto value = exprValue;
                const auto local_kj = rowNode.localIndex(i);
                if (inside_) {
                    elementMatrix[local_kj] += sign * value * 2.0 * ny1/y1 * conormal[coord_];
                } else {
                    elementMatrix[local_kj] += value * conormal[coord_];
                }
            }
        }
    }
};

    template <class LC, class IS>
    struct GridFunctionOperatorRegistry<tag::surface_point_test_trial<IS>, LC>
    {
        static constexpr int degree = 0;
        using type = SurfacePointIntegralOperatorAxiKappa<IS>;
    };

    /** @} **/

} // end namespace AMDiS
