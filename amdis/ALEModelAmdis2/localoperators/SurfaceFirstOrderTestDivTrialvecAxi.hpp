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
        struct surface_test_divtrialvec_axi2 {
            std::size_t inside = 0; // integrate over the shell from inside (0) or outside (1)
            IS const &indexSet = {}; //the indexSet of the leafGridView of the grid!
            std::vector<int> const& partitions = {}; //partition numbers of each element
            std::vector<int> const& facets = {}; //facet number of the surface facet (-1 if element i not touching the surface)
        };
    }


    /// first-order operator \f$ \langle\psi, c\,\nabla\cdot\Phi\rangle \f$
    template<class IS = int>
    class SurfaceFirstOrderTestDivTrialvecAxi2 {
        IS const &indexSet_ = {};
        std::vector<int> const& partitions_;
        std::vector<int> const& facets_;
        std::size_t inside_;
    public:
        SurfaceFirstOrderTestDivTrialvecAxi2(tag::surface_test_divtrialvec_axi2<IS> t)
                : indexSet_(t.indexSet), partitions_(t.partitions), facets_(t.facets), inside_(t.inside) {}

        template<class CG, class RN, class CN, class Quad, class LocalFct, class Mat>
        void assemble(CG const &contextGeo, RN const &rowNode, CN const &colNode,
                      Quad const &quad, LocalFct const &localFct, Mat &elementMatrix) const {
            static_assert(static_size_v<typename LocalFct::Range> == 1, "Expression must be of scalar type.");
            static_assert(RN::isLeaf && CN::isPower,
                          "row-node must be Leaf-Node and col-node must be a Power-Node.");

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

            using RangeFieldType = typename CN::ChildType::LocalBasis::Traits::RangeFieldType;
            using WorldVector = Dune::FieldVector<RangeFieldType, CG::dow>;
            std::vector<WorldVector> colGradients, colGradientsAux;

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

                    // The transposed inverse Jacobian of the map from the reference facet to the real facet
                    const auto jacobian = geo.jacobianInverseTransposed(qp.position());
                    const auto jacobianBulk = geo.jacobianInverseTransposed(qp.position());

                    // The transposed inverse Jacobian of the map from the reference facet to the reference element
                    const auto jacobianRef = localGeo.jacobianTransposed(qp.position());

                    // The multiplicative factor in the integral transformation formula
                    const auto factor = localFct(local) * geo.integrationElement(qp.position()) * qp.weight();

                    // the values of the shape functions on the reference element at the quadrature point
                    auto const &shapeValues = rowNode.localBasisValuesAt(local);

                    // The gradients of the shape functions on the reference element
                    auto const &shapeGradients = colNode.child(0).localBasisJacobiansAt(local);
                    auto const &colShapeValues = colNode.child(1).localBasisValuesAt(local);

                    // Compute the shape function gradients on the real element
                    colGradients.resize(shapeGradients.size());
                    colGradientsAux.resize(shapeGradients.size());

                    // transform local gradients to reference intersection if necessary and compute
                    // gradients on the real intersection
                    auto dim = basisDim(colNode.child(0),
                                        Dune::PriorityTag<2>{}); // = dow-1 if surfaceLagrange node, = dow else
                    int dow = CG::dow;
                    if (dim == dow) {
                        for (std::size_t i = 0; i < colGradients.size(); ++i)
                            jacobianRef.mv(shapeGradients[i][0], colGradientsAux[i]);

                        for (std::size_t i = 0; i < colGradients.size(); ++i)
                            jacobian.mv(colGradientsAux[i], colGradients[i]);
                    } else {
                        for (std::size_t i = 0; i < colGradients.size(); ++i)
                            jacobian.mv(shapeGradients[i][0], colGradients[i]);
                    }


                    auto tangential = (geo.corner(0) - geo.corner(1)) / (geo.corner(0) - geo.corner(1)).two_norm(); // direction of normal doesn't matter here
                    auto normal = tangential;
                    normal[0] = -tangential[1];
                    normal[1] = tangential[0];


                    for (std::size_t i = 0; i < rowSize; ++i) {
                        const auto local_i = rowNode.localIndex(i);
                        const auto value = factor * shapeValues[i];
                        for (std::size_t j = 0; j < colSize; ++j) {
                            for (std::size_t k = 0; k < CG::dow; ++k) {
                                const auto local_kj = colNode.child(k).localIndex(j);
                                if (geo.corner(i)[1] > 1.e-8) {
                                    elementMatrix[local_i][local_kj] += value * colGradients[j][k];
                                } else { // take only the bulk derivative in y-direction (d_r v_r) by division of P[1][1] from the surface gradient
                                    if (k==1)  elementMatrix[local_i][local_kj] += value * colShapeValues[j] / geo.global(qp.position())[1];
                                }
                            }
                        }
                    }
                }
            }

        }
    };

    template<class LC, class IS>
    struct GridFunctionOperatorRegistry<tag::surface_test_divtrialvec_axi2<IS>, LC> {
        static constexpr int degree = 1;
        using type = SurfaceFirstOrderTestDivTrialvecAxi2<IS>;
    };

    /** @} **/

} // end namespace AMDiS
