#pragma once

#include <type_traits>

#include <amdis/GridFunctionOperator.hpp>
#include <amdis/common/StaticSize.hpp>
#include <amdis/typetree/FiniteElementType.hpp>

namespace AMDiS {
    /**
     * \addtogroup operators
     * @{
     **/

    namespace tag {
        template<class IS = int>
        struct surface_test_trial_axi2 {
            std::reference_wrapper<std::pair<Dune::FieldVector<double, 2>, double>> const& valLeft;
            std::reference_wrapper<std::pair<Dune::FieldVector<double, 2>, double>> const& valRight;
            std::size_t inside = 0; // integrate over the shell from inside (0) or outside (1)
            IS const &indexSet = {}; //the indexSet of the leafGridView of the grid!
            std::vector<int> const &partitions = {}; //partition numbers of each element
            std::vector<int> const &facets = {}; //facet number of the surface facet (-1 if element i not touching the surface)
        };
    }


    /// zero-order operator \f$ \langle\psi, c\,\phi\rangle \f$
    template<class IS = int>
    class SurfaceZeroOrderTestTrialAxi2 {
        IS const &indexSet_ = {};
        std::vector<int> const &partitions_;
        std::vector<int> const &facets_;
        std::size_t inside_;
        std::reference_wrapper<std::pair<Dune::FieldVector<double, 2>, double>> const& valLeft_;
        std::reference_wrapper<std::pair<Dune::FieldVector<double, 2>, double>> const& valRight_;

    public:
        SurfaceZeroOrderTestTrialAxi2(tag::surface_test_trial_axi2<IS> t)
                : indexSet_(t.indexSet), partitions_(t.partitions), facets_(t.facets), inside_(t.inside),
                valLeft_(t.valLeft), valRight_(t.valRight) {}

        template<class CG, class RN, class CN, class Quad, class LocalFct, class Mat>
        void assemble(CG const &contextGeo, RN const &rowNode, CN const &colNode,
                      Quad const &quad, LocalFct const &localFct, Mat &elementMatrix) const {
            static_assert(static_size_v<typename LocalFct::Range> == 1,
                          "Expression must be of scalar type.");

            if (rowNode.size() == 0 || colNode.size() == 0)
                return;

            // compute facet number and partition number from surface basis node (if no surface basis node: return -2)
            auto part = std::max(partition(rowNode, Dune::PriorityTag<2>{}),
                                 partition(colNode, Dune::PriorityTag<2>{}));
            auto face = std::max(facet(rowNode, Dune::PriorityTag<2>{}), facet(colNode, Dune::PriorityTag<2>{}));

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
            const bool sameFE = std::is_same_v<FiniteElementType_t<RN>, FiniteElementType_t<CN>>;
            const bool sameNode = rowNode.treeIndex() == colNode.treeIndex();

            if (sameFE && sameNode)
                getElementMatrixOptimized(contextGeo, quad, rowNode, colNode, localFct, elementMatrix, face, part);
            else
                getElementMatrixStandard(contextGeo, quad, rowNode, colNode, localFct, elementMatrix, face, part);
        }

        template<class CG, class Node, class Quad, class LocalFct, class Vec>
        void assemble(CG const &contextGeo, Node const &node,
                      Quad const &quad, LocalFct const &localFct, Vec &elementVector, int face, int part) const {
            static_assert(static_size_v<typename LocalFct::Range> == 1,
                          "Expression must be of scalar type.");
            static_assert(Node::isLeaf,
                          "Operator can be applied to Leaf-Nodes only");

            if (node.size() == 0)
                return;

            std::size_t size = node.size();

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
                    const auto factor = localFct(local) * geo.integrationElement(qp.position()) * qp.weight();

                    auto const &shapeValues = node.localBasisValuesAt(local);
                    for (std::size_t i = 0; i < size; ++i) {
                        const auto local_i = node.localIndex(i);
                        if (geo.corner(i)[1] > 1.e-8)
                            elementVector[local_i] += factor * shapeValues[i];
                    }
                }
            }
        }


    protected:
        template<class CG, class QR, class RN, class CN, class LocalFct, class Mat>
        void getElementMatrixStandard(CG const &contextGeo, QR const &quad,
                                      RN const &rowNode, CN const &colNode,
                                      LocalFct const &localFct, Mat &elementMatrix, int face, int part) const {
            static_assert(RN::isLeaf && CN::isLeaf,
                          "Operator can be applied to Leaf-Nodes only.");

            std::size_t rowSize = rowNode.size();
            std::size_t colSize = colNode.size();

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
                    const auto factor = localFct(local) * geo.integrationElement(qp.position()) * qp.weight();

                    auto const &rowShapeValues = rowNode.localBasisValuesAt(local);
                    auto const &colShapeValues = colNode.localBasisValuesAt(local);

                    for (std::size_t i = 0; i < rowSize; ++i) {
                        const auto local_i = rowNode.localIndex(i);
                        const auto value = factor * rowShapeValues[i];

                        for (std::size_t j = 0; j < colSize; ++j) {
                            const auto local_j = colNode.localIndex(j);
                            if (geo.corner(i)[1] > 1.e-8) {
                                elementMatrix[local_i][local_j] += value * colShapeValues[j];
                            } else {
                                if (std::abs(geo.corner(i)[0] - valLeft_.get().first[0]) < std::abs(geo.corner(i)[0] - valRight_.get().first[0])) {
                                    auto value2 = valLeft_.get().second / geo.global(qp.position())[1];
                                    elementMatrix[local_i][local_j] +=  value2 * colShapeValues[j];
                                }
                                if (std::abs(geo.corner(i)[0] - valLeft_.get().first[0]) > std::abs(geo.corner(i)[0] - valRight_.get().first[0])) {
                                    auto value2 = valRight_.get().second / geo.global(qp.position())[1];
                                    elementMatrix[local_i][local_j] +=  value2 * colShapeValues[j];
                                }
                            }
                        }
                    }
                }
            }
        }


        template<class CG, class QR, class RN, class CN, class LocalFct, class Mat>
        void getElementMatrixOptimized(CG const &contextGeo, QR const &quad,
                                       RN const &node, CN const & /*colNode*/,
                                       LocalFct const &localFct, Mat &elementMatrix, int face, int part) const {
            static_assert(RN::isLeaf && CN::isLeaf,
                          "Operator can be applied to Leaf-Nodes only.");

            std::size_t size = node.size();

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
                    const auto factor = localFct(local) * geo.integrationElement(qp.position()) * qp.weight();

                    auto const &shapeValues = node.localBasisValuesAt(local);

                    for (std::size_t i = 0; i < size; ++i) {
                        const auto local_i = node.localIndex(i);

                        const auto value = factor * shapeValues[i];
                        if (geo.corner(i)[1] > 1.e-8)
                            elementMatrix[local_i][local_i] += value * shapeValues[i];

                        for (std::size_t j = i + 1; j < size; ++j) {
                            const auto local_j = node.localIndex(j);

                            if (geo.corner(i)[1] > 1.e-8) {
                                elementMatrix[local_i][local_j] += value * shapeValues[j];
                                elementMatrix[local_j][local_i] += value * shapeValues[j];
                            }
                        }
                    }
                }
            }
        }
    };

    template<class LC, class IS>
    struct GridFunctionOperatorRegistry<tag::surface_test_trial_axi2<IS>, LC> {
        static constexpr int degree = 0;
        using type = SurfaceZeroOrderTestTrialAxi2<IS>;
    };


    /// Create a zero-order term
    template<class Expr, class IS>
    auto zotAxi2(Expr &&expr, std::vector<int> partitions, IS is, int quadOrder = -1) {
        return makeOperator(tag::surface_test_trial_axi2<IS>{partitions, is}, FWD(expr), quadOrder);
    }

    /** @} **/

} // end namespace AMDiS
