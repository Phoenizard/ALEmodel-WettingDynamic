#pragma once

#include <type_traits>

#include <amdis/GridFunctionOperator.hpp>
#include <amdis/Output.hpp>
#include <amdis/common/StaticSize.hpp>
#include <amdis/common/ValueCategory.hpp>
#include <amdis/typetree/FiniteElementType.hpp>

namespace AMDiS
{
    /**
     * \addtogroup operators
     * @{
     **/

    namespace tag
    {
        struct gradtest_gradtrial_axi {};
    }


    /// second-order operator \f$ \langle\nabla\psi, c\,\nabla\phi\rangle \f$, or \f$
    class SecondOrderGradTestGradTrialAxi
    {
    public:
        SecondOrderGradTestGradTrialAxi(tag::gradtest_gradtrial_axi) {}

        template <class CG, class RN, class CN, class Quad, class LocalFct, class Mat>
        void assemble(CG const& contextGeo, RN const& rowNode, CN const& colNode,
                      Quad const& quad, LocalFct const& localFct, Mat& elementMatrix) const
        {
            using expr_value_type = typename LocalFct::Range;
            static_assert(static_size_v<expr_value_type> == 1 ||
                          (static_num_rows_v<expr_value_type> == CG::dow &&
                           static_num_cols_v<expr_value_type> == CG::dow),
                          "Expression must be of scalar or matrix type.");
            static_assert(RN::isLeaf && CN::isLeaf,
                          "Operator can be applied to Leaf-Nodes only.");

            const bool sameFE = std::is_same_v<FiniteElementType_t<RN>, FiniteElementType_t<CN>>;
            const bool sameNode = rowNode.treeIndex() == colNode.treeIndex();
            using Category = ValueCategory_t<typename LocalFct::Range>;

            if (sameFE && sameNode)
                getElementMatrixOptimized(contextGeo, quad, rowNode, colNode, localFct,  elementMatrix, Category{});
            else if (sameFE)
                getElementMatrixStandard(contextGeo, quad, rowNode, colNode, localFct, elementMatrix);
            else
                error_exit("Not implemented: currently only the implementation for equal fespaces available");
        }

    protected:

        template <class CG, class QR, class RN, class CN, class LocalFct, class Mat>
        void getElementMatrixStandard(CG const& contextGeo, QR const& quad,
                                      RN const& rowNode, CN const& colNode,
                                      LocalFct const& localFct, Mat& elementMatrix) const
        {
            std::size_t size = rowNode.size();

            using RangeFieldType = typename RN::LocalBasis::Traits::RangeFieldType;
            using WorldVector = Dune::FieldVector<RangeFieldType,CG::dow>;
            std::vector<WorldVector> gradients;

            for (auto const& qp : quad) {
                // Position of the current quadrature point in the reference element
                auto&& local = contextGeo.coordinateInElement(qp.position());

                // The transposed inverse Jacobian of the map from the reference element to the element
                const auto jacobian = contextGeo.elementGeometry().jacobianInverseTransposed(local);

                // The multiplicative factor in the integral transformation formula
                const auto factor = contextGeo.integrationElement(qp.position()) * qp.weight();
                const auto exprValue = localFct(local);

                // The gradients of the shape functions on the reference element
                auto const& shapeGradients = rowNode.localBasisJacobiansAt(local);

                // Compute the shape function gradients on the real element
                gradients.resize(shapeGradients.size());

                for (std::size_t i = 0; i < gradients.size(); ++i)
                    jacobian.mv(shapeGradients[i][0], gradients[i]);

                for (std::size_t i = 0; i < size; ++i) {
                    const auto local_i = rowNode.localIndex(i);
                    for (std::size_t j = 0; j < size; ++j) {
                        const auto local_j = colNode.localIndex(j);
                        if (contextGeo.elementGeometry().corner(i)[1] < 1.e-8) {
                            elementMatrix[local_i][local_j] += eval(exprValue, 4.0/3.0*factor, gradients[i], gradients[j]);
                        } else {
                            elementMatrix[local_i][local_j] += eval(exprValue, factor, gradients[i], gradients[j]);
                        }
                    }
                }
            }
        }

        template <class CG, class QR, class RN, class CN, class LocalFct, class Mat>
        void getElementMatrixOptimized(CG const& contextGeo, QR const& quad,
                                       RN const& node, CN const& /*colNode*/,
                                       LocalFct const& localFct, Mat& elementMatrix, tag::scalar) const
        {
            std::size_t size = node.size();

            using RangeFieldType = typename RN::LocalBasis::Traits::RangeFieldType;
            using WorldVector = Dune::FieldVector<RangeFieldType,CG::dow>;
            std::vector<WorldVector> gradients;

            for (auto const& qp : quad) {
                // Position of the current quadrature point in the reference element
                auto&& local = contextGeo.coordinateInElement(qp.position());

                // The transposed inverse Jacobian of the map from the reference element to the element
                const auto jacobian = contextGeo.elementGeometry().jacobianInverseTransposed(local);

                // The multiplicative factor in the integral transformation formula
                const auto factor = localFct(local) * contextGeo.integrationElement(qp.position()) * qp.weight();

                // The gradients of the shape functions on the reference element
                auto const& shapeGradients = node.localBasisJacobiansAt(local);

                // Compute the shape function gradients on the real element
                gradients.resize(shapeGradients.size());
                for (std::size_t i = 0; i < gradients.size(); ++i)
                    jacobian.mv(shapeGradients[i][0], gradients[i]);

                for (std::size_t i = 0; i < size; ++i) {
                    const auto local_i = node.localIndex(i);

                    if (contextGeo.elementGeometry().corner(i)[1] < 1.e-8) {
                        elementMatrix[local_i][local_i] += 3.0/4.0*factor * (gradients[i] * gradients[i]);
                    } else {
                        elementMatrix[local_i][local_i] += factor * (gradients[i] * gradients[i]);
                    }
                    for (std::size_t j = i+1; j < size; ++j) {
                        const auto local_j = node.localIndex(j);
                        const auto value = factor * (gradients[i] * gradients[j]);

                        if (contextGeo.elementGeometry().corner(i)[1] < 1.e-8) {
                            elementMatrix[local_i][local_j] += 3.0/4.0*value;
                        } else {
                            elementMatrix[local_i][local_j] += value;
                        }

                        if (contextGeo.elementGeometry().corner(j)[1] < 1.e-8) {
                            elementMatrix[local_j][local_i] += 3.0/4.0*value;
                        } else {
                            elementMatrix[local_j][local_i] += value;
                        }
                    }
                }
            }
        }

        template <class CG, class QR, class RN, class CN, class LocalFct, class Mat>
        void getElementMatrixOptimized(CG const& contextGeo, QR const& quad,
                                       RN const& node, CN const& /*colNode*/,
                                       LocalFct const& localFct, Mat& elementMatrix, tag::matrix) const
        {
            std::size_t size = node.size();

            using RangeFieldType = typename RN::LocalBasis::Traits::RangeFieldType;
            using WorldVector = Dune::FieldVector<RangeFieldType,CG::dow>;
            std::vector<WorldVector> gradients;

            for (auto const& qp : quad) {
                // Position of the current quadrature point in the reference element
                auto&& local = contextGeo.coordinateInElement(qp.position());

                // The transposed inverse Jacobian of the map from the reference element to the element
                const auto jacobian = contextGeo.elementGeometry().jacobianInverseTransposed(local);

                // The multiplicative factor in the integral transformation formula
                const auto factor = contextGeo.integrationElement(qp.position()) * qp.weight();
                const auto exprValue = localFct(local);

                // The gradients of the shape functions on the reference element
                auto const& shapeGradients = node.localBasisJacobiansAt(local);

                // Compute the shape function gradients on the real element
                gradients.resize(shapeGradients.size());
                for (std::size_t i = 0; i < gradients.size(); ++i)
                    jacobian.mv(shapeGradients[i][0], gradients[i]);

                for (std::size_t i = 0; i < size; ++i) {
                    const auto local_i = node.localIndex(i);
                    for (std::size_t j = 0; j < size; ++j) {
                        const auto local_j = node.localIndex(j);
                        if (contextGeo.elementGeometry().corner(i)[1] < 1.e-8) {
                            elementMatrix[local_i][local_j] += eval(exprValue, 3.0/4.0*factor, gradients[i],
                                                                    gradients[j]);
                        } else {
                            elementMatrix[local_i][local_j] += eval(exprValue, factor, gradients[i],
                                                                    gradients[j]);
                        }
                    }
                }
            }
        }

    protected:

        template <class S, class F, class T, int dow,
                std::enable_if_t<Category::Scalar<S>,int> = 0>
        T eval(S const& scalar, F factor,
               Dune::FieldVector<T,dow> const& grad_test,
               Dune::FieldVector<T,dow> const& grad_trial) const
        {
            return (scalar * factor) * (grad_test * grad_trial);
        }

        template <class M, class F, class T, int dow,
                std::enable_if_t<Category::Matrix<M>,int> = 0>
        T eval(M const& mat, F factor,
               Dune::FieldVector<T,dow> const& grad_test,
               Dune::FieldVector<T,dow> const& grad_trial) const
        {
            return factor * (grad_test * (mat * grad_trial));
        }
    };

    template <class LC>
    struct GridFunctionOperatorRegistry<tag::gradtest_gradtrial_axi, LC>
    {
        static constexpr int degree = 2;
        using type = SecondOrderGradTestGradTrialAxi;
    };

    /// Create a second-order term
    template <class Expr>
    auto sot_axi(Expr&& expr, int quadOrder = -1)
    {
        return makeOperator(tag::gradtest_gradtrial_axi{}, FWD(expr), quadOrder);
    }

    /** @} **/

} // end namespace AMDiS
