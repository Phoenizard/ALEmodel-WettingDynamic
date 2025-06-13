#pragma once

#include <type_traits>

#include <amdis/GridFunctionOperator.hpp>
#include <amdis/Output.hpp>
#include <amdis/common/StaticSize.hpp>
#include <amdis/common/ValueCategory.hpp>
#include <amdis/typetree/FiniteElementType.hpp>
#include <amdis/ALEModelAmdis2/myIndices.hpp>

using namespace MyIndices;
namespace AMDiS
{
  /**
   * \addtogroup operators
   * @{
   **/

  namespace tag
  {
      template <class DF>
      struct cahn_hilliard {
          unsigned int axi = 0; // axisymmetric or not?
          double sigma; // the phase field tension
          double eps; // the phase field interface width
          double const& mobility; // the constant mobility
          double const& invTau; // the inverse of the time step size
          double const& volumeDifference; //for additional droplet volume conservation
          int couplingVelocity = 0; // do you want to add coupling terms to a NS equation?
          FieldVector<double,2> vGrid; // the average grid velocity (to keep the shell centered in ALE simulations)
          DF const& laplaceSolution = {}; // the space dependant grid velocity (for ALE simulations with moving grid)
      };
  }


  /// Operator for the classic Cahn-Hilliard equation with constant mobility and double well
  /// W(\phi)=0.25\phi^(1-\phi)^2
  /// Coupling terms to a Navier-Stokes equation are added.
  /// ALE Operators, subtracting a grid velocity laplaceSolution_ from all convective terms, are also included.
  template <class DF>
  class CahnHilliardOperator {
      unsigned int axi_ = 0;
      double sigma_;
      double eps_;
      double const& mobility_;
      double const& invTau_;
      double const& volumeDifference_;
      int couplingVelocity_ = 0;
      FieldVector<double,2> vGrid_;
      DF const& laplaceSolution_ = {};
      double sigma_nucleus_ = 0;
      double center = -1.625;

  public:
      CahnHilliardOperator(tag::cahn_hilliard<DF> t)
      : axi_(t.axi), sigma_(t.sigma), eps_(t.eps), mobility_(t.mobility), invTau_(t.invTau),
        volumeDifference_(t.volumeDifference), couplingVelocity_(t.couplingVelocity), vGrid_(t.vGrid),
        laplaceSolution_(t.laplaceSolution)
      {
          sigma_nucleus_ = 6.0*sqrt(2)*Parameters::get<double>("parameters->sigma_1").value_or(t.sigma);
          center = Parameters::get<double>("nucleus center").value_or(-1.625);
      }

      template<class CG, class Node, class Quad, class LocalFct, class Mat>
      void assemble(CG const &contextGeo, Node const &tree, Node const /*colNode*/,
                    Quad const &quad, LocalFct const &localFct, Mat &elementMatrix) const {
          using expr_value_type = typename LocalFct::Range;
          static_assert(static_size_v<expr_value_type> == 1,
                        "Expression must be of scalar type.");
          using VelocityNode = Dune::TypeTree::Child<Node, _v.value>;
          using PhiNode = Dune::TypeTree::Child<Node, _phi.value>;
          using MuNode = Dune::TypeTree::Child<Node, _mu.value>;


          static_assert(PhiNode::isLeaf && MuNode::isLeaf && VelocityNode::isPower,
                        "Operator can be applied to correct nodes only.");
          using namespace Dune::Indices;

          auto const &vNode = tree.child(_v);
          auto const &phiNode = tree.child(_phi);
          auto const &muNode = tree.child(_mu);

          std::size_t vSize = vNode.child(0).size();
          std::size_t phiSize = phiNode.size();
          std::size_t muSize = muNode.size();

          using RangeFieldType = typename PhiNode::LocalBasis::Traits::RangeFieldType;
          using WorldVector = Dune::FieldVector<RangeFieldType, CG::dow>;
          std::vector<WorldVector> phiGradients;
          using RangeFieldTypeCol = typename MuNode::LocalBasis::Traits::RangeFieldType;
          using WorldVectorCol = Dune::FieldVector<RangeFieldTypeCol, CG::dow>;
          std::vector<WorldVectorCol> muGradients;

          auto gradPhi = derivativeOf(localFct, tag::gradient{});
          gradPhi.bind(contextGeo.element());

          auto gfU = makeGridFunction(valueOf(*laplaceSolution_), laplaceSolution_->basis().gridView());
          auto u = localFunction(gfU);
          u.bind(contextGeo.element());

          for (auto const &qp: quad) {
              // Position of the current quadrature point in the reference element
              auto &&local = contextGeo.coordinateInElement(qp.position());
              auto &&global = contextGeo.elementGeometry().global(local);

              // The transposed inverse Jacobian of the map from the reference element to the element
              const auto jacobian = contextGeo.elementGeometry().jacobianInverseTransposed(local);

              // The multiplicative factor in the integral transformation formula
              const auto factor = contextGeo.integrationElement(qp.position()) * qp.weight();
              const auto phi = localFct(local);

              // The gradients of the shape functions on the reference element
              auto const &phiShapeGradients = phiNode.localBasisJacobiansAt(local);
              auto const &muShapeGradients = muNode.localBasisJacobiansAt(local);

              // Compute the shape function gradients on the real element
              phiGradients.resize(phiShapeGradients.size());
              muGradients.resize(muShapeGradients.size());

              // assuming that phi and mu have the same polynomial degree
              for (std::size_t i = 0; i < phiGradients.size(); ++i) {
                  jacobian.mv(phiShapeGradients[i][0], phiGradients[i]);
                  jacobian.mv(muShapeGradients[i][0], muGradients[i]);
              }

              const auto factorDoubleWell = (Dune::at(global,0) > center ? sigma_ : sigma_nucleus_) / eps_ * (-0.5 - 3.0 * phi * (phi - 1.0)) * factor;
              const auto aleFactor = vGrid_ - u(local);

              auto const &vValues = vNode.child(0).localBasisValuesAt(local);
              auto const &phiValues = phiNode.localBasisValuesAt(local);
              auto const &muValues = muNode.localBasisValuesAt(local);

              for (std::size_t i = 0; i < phiSize; ++i) {
                  const auto local_i = phiNode.localIndex(i);
                  const auto local_ii = muNode.localIndex(i);
                  const auto valuePhiTime = invTau_ * factor * phiValues[i];
                  const auto valueMu = factor * muValues[i];
                  const auto valueDoubleWell = factorDoubleWell * phiValues[i];
                  double valueALE = aleFactor.dot(phiGradients[i]) * factor;

                  elementMatrix[local_i][local_i] += valuePhiTime * phiValues[i]; // lhs of d_t phi
                  elementMatrix[local_ii][local_ii] += valueMu * muValues[i]; // mu = ...

                  for (std::size_t j = i + 1; j < phiSize; ++j) {
                      const auto local_j = phiNode.localIndex(j);
                      const auto local_jj = muNode.localIndex(j);

                      elementMatrix[local_i][local_j] += valuePhiTime * phiValues[j]; // lhs of d_t phi
                      elementMatrix[local_j][local_i] += valuePhiTime * phiValues[j]; // lhs of d_t phi
                      elementMatrix[local_ii][local_jj] += valueMu * muValues[j]; // mu = ...
                      elementMatrix[local_jj][local_ii] += valueMu * muValues[j]; // mu = ...
                  }
                  
                  for (std::size_t j = 0; j < muSize; ++j) {
                      const auto local_j = muNode.localIndex(j);
                      elementMatrix[local_i][local_j] += eval(mobility_, factor, phiGradients[i],
                                                              muGradients[j]); // sot(mobility) (_phi, _mu)
                      elementMatrix[local_j][local_i] += eval(-(Dune::at(global,0) > center ? sigma_ : sigma_nucleus_) * eps_, factor, muGradients[j],
                                                              phiGradients[i]) // sot(-sigma*eps) (_mu, _phi)
                                                         + valueDoubleWell * muValues[j]; // W'(phi) implicit part (linearized)
                      
                      // axisymmetric terms
                      if (axi_) {
                          const auto valuePhiAxi =
                                  -mobility_ * factor / Dune::at(global,1) * muGradients[j][1];
                          elementMatrix[local_i][local_j] += valuePhiAxi * phiValues[i];

                          const auto valueMuAxi = eps_ * (Dune::at(global,0) > center ? sigma_ : sigma_nucleus_) * factor / Dune::at(global,1) * phiGradients[i][1];
                          elementMatrix[local_j][local_i] += valueMuAxi * muValues[j];
                      }
                  }

                  if (couplingVelocity_) {
                      auto dPhi = gradPhi(local);
                      // coupling term to velocity
                      for (std::size_t j = 0; j < vSize; ++j) {
                          const auto value = factor * phiValues[i] * vValues[j];

                          for (std::size_t k = 0; k < CG::dow; ++k) {
                              const auto local_kj = vNode.child(k).localIndex(j);
                              elementMatrix[local_i][local_kj] += value * Dune::at(dPhi, k);
                          }
                      }
                  }

                  // ALE operator
                  for (std::size_t j = 0; j < phiSize; ++j) {
                      const auto local_j = phiNode.localIndex(j);
                      elementMatrix[local_j][local_i] += valueALE * phiValues[j];
                  }
              }
          }
          gradPhi.unbind();
          u.unbind();
      }

      // assemble method for RHS
      template <class CG, class Node, class Quad, class LocalFct, class Vec>
      void assemble(CG const& contextGeo, Node const& tree, Quad const& quad,
                    LocalFct const& localFct, Vec& elementVector) const
      {
          static_assert(static_size_v<typename LocalFct::Range> == 1,
                        "Expression must be of scalar type." );

          using PhiNode = Dune::TypeTree::Child<Node,_phi.value>;
          using MuNode = Dune::TypeTree::Child<Node,_mu.value>;

          static_assert(PhiNode::isLeaf && MuNode::isLeaf, "");

          auto const& phiNode = tree.child(_phi);
          auto const& muNode = tree.child(_mu);

          std::size_t phiSize = phiNode.size();
          std::size_t muSize = muNode.size();

          auto gradPhi = derivativeOf(localFct,tag::gradient{});
          gradPhi.bind(contextGeo.element());

          for (auto const& qp : quad) {
              // Position of the current quadrature point in the reference element
              auto &&local = contextGeo.coordinateInElement(qp.position());
              auto &&global = contextGeo.elementGeometry().global(local);

              // The multiplicative factor in the integral transformation formula
              const auto factor = contextGeo.integrationElement(qp.position()) * qp.weight();
              const auto phi = localFct(local);
              auto phiFactor = (invTau_ * phi) * factor;
              if (std::abs(volumeDifference_) > 1.e-10)
                  phiFactor += volumeDifference_ * two_norm(gradPhi(local)) * factor;

              const auto muFactor = ((Dune::at(global,0) > center ? sigma_ : sigma_nucleus_) / eps_) * phi * phi * (1.5 - 2.0 * phi) * factor;

              auto const &shapeValuesPhi = phiNode.localBasisValuesAt(local);
              auto const &shapeValuesMu = muNode.localBasisValuesAt(local);

              // assuming phi and mu both have the same polynomial degree
              for (std::size_t i = 0; i < phiSize; ++i) {
                  const auto local_i = phiNode.localIndex(i);
                  const auto local_ii = muNode.localIndex(i);
                  elementVector[local_i] += phiFactor * shapeValuesPhi[i];
                  elementVector[local_ii] += muFactor * shapeValuesMu[i];
              }
          }
          gradPhi.unbind();
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

  template <class LC, class DF>
  struct GridFunctionOperatorRegistry<tag::cahn_hilliard<DF>, LC>
  {
    static constexpr int degree = 2;
    using type = CahnHilliardOperator<DF>;
  };

    /** @} **/
} // end namespace AMDiS
