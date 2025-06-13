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
      template <class DF, class DF2, class DF3>
      struct navier_stokes {
          unsigned int axi = 0; // axisymmetric or not?
          double rho; // fluid density (constant)
          DF const& nuPhase; // fluid viscosity (can be non-constant)
          double const& invTau; // inverse of the time step size
          double const& volumeDifference = 0; //for additional droplet volume conservation
          int constViscosity = 0; // is viscosity constant?
          int couplingPhaseField = 0; // do you want to add coupling terms?
          DF2 const& phi = {}; // the phase field
          double const& sigma = 0; // the phase field tension
          double const& eps = 0; // the phase field interface width
          FieldVector<double,2> const& vGrid; // the average grid velocity (to keep the shell centered in ALE simulations)
          DF3 const& laplaceSolution; // the space dependant grid velocity (for ALE simulations with moving grid)
      };
  }


  /// Operator for the classic Navier-Stokes equation with constant density. Coupling terms to a Cahn-Hilliard equation
  /// are added, including the possibility to use phase-dependant viscosity. The phase field is coupled to the NS
  /// equations by a capillary stress S_{\sigma_f} = \sigma\eps (\nabla\phi\otimes\nabla\phi)
  /// ALE Operators, subtracting a grid velocity laplaceSolution_ from all convective terms, are also included.
  template <class DF, class DF2, class DF3>
  class NavierStokesOperator {
      unsigned int axi_ = 0;
      double rho_;
      DF const& nuPhase_;
      double const& invTau_;
      double const& volumeDifference_;
      int constViscosity_;
      int couplingPhaseField_ = 0;
      DF2 const& phi_ = {};
      double const& sigma_ = 0;
      double const& eps_ = 0;
      FieldVector<double,2> const& vGrid_;
      DF3 const& laplaceSolution_;
      double sigma_nucleus_ = 0;
      double center = -1.625;

  public:
      NavierStokesOperator(tag::navier_stokes<DF,DF2,DF3> t)
      : axi_(t.axi), rho_(t.rho), nuPhase_(FWD(t.nuPhase)), invTau_(t.invTau), volumeDifference_(t.volumeDifference),
        constViscosity_(t.constViscosity), couplingPhaseField_(t.couplingPhaseField),
        phi_(FWD(t.phi)), sigma_(t.sigma), eps_(t.eps), vGrid_(t.vGrid), laplaceSolution_(t.laplaceSolution)
        {
          sigma_nucleus_ = 6.0*sqrt(2)*Parameters::get<double>("parameters->sigma_1").value_or(t.sigma);
          center = Parameters::get<double>("nucleus center").value_or(-1.625);
        }

      template<class CG, class Node, class Quad, class LocalFct, class Mat>
      void assemble(CG const &contextGeo, Node const &tree, Node const /*colNode*/,
                    Quad const &quad, LocalFct const &localFct, Mat &elementMatrix) const {
          using expr_value_type = typename LocalFct::Range;
          static_assert(static_size_v<expr_value_type> == CG::dow,
                        "Expression must be of vector type.");
          using VelocityNode = Dune::TypeTree::Child<Node, _v.value>;
          using PNode = Dune::TypeTree::Child<Node, _p.value>;
          using PhiNode = Dune::TypeTree::Child<Node, _phi.value>;
          using MuNode = Dune::TypeTree::Child<Node, _mu.value>;


          static_assert(PhiNode::isLeaf &&
                        MuNode::isLeaf &&
                        PNode::isLeaf &&
                        VelocityNode::isPower,
                        "Operator can be applied to correct nodes only.");
          using namespace Dune::Indices;

          auto const &vNode = tree.child(_v);
          auto const &pNode = tree.child(_p);
          auto const &phiNode = tree.child(_phi);
          auto const &muNode = tree.child(_mu);

          std::size_t vSize = vNode.child(0).size();
          std::size_t pSize = pNode.size();
          std::size_t phiSize = phiNode.size();
          std::size_t muSize = muNode.size();

          using RangeFieldType = typename VelocityNode::ChildType::LocalBasis::Traits::RangeFieldType;
          using WorldVector = Dune::FieldVector<RangeFieldType,CG::dow>;
          std::vector<WorldVector> vGradients;

          using RangeFieldTypePhi = typename PhiNode::LocalBasis::Traits::RangeFieldType;
          using WorldVectorPhi = Dune::FieldVector<RangeFieldTypePhi, CG::dow>;
          std::vector<WorldVectorPhi> phiGradients;

          auto gfNu = makeGridFunction(valueOf(*nuPhase_), nuPhase_->basis().gridView());
          auto nu = localFunction(gfNu);
          nu.bind(contextGeo.element());

          auto gfPhi = makeGridFunction(valueOf(*phi_), phi_->basis().gridView());
          auto phi = localFunction(gfPhi);
          phi.bind(contextGeo.element());

          auto gfU = makeGridFunction(valueOf(*laplaceSolution_), laplaceSolution_->basis().gridView());
          auto u = localFunction(gfU);
          u.bind(contextGeo.element());

          auto gradPhi = derivativeOf(phi,tag::gradient{});

          for (auto const &qp: quad) {
              // Position of the current quadrature point in the reference element
              auto &&local = contextGeo.coordinateInElement(qp.position());
              auto &&global = contextGeo.elementGeometry().global(local);

              // The transposed inverse Jacobian of the map from the reference element to the element
              const auto jacobian = contextGeo.elementGeometry().jacobianInverseTransposed(local);

              // The multiplicative factor in the integral transformation formula
              const auto factor = contextGeo.integrationElement(qp.position()) * qp.weight();
              const auto v = localFct(local);
              const auto vNu = nu(local) * factor;
              const auto vRho = invTau_ * rho_ * factor;
              const auto couplingFactor = -eps_ * (Dune::at(global,0) > center ? sigma_ : sigma_nucleus_) * factor;
              const auto grdPhi = gradPhi(local);
              const auto aleFactor = vGrid_ - u(local);

              // The gradients of the shape functions on the reference element
              auto const &vShapeGradients = vNode.child(0).localBasisJacobiansAt(local);
              auto const &phiShapeGradients = phiNode.localBasisJacobiansAt(local);

              // Compute the shape function gradients on the real element
              vGradients.resize(vShapeGradients.size());
              phiGradients.resize(phiShapeGradients.size());

              for (std::size_t i = 0; i < vGradients.size(); ++i)
                  jacobian.mv(vShapeGradients[i][0], vGradients[i]);

              for (std::size_t i = 0; i < phiGradients.size(); ++i)
                  jacobian.mv(phiShapeGradients[i][0], phiGradients[i]);

              auto const &vValues = vNode.child(0).localBasisValuesAt(local);
              auto const &pValues = pNode.localBasisValuesAt(local);
              auto const &phiValues = phiNode.localBasisValuesAt(local);

              // <viscosity * grad(u), grad(v)>
              // + <u^new*rho/tau,v> (time derivative impl part)
              // + <u^old grad u,v> (advective term)
              for (std::size_t i = 0; i < vSize; ++i) {
                  const auto value_ii = vNu * unary_dot(vGradients[i]);
                  const auto value_dt = vRho * vValues[i];
                  const auto value_ad = rho_ * factor * v.dot(vGradients[i]);
                  const auto value_al = rho_ * factor * aleFactor.dot(vGradients[i]);

                  for (std::size_t k = 0; k < CG::dow; ++k) {
                      const auto local_ki = vNode.child(k).localIndex(i);
                      elementMatrix[local_ki][local_ki] += value_ii + value_dt * vValues[i];

                      for (std::size_t j = 0; j < vSize; ++j) {
                          const auto local_kj = vNode.child(k).localIndex(j);
                          elementMatrix[local_kj][local_ki] += value_ad * vValues[j];

                          // ALE Operator
                          elementMatrix[local_kj][local_ki] += value_al * vValues[j];
                      }
                  }

                  for (std::size_t j = i+1; j < vSize; ++j) {
                      const auto value = vNu * vGradients[i].dot(vGradients[j]) + value_dt * vValues[j];

                      for (std::size_t k = 0; k < CG::dow; ++k) {
                          const auto local_ki = vNode.child(k).localIndex(i);
                          const auto local_kj = vNode.child(k).localIndex(j);
                          elementMatrix[local_ki][local_kj] += value;
                          elementMatrix[local_kj][local_ki] += value;
                      }
                  }
              }

              if (!constViscosity_) {
                  // <viscosity * d_i u_j, d_j v_i>
                  for (std::size_t i = 0; i < vSize; ++i) {
                      for (std::size_t kj = 0; kj < CG::dow; ++kj) {
                          const auto value_i_kj = vNu * vGradients[i][kj];
                          for (std::size_t j = 0; j < vSize; ++j) {
                              for (std::size_t ki = 0; ki < CG::dow; ++ki) {
                                  const auto local_ki = vNode.child(ki).localIndex(i);
                                  const auto local_kj = vNode.child(kj).localIndex(j);
                                  elementMatrix[local_ki][local_kj] += value_i_kj * vGradients[j][ki];
                              }
                          }
                      }
                  }
              }

              // <p, div(v)> + <div(u), q>
              for (std::size_t i = 0; i < vSize; ++i) {
                  for (std::size_t j = 0; j < pSize; ++j) {
                      const auto value = factor * pValues[j];
                      for (std::size_t k = 0; k < CG::dow; ++k) {
                          const auto local_v = vNode.child(k).localIndex(i);
                          const auto local_p = pNode.localIndex(j);

                          elementMatrix[local_v][local_p] += vGradients[i][k] * value;
                          elementMatrix[local_p][local_v] += vGradients[i][k] * value;
                      }
                  }
              }

              //coupling to phase field: <grad(phi)\otimes grad(phi), grad(v)>
              if (couplingPhaseField_) {
                  for (std::size_t j = 0; j < phiSize; ++j) {
                      const auto local_j = phiNode.localIndex(j);
                      for (std::size_t k = 0; k < CG::dow; ++k) {
                          const auto value = couplingFactor * Dune::at(grdPhi,k);
                          for (std::size_t i = 0; i < vSize; ++i) {
                              const auto local_ki = vNode.child(k).localIndex(i);
                              for (std::size_t l = 0; l < CG::dow; ++l) {
                                  elementMatrix[local_ki][local_j] += value * vGradients[i][l] * phiGradients[j][l];
                              }
                          }
                      }
                  }
              }

              // axisymmetric terms
              if (axi_) {
                  // coupling terms to phase field
                  if (couplingPhaseField_) {
                      for (std::size_t j = 0; j < phiSize; ++j) {
                          const auto local_j = phiNode.localIndex(j);

                          for (std::size_t i = 0; i < vSize; ++i) {
                              for (std::size_t k = 0; k < CG::dow; ++k) {
                                  const auto value = -couplingFactor * Dune::at(grdPhi,1) / Dune::at(global,1) * phiGradients[j][k];
                                  const auto local_ki = vNode.child(k).localIndex(i);
                                  elementMatrix[local_ki][local_j] += value * vValues[i];
                              }
                          }
                      }
                  }
                  //axi terms due to NS Equation
                  for (std::size_t j = 0; j < vSize; ++j) {
                      const auto local_j0 = vNode.child(0).localIndex(j);
                      const auto local_j1 = vNode.child(1).localIndex(j);
                      const auto value00 = -vNu / Dune::at(global,1) * vGradients[j][1];
                      const auto value01 = -vNu / Dune::at(global,1) * vGradients[j][0];
                      const auto value11 = -2.0 * vNu / Dune::at(global,1) * vGradients[j][1];
                      const auto value2 = 2.0 * vNu / (Dune::at(global,1)*Dune::at(global,1)) * vValues[j];
                      const auto valueP = factor / Dune::at(global,1) * vValues[j];

                      for (std::size_t i = 0; i < vSize; ++i) {
                          const auto local_i0 = vNode.child(0).localIndex(i);
                          const auto local_i1 = vNode.child(1).localIndex(i);
                          elementMatrix[local_i0][local_j0] += value00 * vValues[i];
                          elementMatrix[local_i0][local_j1] += value01 * vValues[i];
                          elementMatrix[local_i1][local_j1] += (value11 + value2) * vValues[i];
                      }

                      for (std::size_t i = 0; i < pSize; ++i) {
                          const auto local_ip = pNode.localIndex(i);
                          elementMatrix[local_ip][local_j1] += valueP * pValues[i];
                      }
                  }

              }
          }
          nu.unbind();
          phi.unbind();
          u.unbind();
      }
      
      // assemble method for RHS
      template <class CG, class Node, class Quad, class LocalFct, class Vec>
      void assemble(CG const& contextGeo, Node const& tree, Quad const& quad,
                    LocalFct const& localFct, Vec& elementVector) const
      {
          static_assert(static_size_v<typename LocalFct::Range> == CG::dow,
                        "Expression must be of vector type." );

          using VNode = Dune::TypeTree::Child<Node,_v.value>;
          using PNode = Dune::TypeTree::Child<Node,_p.value>;

          static_assert(VNode::isPower && PNode::isLeaf, "");

          auto const& vNode = tree.child(_v);
          auto const& pNode = tree.child(_p);

          std::size_t vSize = vNode.child(0).size();
          std::size_t pSize = pNode.size();

          for (auto const& qp : quad) {
              // Position of the current quadrature point in the reference element
              auto &&local = contextGeo.coordinateInElement(qp.position());

              // The multiplicative factor in the integral transformation formula
              auto factor = contextGeo.integrationElement(qp.position()) * qp.weight();
              auto exprValue = localFct(local);

              auto const &shapeValues = vNode.child(0).localBasisValuesAt(local);

              for (std::size_t i = 0; i < vSize; ++i) {
                  const auto value = invTau_ * rho_ * (factor * shapeValues[i]);
                  for (std::size_t k = 0; k < CG::dow; ++k) {
                      const auto local_ki = vNode.child(k).localIndex(i);
                      elementVector[local_ki] += value * Dune::at(exprValue, k);
                  }
              }

              // volume conservation term
              auto const& pShapeValues = pNode.localBasisValuesAt(local);
              for (std::size_t i = 0; i < pSize; ++i) {
                  const auto local_i = pNode.localIndex(i);
                  elementVector[local_i] += volumeDifference_ * factor * pShapeValues[i];
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

  template <class LC, class DF, class DF2, class DF3>
  struct GridFunctionOperatorRegistry<tag::navier_stokes<DF, DF2, DF3>, LC>
  {
    static constexpr int degree = 2;
    using type = NavierStokesOperator<DF, DF2, DF3>;
  };

    /** @} **/
} // end namespace AMDiS
