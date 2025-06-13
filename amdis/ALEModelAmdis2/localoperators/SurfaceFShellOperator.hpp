#pragma once

#include <type_traits>

#include <amdis/GridFunctionOperator.hpp>
#include <amdis/common/StaticSize.hpp>

using namespace MyIndices;
namespace AMDiS {
    /**
     * \addtogroup operators
     * @{
     **/

    namespace tag {
        template<class DV, class DV2, class DV3, class IS = int>
        struct surface_fShell {
            unsigned int axi = 0; // axisymmetric simulations (=1) or classic (=0)
            double surfaceTension = 0; // constant surface tension force
            double Ka = 0; // area dilation modulus
            double Ks = 0; // area shear modulus
            double Kb = 0; // bending stiffness
            double sigma0 = 0; // surface tension of membrane in contact with ambient (phi = 0)
            double sigma1 = 0; // surface tension of membrane in contact with droplet (phi = 1)
            double kappa0 = 0; // spontaneous curvature of membrane
            DV const& normalVec = {}; // phase field solution from last time step
            DV2 const& oldKappaVec = {}; // curvature vector on current grid
            DV3 const& oldPhi = {}; // phase field solution from last time step
            int inside = 0; // integrate over the shell from inside?
            int outside = 1; // integrate over the shell from  outside?
            int incompressibleShell = 0; // is shell assumed strictly incompressible?
            double lineTension = 0; // line tension (only used if axi = 1)
            IS const &indexSet = {}; //the indexSet of the leafGridView of the grid!
            std::vector<int> const &partitions = {}; //partition numbers of each element
            std::vector<int> const &facets = {}; //facet number of the surface facet (-1 if element i not touching the surface)
        };
    }


    /// Sum of all forces on the shell, including surface tension, bending stiffness, stretching,
    /// force due to surface phase field and also a contact angle condition for the contact angle of a droplet on the
    /// shell.
    template<class DV, class DV2, class DV3, class IS = int>
    class SurfaceFShell {
        unsigned int axi_ = 0;
        double surfaceTension;
        DV const& normalVec_;
        DV2 const& oldKappaVec_;
        DV3 const& oldPhi_;
        double sigma0_ = 0;
        double sigma1_ = 0;
        double kappa0_ = 0;
        double Ka_ = 0;
        double Ks_ = 0;
        double Kb_ = 0;
        int inside_;
        int outside_;
        int incompressible_;
        double lineTension_;
        IS const &indexSet_ = {};
        std::vector<int> const &partitions_;
        std::vector<int> const &facets_;
        double sigma0_nucleus_ = 0;
        double sigma1_nucleus_ = 0;
        double center = -1.625;

    public:
        SurfaceFShell(tag::surface_fShell<DV,DV2,DV3,IS> t)
                : axi_(t.axi), surfaceTension(t.surfaceTension), normalVec_(FWD(t.normalVec)),
                  oldKappaVec_(FWD(t.oldKappaVec)), oldPhi_(FWD(t.oldPhi)), sigma0_(t.sigma0), sigma1_(t.sigma1), kappa0_(t.kappa0),
                  partitions_(t.partitions), Ka_(t.Ka), Ks_(t.Ks), Kb_(t.Kb), facets_(t.facets), indexSet_(t.indexSet),
                  inside_(t.inside), outside_(t.outside), incompressible_(t.incompressibleShell),
                  lineTension_(t.lineTension)
        {
            sigma0_nucleus_ = Parameters::get<double>("parameters->sigma0_1").value_or(t.sigma0);
            sigma1_nucleus_ = Parameters::get<double>("parameters->sigma1_1").value_or(t.sigma1);
            center = Parameters::get<double>("nucleus center").value_or(-1.e30);
        }

        template<class CG, class Node, class Quad, class LocalFct, class Mat>
        void assemble(CG const &contextGeo, Node const &tree, Node const /*colNode*/,
                      Quad const &quad, LocalFct const &localFct, Mat &elementMatrix) const {
            static_assert(static_size_v<typename LocalFct::Range> == 1,
                          "Expression must be of scalar type.");
            
            using FShellNode = Dune::TypeTree::Child<Node, _fShell.value>;
            using VNode = Dune::TypeTree::Child<Node, _v.value>;
            using PhiNode = Dune::TypeTree::Child<Node, _phi.value>;
            using MuNode = Dune::TypeTree::Child<Node, _mu.value>;
            using PhiSNode = Dune::TypeTree::Child<Node, _phiS.value>;
            using KappaVecNode = Dune::TypeTree::Child<Node, _kappaVec.value>;
            using Lambda1Node = Dune::TypeTree::Child<Node, _lambda1.value>;
            using Lambda2Node = Dune::TypeTree::Child<Node, _lambda2.value>;

            static_assert(FShellNode::isPower && VNode::isPower && KappaVecNode::isPower 
                           && PhiNode::isLeaf && PhiSNode::isLeaf && MuNode::isLeaf && Lambda1Node::isLeaf
                          && Lambda2Node::isLeaf,
                          "Nodes must have correct dimensions.");
            
            auto const &vNode = tree.child(_v);
            auto const &fShellNode = tree.child(_fShell);
            auto const &phiNode = tree.child(_phi);
            auto const &phiSNode = tree.child(_phiS);
            auto const &muNode = tree.child(_mu);
            auto const &kappaVecNode = tree.child(_kappaVec);
            auto const &lambda1Node = tree.child(_lambda1);
            auto const &lambda2Node = tree.child(_lambda2);


            if (fShellNode.child(0).size() == 0)
                return;


            // compute facet number and partition number from surface basis fShellNode (if no surface basis fShellNode: return -2)
            auto part = partition(fShellNode.child(0), Dune::PriorityTag<2>{});
            auto face = facet(fShellNode.child(0), Dune::PriorityTag<2>{});

            if (face == -1) // do nothing if element has no facet on the shell
                return;

            std::size_t fShellSize = fShellNode.child(0).size();
            std::size_t vSize = vNode.child(0).size();
            std::size_t phiSize = phiNode.size();
            std::size_t muSize = muNode.size();
            std::size_t phiS_Size = phiSNode.size();
            std::size_t kappaVecSize = kappaVecNode.child(0).size();
            std::size_t lambda1Size = lambda1Node.size();
            std::size_t lambda2Size = lambda2Node.size();

            using RangeFieldType = typename FShellNode::ChildType::LocalBasis::Traits::RangeFieldType;
            using WorldVector = Dune::FieldVector<RangeFieldType, CG::dow>;
            using WorldMatrix = Dune::FieldMatrix<RangeFieldType,CG::dow,CG::dow>;
            std::vector<WorldVector> fShellGradients, lambda1Gradients, lambda2Gradients,
                                     kappaVecGradients, phiSGradients;

            auto element = contextGeo.context();
            auto intersection = element.template subEntity<1>(face);

            constexpr auto intersectionDim = std::decay_t<decltype(intersection)>::mydimension;
            const auto &quadShell = Dune::QuadratureRules<double, intersectionDim>::rule(intersection.type(), 2);

            // apply the operator only to intersections on the shell
            if (part == 1 && inside_ || part == 0 && outside_) {
                auto gfPhi = makeGridFunction(valueOf(*oldPhi_), oldPhi_->basis().gridView());
                auto phi = localFunction(gfPhi);
                phi.bind(contextGeo.element());

                auto normalX = localFunction(valueOf(normalVec_[0]));
                auto normalY = localFunction(valueOf(normalVec_[1]));
                normalX.bind(contextGeo.element());
                normalY.bind(contextGeo.element());
                auto geo = intersection.geometry();

                auto normal = geo.corner(0);
                for (auto const &qp: quadShell) {
                    auto localGeo = referenceElement(element).template geometry<1>(face);

                    // Position of the current quadrature point in the reference element
                    auto &&local = localGeo.global(qp.position());
                    auto &&global = contextGeo.elementGeometry().global(local);

                    normal[0] = normalX(local);
                    normal[1] = normalY(local);

                    // The transposed inverse Jacobian of the map from the reference facet to the real facet
                    const auto jacobian = geo.jacobianInverseTransposed(qp.position());

                    // The multiplicative factor in the integral transformation formula
                    auto dx = geo.integrationElement(qp.position()) * qp.weight();
                    auto dxAxi = (axi_ * Dune::at(global,1) + 1.0 - axi_) * dx;

                    // The gradients of the shape functions on the reference element
                    auto const &shapeGradientsFShell = fShellNode.child(0).localBasisJacobiansAt(local);
                    auto const &shapeGradientsPhiS = phiSNode.localBasisJacobiansAt(local);

                    fShellGradients.resize(shapeGradientsFShell.size());
                    phiSGradients.resize(shapeGradientsPhiS.size());

                    // transform local gradients to reference intersection if necessary and the compute
                    // gradients on the real intersection
                    int dow = CG::dow;

                    for (std::size_t i = 0; i < phiSGradients.size(); ++i)
                        jacobian.mv(shapeGradientsPhiS[i][0], phiSGradients[i]);


                    auto const &fShellValues = fShellNode.child(0).localBasisValuesAt(local);
                    auto const &phiValues = phiNode.localBasisValuesAt(local);
                    auto const &muValues = muNode.localBasisValuesAt(local);
                    auto const &kappaVecValues = kappaVecNode.child(0).localBasisValuesAt(local);

                    // assemble the following operators only from one side!
                    if (part == 0 || (!outside_ && part == 1)) {
                        auto gfKappaVec = makeGridFunction(valueOf(*oldKappaVec_), oldKappaVec_->basis().gridView());
                        auto kappaVec = localFunction(gfKappaVec);
                        kappaVec.bind(contextGeo.element());

                        auto const &shapeGradientsLambda1 = lambda1Node.localBasisJacobiansAt(local);
                        lambda1Gradients.resize(shapeGradientsLambda1.size());

                        auto const &shapeGradientsLambda2 = lambda2Node.localBasisJacobiansAt(local);
                        lambda2Gradients.resize(shapeGradientsLambda2.size());

                        auto const &shapeGradientsKappa = kappaVecNode.child(0).localBasisJacobiansAt(local);
                        kappaVecGradients.resize(shapeGradientsKappa.size());

                        for (std::size_t i = 0; i < fShellGradients.size(); ++i)
                            jacobian.mv(shapeGradientsFShell[i][0], fShellGradients[i]);

                        for (std::size_t i = 0; i < kappaVecGradients.size(); ++i)
                            jacobian.mv(shapeGradientsKappa[i][0], kappaVecGradients[i]);

                        auto kappaVecFactor = Dune::FieldVector<double,CG::dow>{0.0};
                        auto stretchingMarangoniFactor = 0.0;
                        auto kappaVecFactor2 = Dune::FieldVector<double,CG::dow>{0.0};
                        auto stretchingMarangoniFactor2 = 0.0;
                        auto kV = kappaVec(local);

                        if (!incompressible_) {
                            for (std::size_t i = 0; i < lambda1Gradients.size(); ++i)
                                jacobian.mv(shapeGradientsLambda1[i][0], lambda1Gradients[i]);

                            for (std::size_t i = 0; i < lambda2Gradients.size(); ++i)
                                jacobian.mv(shapeGradientsLambda2[i][0], lambda2Gradients[i]);
                           
                            for (std::size_t i = 0; i < CG::dow; ++i) {
                                kappaVecFactor[i] = -(Ka_ + Ks_)  * (Dune::at(kV,i) + kappa0_ * Dune::at(normal,i)) * dxAxi;
                            }
                            stretchingMarangoniFactor = -(Ka_ + Ks_) * dxAxi;
                            if (axi_ && Ks_) {
                                for (std::size_t i = 0; i < CG::dow; ++i) {
                                    kappaVecFactor2[i] = -(Ka_ - Ks_)  * (Dune::at(kV,i) + kappa0_ * Dune::at(normal,i)) * dxAxi;
                                }
                                stretchingMarangoniFactor2 = -(Ka_ - Ks_) * dxAxi;
                            }
                        }


                        // Projection Matrix
                        WorldMatrix P = Dune::ScaledIdentityMatrix<RangeFieldType,CG::dow>(1);
                        P -= Dune::outer(normal,normal);

                        auto const &vValues = vNode.child(0).localBasisValuesAt(local);
                        auto const &lambda1Values = lambda1Node.localBasisValuesAt(local);
                        auto const &lambda2Values = lambda2Node.localBasisValuesAt(local);

                        const auto bendingFactor = (-0.5 * Kb_ * dxAxi);
                        auto bendingFactorAxi = bendingFactor;
                        if (axi_) bendingFactorAxi = (-0.5 * Kb_ * dx);
                        const auto bendingFactor2D = (dow == 3 || axi_ ? -1.0 : 1.0) * Kb_ * dxAxi;
                        const auto spontFac = Kb_ * kappa0_ * dxAxi;
                        const auto spontFac2 = Kb_ * kappa0_ * dx;
                        const auto stretchKSFac = 2.0 * Ks_ * dx;
                        auto kbdx = 0.0;
                        if (axi_) kbdx = -Kb_ * dx;

                        // v = fShell on Gamma
                        for (std::size_t i = 0; i < fShellSize; ++i) {
                            const auto value = dxAxi * fShellValues[i];
                            const auto valueV = dx * fShellValues[i];

                            // fShell = ...
                            for (std::size_t k = 0; k < fShellNode.degree(); ++k) {
                                const auto local_ki = fShellNode.child(k).localIndex(i);
                                elementMatrix[local_ki][local_ki] += value * fShellValues[i];
                            }
                            for (std::size_t j = i + 1; j < fShellSize; ++j) {
                                for (std::size_t k = 0; k < fShellNode.degree(); ++k) {
                                    const auto local_ki = fShellNode.child(k).localIndex(i);
                                    const auto local_kj = fShellNode.child(k).localIndex(j);

                                    elementMatrix[local_ki][local_kj] += value * fShellValues[j];
                                    elementMatrix[local_kj][local_ki] += value * fShellValues[j];
                                }
                            }

                            // v = fShell on Gamma
                            auto ignoreLaplace = Parameters::get<int>("ignore laplace problem").value_or(0);
                            if (!ignoreLaplace) { //only couple force to NS system if the surface is not stiff
                                for (std::size_t j = 0; j < vSize; ++j) {
                                    for (std::size_t k = 0; k < CG::dow; ++k) {
                                        const auto local_kj = vNode.child(k).localIndex(j);
                                        const auto local_ki = fShellNode.child(k).localIndex(i);
                                        elementMatrix[local_kj][local_ki] += -valueV * vValues[j];
                                    }
                                }
                            }
                            // stretching force operators
                            if (!incompressible_) {
                                const auto valueStretchMar = stretchingMarangoniFactor * fShellValues[i];
                                const auto valueStretchMar2 = stretchingMarangoniFactor2 * fShellValues[i];
                                const auto axiFactor = stretchKSFac * fShellValues[i];

                                for (std::size_t j = 0; j < lambda1Size; ++j) {
                                    const auto stretchST = fShellValues[i] * lambda1Values[j];
                                    const auto stretchST2 = fShellValues[i] * lambda2Values[j];
                                    for (std::size_t k = 0; k < CG::dow; ++k) {
                                        const auto stretchMa = valueStretchMar * lambda1Gradients[j][k];
                                        const auto local_ki = fShellNode.child(k).localIndex(i);
                                        const auto local_j = lambda1Node.localIndex(j);

                                        elementMatrix[local_ki][local_j] += stretchST * Dune::at(kappaVecFactor,k) + stretchMa;
                                        if (axi_ && Ks_) {
                                            const auto local_jj = lambda2Node.localIndex(j);
                                            const auto stretchMa2 = valueStretchMar2 * lambda2Gradients[j][k];
                                            const auto axiValue = -axiFactor * lambda1Values[j];
                                            const auto axiValue2 = axiFactor * lambda2Values[j];

                                            elementMatrix[local_ki][local_jj] += stretchST2* Dune::at(kappaVecFactor2,k) + stretchMa2;
                                            if (k == 1) {
                                                elementMatrix[local_ki][local_j] += axiValue;
                                                elementMatrix[local_ki][local_jj] += axiValue2;
                                            }
                                        }
                                    }
                                }
                            }

                            // bending force in weak form
                            for (std::size_t l = 0; l < CG::dow; ++l) {
                                const auto local_li = fShellNode.child(l).localIndex(i);
                                for (std::size_t j = 0; j < kappaVecSize; ++j) {
                                    const auto local_lj = kappaVecNode.child(l).localIndex(j);

                                    // gradtest_gradtrial
                                    elementMatrix[local_li][local_lj] += eval(-Kb_, dxAxi, fShellGradients[i],
                                                                              kappaVecGradients[j]);
                                    if (axi_) {
                                        const auto local_1i = fShellNode.child(1).localIndex(i);
                                        const auto local_1j = kappaVecNode.child(1).localIndex(j);
                                        elementMatrix[local_1i][local_lj] +=
                                                fShellValues[i] * bendingFactorAxi * Dune::at(kV,l) * kappaVecValues[j]
                                                + kbdx * fShellValues[i] * kappaVecGradients[j][l];

                                        elementMatrix[local_li][local_1j] +=
                                                kbdx * fShellGradients[i][l] * kappaVecValues[j];

                                        if (kappa0_ && axi_) {
                                            const auto local_1i = fShellNode.child(1).localIndex(i);
                                            elementMatrix[local_1i][local_lj] += -spontFac2 * Dune::at(normal,l) * fShellValues[i] * kappaVecValues[j];
                                        }
                                    }
                                    for (std::size_t k = 0; k < CG::dow; ++k) {
                                        const auto local_kj = kappaVecNode.child(k).localIndex(j);
                                        // divtestvec_norm(kappaVec) + divtestvec_divtrialvec
                                        // (including distinction between 2D and 3D/axi case)
                                        elementMatrix[local_li][local_kj] +=
                                                  bendingFactor * Dune::at(kV,k) * fShellGradients[i][l] * kappaVecValues[j]
                                                + bendingFactor2D * fShellGradients[i][l] * kappaVecGradients[j][k];

                                        // rate of deformation (only in axi or 3D case)
                                        if (dow == 3 || axi_) {
                                            elementMatrix[local_li][local_kj] +=
                                                    -bendingFactor2D * ((P[k][l] * fShellGradients[i] * kappaVecGradients[j])
                                                                        + fShellGradients[i][k] * kappaVecGradients[j][l]);
                                        }

                                        // correction if spontaneous curvature is nonzero
                                        if (kappa0_) {
                                            elementMatrix[local_li][local_kj] +=
                                                    -spontFac * Dune::at(normal,k) * fShellGradients[i][l] * kappaVecValues[j]
                                                    +spontFac * Dune::at(normal,l) * fShellGradients[i][k] * kappaVecValues[j];
                                        }
                                    }
                                }
                            }
                        }
                        kappaVec.unbind();
                    }
                    
                    // assemble the following operators from outside, if outside_ == 1, and from inside if inside_ == 1
                    auto phiVal = phi(local);

                    /*
                    // contact angle condition due to phase field in contact with the surface
                    auto contactFactor = -6.0 * (1.0 - phiVal) * (sigma1_ - sigma0_) * dx;
                    for (std::size_t i = 0; i < muSize; ++i) {
                        const auto local_i = muNode.localIndex(i);
                        const auto value = contactFactor * muValues[i];

                        for (std::size_t j = 0; j < phiSize; ++j) {
                            const auto local_j = phiNode.localIndex(j);
                            elementMatrix[local_i][local_j] += value * phiValues[j];
                        }
                    }
*/
                    if (!incompressible_) {
                        // surface tension due to phase field on the shell
                        auto surfTenFactor = -(surfaceTension + ((Dune::at(global,0) > center ? sigma1_ - sigma0_ : sigma1_nucleus_ - sigma0_nucleus_) * (phiVal * phiVal) *
                                               (3.0 - 2.0 * phiVal) + (Dune::at(global,0) > center ? sigma0_ : sigma0_nucleus_))) * dxAxi;
                        for (std::size_t i = 0; i < fShellSize; ++i) {
                            const auto value = surfTenFactor * fShellValues[i];

                            for (std::size_t j = 0; j < kappaVecSize; ++j) {
                                const auto value0 = value * kappaVecValues[j];
                                for (std::size_t k = 0; k < CG::dow; ++k) {
                                    const auto local_ki = fShellNode.child(k).localIndex(i);
                                    const auto local_kj = kappaVecNode.child(k).localIndex(j);

                                    elementMatrix[local_ki][local_kj] += value0;
                                }
                            }

                            // marangoni term of surface tension due to phase field
                            auto gradSigmaSFactor = -6.0 * phiVal * (1.0 - phiVal) * (Dune::at(global,0) > center ? sigma1_ - sigma0_ : sigma1_nucleus_ - sigma0_nucleus_)
                                                  * dxAxi * fShellValues[i];
                            for (std::size_t j = 0; j < phiS_Size; ++j) {
                                const auto local_j = phiSNode.localIndex(j);
                                for (std::size_t k = 0; k < dow; ++k) {
                                    const auto local_ki = fShellNode.child(k).localIndex(i);
                                    elementMatrix[local_ki][local_j] += gradSigmaSFactor * phiSGradients[j][k];
                                }
                            }
                        }
                    }
                }
                normalX.unbind();
                normalY.unbind();
                phi.unbind();
            }
        }
        
        // assemble method for RHS
        template<class CG, class Node, class Quad, class LF, class Vec>
        void assemble(CG const &contextGeo, Node const &tree, Quad const &quad,
                      LF const &localFct, Vec &elementVector) const {
            static_assert(static_size_v<typename LF::Range> == 1,
                          "Expression must be of scalar type.");

            using FShellNode = Dune::TypeTree::Child<Node, _fShell.value>;
            using MuNode = Dune::TypeTree::Child<Node, _mu.value>;
            using PhiNode = Dune::TypeTree::Child<Node, _phi.value>;

            static_assert(FShellNode::isPower && MuNode::isLeaf && PhiNode::isLeaf,
                          "Nodes must have correct dimensions.");

            auto const &fShellNode = tree.child(_fShell);
            auto const &muNode = tree.child(_mu);
            auto const &phiNode = tree.child(_phi);

            if (fShellNode.child(0).size() == 0)
                return;

            // compute facet number and partition number from surface basis fShellNode (if no surface basis fShellNode: return -2)
            auto part = partition(fShellNode.child(0), Dune::PriorityTag<2>{});
            auto face = facet(fShellNode.child(0), Dune::PriorityTag<2>{});

            if (face == -1) // do nothing if element has no facet on the shell
                return;

            std::size_t fShellSize = fShellNode.child(0).size();
            std::size_t muSize = muNode.size();
            std::size_t phiSize = phiNode.size();

            using RangeFieldType = typename FShellNode::ChildType::LocalBasis::Traits::RangeFieldType;
            using WorldVector = FieldVector<RangeFieldType, CG::dow>;
            using WorldMatrix = Dune::FieldMatrix<RangeFieldType,CG::dow,CG::dow>;
            std::vector<WorldVector> gradients;

            auto element = contextGeo.context();

            auto intersection = element.template subEntity<1>(face);

            constexpr auto intersectionDim = std::decay_t<decltype(intersection)>::mydimension;
            const auto &quadShell = Dune::QuadratureRules<double, intersectionDim>::rule(intersection.type(), 2);

            assert(inside_ != outside_); //not implemented, use classical operator instead

            // apply the operator only to intersections on the shell
            if (part == 1 && inside_ || part == 0 && outside_) {
                auto gfPhi = makeGridFunction(valueOf(*oldPhi_), oldPhi_->basis().gridView());
                auto phi = localFunction(gfPhi);
                phi.bind(contextGeo.element());
                auto gradPhi = derivativeOf(phi, tag::gradient{});
                gradPhi.bind(contextGeo.element());

                auto normalX = localFunction(valueOf(normalVec_[0]));
                auto normalY = localFunction(valueOf(normalVec_[1]));
                normalX.bind(contextGeo.element());
                normalY.bind(contextGeo.element());

                auto geo = intersection.geometry();
                auto normal = geo.corner(0);

                for (auto const &qp: quadShell) {
                    auto localGeo = referenceElement(element).template geometry<1>(face);

                    // Position of the current quadrature point in the reference element
                    auto &&local = localGeo.global(qp.position());
                    auto &&global = contextGeo.elementGeometry().global(local);

                    // The transposed inverse Jacobian of the map from the reference facet to the real facet
                    const auto jacobian = geo.jacobianInverseTransposed(qp.position());


                    // The multiplicative factor in the integral transformation formula
                    const auto dx = geo.integrationElement(qp.position()) * qp.weight();
                    const auto dxAxi = (axi_ * Dune::at(global,1) + 1.0 - axi_) * dx;

                    // The gradients of the shape functions on the reference element
                    auto const &shapeGradients = fShellNode.child(0).localBasisJacobiansAt(local);

                    // Compute the shape function gradients on the real element
                    gradients.resize(shapeGradients.size());

                    // transform local gradients to reference intersection if necessary and compute
                    // gradients on the real intersection
                    int dow = CG::dow;

                    for (std::size_t i = 0; i < gradients.size(); ++i)
                        jacobian.mv(shapeGradients[i][0], gradients[i]);

                    normal[0] = normalX(local);
                    normal[1] = normalY(local);

                    double lineTensionFactor = 0.0;
                    if (lineTension_ && axi_) {
                        auto gradOld = gradPhi(local);
                        auto aux = gradOld;
                        WorldMatrix P = Dune::ScaledIdentityMatrix<RangeFieldType,CG::dow>(1);
                        P -= Dune::outer(normal,normal);
                        P.mv(aux,gradOld);

                        lineTensionFactor = lineTension_ * gradOld.two_norm() * dx;
                    }
                    double kappa0Factor = 0.0;
                    auto const &shapeValues = fShellNode.child(1).localBasisValuesAt(local);
                    auto const &muValues = muNode.localBasisValuesAt(local);
                    double spontVal = 0.0;
                    if (kappa0_) {
                        auto phiVal = phi(local);
                        spontVal = -Kb_ * kappa0_ * kappa0_ * dxAxi;
                        kappa0Factor = (kappa0_ * (surfaceTension + ((Dune::at(global,0) > center ? sigma1_ - sigma0_ : sigma1_nucleus_ - sigma0_nucleus_) * phiVal * phiVal *
                                                                      (3.0 - 2.0 * phiVal) + (Dune::at(global,0) > center ? sigma0_ : sigma0_nucleus_))) * dxAxi);
                    }
                    for (std::size_t j = 0; j < fShellSize; ++j) {
                        for (std::size_t k = 0; k < CG::dow; ++k) {
                            const auto local_kj = fShellNode.child(k).localIndex(j);
                            if (kappa0_) {
                                for (std::size_t l = 0; l < CG::dow; ++l) {
                                    elementVector[local_kj] += spontVal * Dune::at(normal,k) * Dune::at(normal,l) * gradients[j][l];

                                }
                                if (!incompressible_)
                                    elementVector[local_kj] += kappa0Factor * Dune::at(normal,k) * shapeValues[j];
                            }
                        }
                        if (axi_) {
                            const auto local_j = fShellNode.child(1).localIndex(j);
                            if (lineTension_) {
                                elementVector[local_j] += -lineTensionFactor * shapeValues[j];
                            }
                        }
                    }

                    // contact angle condition due to phase field in contact with the surface
                    // assemble the following operators from outside, if outside_ == 1, and from inside if inside_ == 1
                    auto phiVal = phi(local);
                    auto contactFactor = 6.0 * phiVal * (1.0 - phiVal) * (Dune::at(global,0) > center ? sigma1_ - sigma0_ : sigma1_nucleus_ - sigma0_nucleus_) * dx;
                    for (std::size_t i = 0; i < muSize; ++i) {
                        const auto local_i = muNode.localIndex(i);
                        const auto value = contactFactor * muValues[i];

                        elementVector[local_i] += value;
                    }
                }
                phi.unbind();
                gradPhi.unbind();
                normalX.unbind();
                normalY.unbind();
            }

        }
    protected:

        template<class S, class F, class T, int dow,
                std::enable_if_t<Category::Scalar < S>, int> = 0>

        T eval(S const &scalar, F factor,
               Dune::FieldVector<T, dow> const &grad_test,
               Dune::FieldVector<T, dow> const &grad_trial) const {
            return (scalar * factor) * (grad_test * grad_trial);
        }

        template<class M, class F, class T, int dow,
                std::enable_if_t<Category::Matrix < M>, int> = 0>

        T eval(M const &mat, F factor,
               Dune::FieldVector<T, dow> const &grad_test,
               Dune::FieldVector<T, dow> const &grad_trial) const {
                 return factor * (grad_test * (mat * grad_trial));
        }
    };

    template<class LC, class DV, class DV2, class DV3, class IS>
    struct GridFunctionOperatorRegistry<tag::surface_fShell<DV,DV2,DV3,IS>, LC> {
        static constexpr int degree = 2;
        using type = SurfaceFShell<DV,DV2,DV3,IS>;
    };
    /** @} **/

} // end namespace AMDiS
