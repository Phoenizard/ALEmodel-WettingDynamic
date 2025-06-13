#pragma once

#include <type_traits>

#include <amdis/GridFunctionOperator.hpp>
#include <amdis/common/StaticSize.hpp>
#include <amdis/ALEModelAmdis2/myIndices.hpp>

using namespace MyIndices;
namespace AMDiS {
    /**
     * \addtogroup operators
     * @{
     **/

    namespace tag {
        template<class DF, class IS = int>
        struct surface_cahn_hilliard {
            unsigned int axi = 0; // axisyppetric or not?
            double const& sigma = 0; // the surface phase field tension
            double const& eps = 0; // the surface phase field interface width
            double const& mobility = 0; // the constant surface mobility
            DF const& phi = {}; // the phase field
            double const& invTau = 0; // the inverse of the time step size
            double const& k0 = 0; // prefactor \frac{k_0\nu}{k_bT} of the source term
            double const& kOn = 0; // binding factor for binding of molecules from bulk onto surface
            double const& kOff = 0; // binding factor for binding of molecules from surface into buld
            double const& molBulkVol = 0; // molecular bulk volume on the membrane
            double const& molSurfArea = 0; // molecular surface area on the membrane
            std::size_t inside = 0; // integrate over the shell from inside (0) or outside (1)
            IS const &indexSet = {}; //the indexSet of the leafGridView of the grid!
            std::vector<int> const &partitions = {}; //partition numbers of each element
            std::vector<int> const &facets = {}; //facet number of the surface facet (-1 if element i not touching the surface)
        };
    }


    /// Surface Cahn Hilliard equation with constant mobility and single well potential W = 1/4(\phi_\Gamma - 1/2)^4
    /// Including source terms to couple the surface equations to a bulk NSCH equation
    /// The source term reads: s = \frac{k_0\nu}{k_bT}\left(-\frac{\overline{\nu}}{\nu}\mu_\Gamma + \mu\right)
    ///                            + k_{on}\phi - k_{off}\phi_\Gamma
    /// This term is also added as an internal BC to the phase field equation in the bulk.
    template<class DF, class IS = int>
    class SurfaceCahnHilliard {
        unsigned int axi_ = 0;
        double sigma_;
        double eps_;
        double const& mobility_;
        DF const& phi_ = {}; // the phase field
        double const& invTau_;
        double k0_ = 0; //
        double kOn_ = 0;
        double kOff_ = 0;
        double molBulkVol_ = 0;
        double molSurfArea_ = 0;
        IS const &indexSet_ = {};
        std::vector<int> const &partitions_;
        std::vector<int> const &facets_;
        std::size_t inside_;

    public:
        SurfaceCahnHilliard(tag::surface_cahn_hilliard<DF,IS> t)
                : axi_(t.axi), sigma_(t.sigma), eps_(t.eps), mobility_(t.mobility),
                  invTau_(t.invTau), k0_(t.k0), kOn_(t.kOn), kOff_(t.kOff),
                  molBulkVol_(t.molBulkVol), molSurfArea_(t.molSurfArea),
                  partitions_(t.partitions), facets_(t.facets),
                  indexSet_(t.indexSet), inside_(t.inside), phi_(FWD(t.phi)) {}

        template<class CG, class Node, class Quad, class LocalFct, class Mat>
        void assemble(CG const &contextGeo, Node const &tree, Node const /*colNode*/,
                      Quad const &quad, LocalFct const &localFct, Mat &elementMatrix) const {
            static_assert(static_size_v<typename LocalFct::Range> == 1,
                          "Expression must be of scalar type.");

            using VelocityNode = Dune::TypeTree::Child<Node, _vS.value>;
            using PhiBulkNode = Dune::TypeTree::Child<Node, _phi.value>;
            using MuBulkNode = Dune::TypeTree::Child<Node, _mu.value>;
            using PhiNode = Dune::TypeTree::Child<Node, _phiGamma.value>;
            using MuNode = Dune::TypeTree::Child<Node, _muGamma.value>;
            
            static_assert(PhiNode::isLeaf && MuNode::isLeaf && 
                          PhiBulkNode::isLeaf && MuBulkNode::isLeaf &&
                          VelocityNode::isPower,
                          "Operator can be applied to correct nodes only.");

            auto newS = Parameters::get<int>("use new binding flux").value_or(0);

            auto const &vNode = tree.child(_vS);
            auto const &phiNode = tree.child(_phiGamma);
            auto const &muNode = tree.child(_muGamma);
            auto const &phiBulkNode = tree.child(_phi);
            auto const &muBulkNode = tree.child(_mu);
            auto const &fShellNode = tree.child(_fShell);

            std::size_t vSize = vNode.child(0).size();
            std::size_t phiSize = phiNode.size();
            std::size_t muSize = muNode.size();
            std::size_t phiBulkSize = phiBulkNode.size();
            std::size_t muBulkSize = muBulkNode.size();
            std::size_t fShellSize = fShellNode.child(0).size();

            if (phiNode.size() == 0 || muNode.size() == 0)
                return;

            // compute facet number and partition number from surface basis node (if no surface basis node: return -2)
            auto part = std::max(partition(phiNode, Dune::PriorityTag<2>{}),
                                 partition(muNode, Dune::PriorityTag<2>{}));
            auto face = std::max(facet(phiNode, Dune::PriorityTag<2>{}), facet(muNode, Dune::PriorityTag<2>{}));

            if (face == -2) {
                auto index = idx(indexSet_, contextGeo.context(), Dune::PriorityTag<2>{});
                test_exit(!(partitions_.size() == 0 || facets_.size() == 0 || index == -1),
                          "please provide partitions, facets and indexSet of the leafGridView when test and trial "
                          "functions are not from the surfaceLagrange basis!");
                face = facets_[index];
                part = partitions_[index];
            }

            if (face == -1 || part != inside_) // do nothing if element has no facet on the shell
                return;

            using RangeFieldType = typename PhiNode::LocalBasis::Traits::RangeFieldType;
            using WorldVector = Dune::FieldVector<RangeFieldType, CG::dow>;
            using WorldMatrix = Dune::FieldMatrix<RangeFieldType,CG::dow,CG::dow>;
            std::vector<WorldVector> phiGradients, phiGradientsAux;

            using MuFieldType = typename MuNode::LocalBasis::Traits::RangeFieldType;
            using MuWorldVector = Dune::FieldVector<MuFieldType, CG::dow>;
            std::vector<MuWorldVector> muGradients, muGradientsAux;

            using VFieldType = typename VelocityNode::ChildType::LocalBasis::Traits::RangeFieldType;
            using VWorldVector = Dune::FieldVector<VFieldType, CG::dow>;
            std::vector<VWorldVector> vGradients, vGradientsAux, fShellGradients, fShellGradientsAux;

            auto phiDim = basisDim(phiNode, Dune::PriorityTag<2>{});
            auto muDim = basisDim(muNode, Dune::PriorityTag<2>{});
            auto vDim = basisDim(vNode.child(0), Dune::PriorityTag<2>{});

            auto gfPhi = makeGridFunction(valueOf(*phi_), phi_->basis().gridView());
            auto phiBulk = localFunction(gfPhi);
            phiBulk.bind(contextGeo.element());

            auto element = contextGeo.context();
            auto intersection = element.template subEntity<1>(face);

            constexpr auto intersectionDim = std::decay_t<decltype(intersection)>::mydimension;
            const auto &quadShell = Dune::QuadratureRules<double, intersectionDim>::rule(intersection.type(), 2);

            auto gradPhiGammaLF = derivativeOf(localFct,tag::gradient{});
            gradPhiGammaLF.bind(contextGeo.element());

            // apply the operator only to intersections on the shell
            if (part == inside_) {
                auto geo = intersection.geometry();

                for (auto const &qp: quadShell) {
                    auto localGeo = referenceElement(element).template geometry<1>(face);

                    // Position of the current quadrature point in the reference element
                    auto &&local = localGeo.global(qp.position());
                    auto &&global = contextGeo.elementGeometry().global(local);

                    // The transposed inverse Jacobian of the map from the reference facet to the real facet
                    const auto jacobian = geo.jacobianInverseTransposed(qp.position());

                    // The transposed inverse Jacobian of the map from the reference facet to the reference element
                    const auto jacobianRef = localGeo.jacobianTransposed(qp.position());

                    // The multiplicative factor in the integral transformation formula
                    const auto factor = geo.integrationElement(qp.position()) * qp.weight();
                    const auto phi = localFct(local);
                    auto phiB = phiBulk(local);
                    auto axiFactor = (axi_ * Dune::at(global,1) + 1.0 - axi_);
                    auto gradPhiGamma = gradPhiGammaLF(local);
                    auto gradOld = gradPhiGamma;
                    //phiB = (double)phiB > 1.0 ? 1.0 : ((double)phiB < 0.0 ? 0.0 : (double)phiB); // clamp phiB between 0 and 1

                    // The gradients of the shape functions on the reference element
                    auto const &phiShapeGradients = phiNode.localBasisJacobiansAt(local);
                    auto const &muShapeGradients = muNode.localBasisJacobiansAt(local);
                    auto const &vShapeGradients = vNode.child(0).localBasisJacobiansAt(local);
                    auto const &fShellShapeGradients = fShellNode.child(0).localBasisJacobiansAt(local);

                    phiGradients.resize(phiShapeGradients.size());
                    muGradients.resize(muShapeGradients.size());
                    phiGradientsAux.resize(phiShapeGradients.size());
                    muGradientsAux.resize(muShapeGradients.size());
                    vGradients.resize(vShapeGradients.size());
                    vGradientsAux.resize(vShapeGradients.size());
                    fShellGradients.resize(fShellShapeGradients.size());
                    fShellGradientsAux.resize(fShellShapeGradients.size());

                    // transform local gradients to reference intersection if necessary and the compute
                    // gradients on the real intersection
                    int dow = CG::dow;
                    if (phiDim == dow) {
                        for (std::size_t i = 0; i < phiGradients.size(); ++i)
                            jacobianRef.mv(phiShapeGradients[i][0], phiGradientsAux[i]);

                        for (std::size_t i = 0; i < phiGradients.size(); ++i)
                            jacobian.mv(phiGradientsAux[i], phiGradients[i]);
                    } else {
                        for (std::size_t i = 0; i < phiGradients.size(); ++i)
                            jacobian.mv(phiShapeGradients[i][0], phiGradients[i]);
                    }
                    if (muDim == dow) {
                        for (std::size_t i = 0; i < muGradients.size(); ++i)
                            jacobianRef.mv(muShapeGradients[i][0], muGradientsAux[i]);

                        for (std::size_t i = 0; i < muGradients.size(); ++i)
                            jacobian.mv(muGradientsAux[i], muGradients[i]);
                    } else {
                        for (std::size_t i = 0; i < muGradients.size(); ++i)
                            jacobian.mv(muShapeGradients[i][0], muGradients[i]);
                    }
                    if (vDim == dow) {
                        for (std::size_t i = 0; i < vGradients.size(); ++i)
                            jacobianRef.mv(vShapeGradients[i][0], vGradientsAux[i]);

                        for (std::size_t i = 0; i < vGradients.size(); ++i)
                            jacobian.mv(vGradientsAux[i], vGradients[i]);
                    } else {
                        for (std::size_t i = 0; i < vGradients.size(); ++i)
                            jacobian.mv(vShapeGradients[i][0], vGradients[i]);
                    }

                    for (std::size_t i = 0; i < fShellGradients.size(); ++i)
                        jacobian.mv(fShellShapeGradients[i][0], fShellGradients[i]);


                    auto factorDoubleWell = sigma_ / eps_ * (-3.0 * std::pow(phi - 0.5, 2));

                    auto const &phiValues = phiNode.localBasisValuesAt(local);
                    auto const &muValues = muNode.localBasisValuesAt(local);
                    auto const &phiBValues = phiBulkNode.localBasisValuesAt(local);
                    auto const &muBValues = muBulkNode.localBasisValuesAt(local);
                    auto const &vValues = vNode.child(0).localBasisValuesAt(local);
                    auto const &fShellValues = fShellNode.child(0).localBasisValuesAt(local);
                    for (std::size_t i = 0; i < phiSize; ++i) {
                        const auto local_i = phiNode.localIndex(i);
                        const auto local_ii = muNode.localIndex(i);
                        const auto valuePhi = factor * phiValues[i];
                        const auto valueMu = factor * muValues[i];
                        const auto valueDoubleWell = factorDoubleWell * factor * phiValues[i];
                        const auto valueV = phi * factor * phiValues[i];
                        const auto valueMuSource = k0_ * molSurfArea_ / molBulkVol_ * factor * phiValues[i];

                        elementMatrix[local_i][local_i] +=
                                (invTau_ + kOff_ * (newS ? std::pow(phiB[0],2) : 1.0)) * valuePhi * phiValues[i]; // lhs of d_t phi
                        elementMatrix[local_ii][local_ii] += valueMu * muValues[i]; // mu = ...
                        if (newS) elementMatrix[local_i][local_i] += kOn_ * valuePhi * phiValues[i]; // source term

                        for (std::size_t j = i + 1; j < phiSize; ++j) {
                            const auto local_j = phiNode.localIndex(j);
                            const auto local_jj = muNode.localIndex(j);

                            elementMatrix[local_i][local_j] +=
                                    (invTau_ + kOff_ * (newS ? std::pow(phiB[0],2) : 1.0)) * valuePhi * phiValues[j]; // lhs of d_t phi
                            elementMatrix[local_j][local_i] +=
                                    (invTau_ + kOff_ * (newS ? std::pow(phiB[0],2) : 1.0)) * valuePhi * phiValues[j]; // lhs of d_t phi

                            if (newS) {
                                elementMatrix[local_i][local_j] += kOn_ * valuePhi * phiValues[j]; // source term
                                elementMatrix[local_j][local_i] += kOn_ * valuePhi * phiValues[j]; // source term
                            }
                            elementMatrix[local_ii][local_jj] += valueMu * muValues[j]; // mu = ...
                            elementMatrix[local_jj][local_ii] += valueMu * muValues[j]; // mu = ...
                        }

                        for (std::size_t j = 0; j < muSize; ++j) {
                            const auto local_j = muNode.localIndex(j);
                            elementMatrix[local_i][local_j] += eval(mobility_, factor, phiGradients[i],
                                                                    muGradients[j]); // sot(mobility) (_phi, _mu)
                            elementMatrix[local_j][local_i] += eval(-sigma_ * eps_, factor, muGradients[j],
                                                                    phiGradients[i]); // sot(-sigma*eps) (_mu, _phi)

                            elementMatrix[local_j][local_i] +=
                                    valueDoubleWell * muValues[j]; // W'(phi) implicit part (linearized)

                            elementMatrix[local_i][local_j] += valueMuSource * muValues[j]; // source term

                            // axisymmetric terms
                            if (axi_) {
                                const auto valuePhi =
                                        -mobility_ * factor / Dune::at(global,1) * muGradients[j][1];
                                elementMatrix[local_i][local_j] += valuePhi * phiValues[i];

                                const auto valueMu = eps_ * sigma_ * factor / Dune::at(global,1) * phiGradients[i][1];
                                elementMatrix[local_j][local_i] += valueMu * muValues[j];
                            }
                        }

                        // coupling terms to velocity
                        auto ignoreLaplace = Parameters::get<int>("ignore laplace problem").value_or(0);
                        if (!ignoreLaplace) {
                            for (std::size_t j = 0; j < vSize; ++j) {
                                for (std::size_t k = 0; k < CG::dow; ++k) {
                                    const auto local_vk = vNode.child(k).localIndex(j);
                                    elementMatrix[local_i][local_vk] +=
                                            valueV * vGradients[j][k]; // <phi^old div(u), v>

                                    if (axi_ && k == 1)
                                        elementMatrix[local_i][local_vk] +=
                                                valueV / Dune::at(global,1) * vValues[j]; // axi coupling term
                                }
                            }
                        }

                        // source terms with phiBulk and muBulk
                        for (std::size_t j = 0; j < phiBulkSize; ++j) { //it is phiBulkSize == muBulkSize
                            const auto local_jB = phiBulkNode.localIndex(j);
                            const auto local_jjB = muBulkNode.localIndex(j);

                            if (!newS) elementMatrix[local_i][local_jB] += -kOn_ * valuePhi * phiBValues[j]; // source term
                            elementMatrix[local_i][local_jjB] += -k0_ * valuePhi * muBValues[j]; // source term
                        }
                    }

                    for (std::size_t i = 0; i < phiBulkSize; ++i) {
                        const auto local_i = phiBulkNode.localIndex(i);
                        const auto local_ii = muBulkNode.localIndex(i);

                        const auto value = kOn_  * molBulkVol_ / molSurfArea_ * factor * phiBValues[i];
                        const auto valueK0 = k0_ * molBulkVol_ / molSurfArea_ * factor * phiBValues[i];
                        const auto valueMinusK0 = -k0_ * factor * phiBValues[i];
                        const auto valueKOff = -kOff_ * molBulkVol_ / molSurfArea_ * (newS ? std::pow(phiB[0],2) : 1.0) * factor * phiBValues[i];
                        if (!newS) elementMatrix[local_i][local_i] += value * phiBValues[i];

                        for (std::size_t j = i + 1; j < phiBulkSize; ++j) {
                            const auto local_j = phiBulkNode.localIndex(j);
                            const auto local_jj = muBulkNode.localIndex(j);

                            if (!newS) {
                                elementMatrix[local_i][local_j] += value * phiBValues[j];
                                elementMatrix[local_j][local_i] += value * phiBValues[j];
                            }
                        }
                        for (std::size_t j = 0; j < muBulkSize; ++j) {
                            const auto local_j = muBulkNode.localIndex(j);
                            elementMatrix[local_i][local_j] += valueK0 * muBValues[j];

                        }
                        if (muSize == phiSize) {
                            for (std::size_t j = 0; j < phiSize; ++j) {
                                const auto local_j = phiNode.localIndex(j);
                                const auto local_jj = muNode.localIndex(j);
                                const auto value = -kOn_ * newS * molBulkVol_ / molSurfArea_ * factor * phiBValues[i];

                                elementMatrix[local_i][local_j] += (valueKOff + value) * phiValues[j];
                                elementMatrix[local_i][local_jj] += valueMinusK0 * muValues[j];
                            }
                        } else {
                            for (std::size_t j = 0; j < phiSize; ++j) {
                                const auto local_j = phiNode.localIndex(j);
                                const auto value = -kOn_ * newS * molBulkVol_ / molSurfArea_ * factor * phiBValues[i];

                                elementMatrix[local_i][local_j] += (valueKOff + value) * phiValues[j];
                            }
                            for (std::size_t j = 0; j < muSize; ++j) {
                                const auto local_jj = muNode.localIndex(j);
                                elementMatrix[local_i][local_jj] += valueMinusK0 * muValues[j];
                            }
                        }
                    }

                    // contribution to shell force
                    for (int ii = 0; ii < dow; ++ii) {
                        for (std::size_t j = 0; j < phiSize; ++j) {
                            const auto local_j = phiNode.localIndex(j);
                            const auto value =
                                    (eps_ * sigma_) * 0.5 * factor * axiFactor * gradOld[ii] * phiGradients[j][ii];
                            for (std::size_t i = 0; i < fShellSize; ++i) {
                                for (int kk = 0; kk < dow; ++kk) {
                                    const auto value2 =
                                            -(eps_ * sigma_) * factor * axiFactor * gradOld[ii] * phiGradients[j][kk];
                                    const auto local_i = fShellNode.child(kk).localIndex(i);
                                    elementMatrix[local_i][local_j] += value * fShellGradients[i][kk];
                                    elementMatrix[local_i][local_j] += value2 * fShellGradients[i][ii];
                                }
                                if (axi_) {
                                    const auto local_i = fShellNode.child(1).localIndex(i);
                                    elementMatrix[local_i][local_j] += value / Dune::at(global,1) * fShellValues[i];
                                }
                            }
                        }
                    }
                }
            }
            gradPhiGammaLF.unbind();
        }

        template<class CG, class Node, class Quad, class LF, class Vec>
        void assemble(CG const &contextGeo, Node const &tree, Quad const &quad,
                      LF const &localFct, Vec &elementVector) const {
            static_assert(static_size_v<typename LF::Range> == 1,
                          "Expression must be of scalar type.");

            using PhiNode = Dune::TypeTree::Child<Node,_phiGamma.value>;
            using PhiBulkNode = Dune::TypeTree::Child<Node,_phi.value>;
            using MuNode = Dune::TypeTree::Child<Node,_muGamma.value>;
            using FShellNode = Dune::TypeTree::Child<Node, _fShell.value>;

            static_assert(PhiNode::isLeaf && MuNode::isLeaf, "");

            auto newS = Parameters::get<int>("use new binding flux").value_or(0);

            auto const& phiNode = tree.child(_phiGamma);
            auto const& muNode = tree.child(_muGamma);
            auto const &phiBulkNode = tree.child(_phi);
            auto const &fShellNode = tree.child(_fShell);

            std::size_t phiSize = phiNode.size();
            std::size_t phiBulkSize = phiBulkNode.size();
            std::size_t muSize = muNode.size();
            std::size_t fShellSize = fShellNode.child(0).size();

            if (phiNode.size() == 0)
                return;

            // compute facet number and partition number from surface basis node (if no surface basis node: return -2)
            auto part = partition(phiNode, Dune::PriorityTag<2>{});
            auto face = facet(phiNode, Dune::PriorityTag<2>{});

            if (face == -2) {
                auto index = idx(indexSet_, contextGeo.context(), Dune::PriorityTag<2>{});
                test_exit(!(partitions_.size() == 0 || facets_.size() == 0 || index == -1),
                          "please provide partitions, facets and indexSet of the leafGridView when test and trial "
                          "functions are not from the surfaceLagrange basis!");
                face = facets_[index];
                part = partitions_[index];
            }

            if (face == -1 || part != inside_) // do nothing if element has no facet on the shell
                return;

            auto element = contextGeo.context();
            auto intersection = element.template subEntity<1>(face);

            constexpr auto intersectionDim = std::decay_t<decltype(intersection)>::mydimension;
            const auto &quadShell = Dune::QuadratureRules<double, intersectionDim>::rule(intersection.type(), 2);

            using RangeFieldType = typename FShellNode::ChildType::LocalBasis::Traits::RangeFieldType;
            using WorldVector = FieldVector<RangeFieldType, CG::dow>;
            std::vector<WorldVector> fShellGradients;

            // apply the operator only to intersections on the shell
            if (part == inside_) {
                auto geo = intersection.geometry();

                for (auto const &qp: quadShell) {
                    auto localGeo = referenceElement(element).template geometry<1>(face);

                    // The transposed inverse Jacobian of the map from the reference facet to the real facet
                    const auto jacobian = geo.jacobianInverseTransposed(qp.position());

                    // Position of the current quadrature point in the reference element
                    auto &&local = localGeo.global(qp.position());
                    auto &&global = contextGeo.elementGeometry().global(local);

                    // The multiplicative factor in the integral transformation formula
                    const auto factor = geo.integrationElement(qp.position()) * qp.weight();

                    auto phi = localFct(local);
                    auto phiFactor = (invTau_ * phi + kOn_ * newS) * factor;
                    auto muFactor = sigma_ / eps_ * std::pow(phi - 0.5, 2) * (-2.0 * phi - 0.5) * factor;
                    const auto dx = geo.integrationElement(qp.position()) * qp.weight();
                    const auto dxAxi = (axi_ * Dune::at(global,1) + 1.0 - axi_) * dx;
                    const auto fShellFactor = -sigma_ / eps_ * 0.25 * std::pow(phi - 0.5,4) * dxAxi;

                    // The gradients of the shape functions on the reference element
                    auto const &shapeGradients = fShellNode.child(0).localBasisJacobiansAt(local);

                    // Compute the shape function gradients on the real element
                    fShellGradients.resize(shapeGradients.size());

                    for (std::size_t i = 0; i < fShellGradients.size(); ++i)
                        jacobian.mv(shapeGradients[i][0], fShellGradients[i]);

                    auto const &phiValues = phiNode.localBasisValuesAt(local);
                    auto const &phiBValues = phiBulkNode.localBasisValuesAt(local);
                    auto const &muValues = muNode.localBasisValuesAt(local);
                    if (phiSize == muSize) {
                        for (std::size_t i = 0; i < phiSize; ++i) {
                            const auto local_i = phiNode.localIndex(i);
                            const auto local_ii = muNode.localIndex(i);
                            elementVector[local_i] += phiFactor * phiValues[i];
                            elementVector[local_ii] += muFactor * muValues[i];
                        }
                    } else {
                        for (std::size_t i = 0; i < phiSize; ++i) {
                            const auto local_i = phiNode.localIndex(i);
                            elementVector[local_i] += phiFactor * phiValues[i];
                        }
                        for (std::size_t i = 0; i < muSize; ++i) {
                            const auto local_ii = muNode.localIndex(i);
                            elementVector[local_ii] += muFactor * muValues[i];
                        }
                    }
                    if (newS) {
                        for (std::size_t i = 0; i < phiBulkSize; ++i) {
                            const auto local_i = phiBulkNode.localIndex(i);
                            elementVector[local_i] += -kOn_ * molBulkVol_/molSurfArea_ * factor * phiBValues[i];
                        }
                    }

                    // shell force contribution
                    auto const &shapeValues = fShellNode.child(1).localBasisValuesAt(local);

                    for (std::size_t j = 0; j < fShellSize; ++j) {
                        for (std::size_t k = 0; k < CG::dow; ++k) {
                            const auto local_kj = fShellNode.child(k).localIndex(j);
                            elementVector[local_kj] += fShellFactor * fShellGradients[j][k];
                        }
                        if (axi_) {
                            const auto local_j = fShellNode.child(1).localIndex(j);
                            elementVector[local_j] += fShellFactor / geo.corner(j)[1] * shapeValues[j];
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

    template<class LC, class DF, class IS>
    struct GridFunctionOperatorRegistry<tag::surface_cahn_hilliard<DF, IS>, LC> {
        static constexpr int degree = 2;
        using type = SurfaceCahnHilliard<DF,IS>;
    };


    /** @} **/

} // end namespace AMDiS
