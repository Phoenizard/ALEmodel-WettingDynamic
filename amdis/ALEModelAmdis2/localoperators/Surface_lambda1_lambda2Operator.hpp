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
        template<class DV, class IS = int>
        struct surface_lambda1_lambda2 {
            unsigned int axi = 0; // axisymmetric or not?
            double const& tau; // the time step size
            double const& Ks = 0; // area shear modulus on the surface
            DV const& lambda2; // second principal stretch in axisymmetric case (\lambda2 = R/R_0, with R the distance to the symmetry axis)
            std::size_t inside = 0; // integrate over the shell from inside (0) or outside (1)
            IS const &indexSet = {}; //the indexSet of the leafGridView of the grid!
            std::vector<int> const &partitions = {}; //partition numbers of each element
            std::vector<int> const &facets = {}; //facet number of the surface facet (-1 if element i not touching the surface)
        };        
    }


    /// An operator for the computation of the two principal stretches of a surface (in axisymmetric and general case)
    /// needed for the stretching force of a shell in ALE simulations.
    template<class DV, class IS = int>
    class Surface_lambda1_lambda2 {
        unsigned int axi_;
        double const &tau_;
        double const &Ks_;
        DV const &lambda2_;
        IS const &indexSet_ = {};
        std::vector<int> const &partitions_;
        std::vector<int> const &facets_;
        std::size_t inside_;
    public:
        Surface_lambda1_lambda2(tag::surface_lambda1_lambda2<DV, IS> t)
                : axi_(t.axi), tau_(t.tau), Ks_(t.Ks), lambda2_(FWD(t.lambda2)), indexSet_(t.indexSet),
                  partitions_(t.partitions), facets_(t.facets), inside_(t.inside) {}

        // assemble method for LHS
        template<class CG, class Node, class Quad, class LocalFct, class Mat>
        void assemble(CG const &contextGeo, Node const &tree, Node const /*colNode*/,
                      Quad const &quad, LocalFct const &localFct, Mat &elementMatrix) const {
            static_assert(static_size_v<typename LocalFct::Range> == 1,
                          "Expression must be of scalar type.");

            using VelocityNode = Dune::TypeTree::Child<Node, _vS.value>;
            using Lambda1Node = Dune::TypeTree::Child<Node, _lambda1.value>;
            using Lambda2Node = Dune::TypeTree::Child<Node, _lambda2.value>;


            static_assert(VelocityNode::isPower && Lambda1Node::isLeaf && Lambda2Node::isLeaf,
                          "Nodes must have correct dimensions.");

            auto const &vSNode = tree.child(_vS);
            auto const &lambda1Node = tree.child(_lambda1);
            auto const &lambda2Node = tree.child(_lambda2);

            if (lambda1Node.size() == 0 || lambda2Node.size() == 0 || vSNode.child(0).size() == 0)
                return;

            // compute facet number and partition number from surface basis node (if no surface basis node: return -2)
            auto part = std::max(partition(lambda1Node, Dune::PriorityTag<2>{}),
                                 partition(vSNode.child(0), Dune::PriorityTag<2>{}));
            auto face = std::max(facet(lambda1Node, Dune::PriorityTag<2>{}),
                                 facet(vSNode.child(0), Dune::PriorityTag<2>{}));

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

            std::size_t lambda1Size = lambda1Node.size();
            std::size_t lambda2Size = lambda2Node.size();
            std::size_t vSSize = vSNode.child(0).size();

            using FieldType = typename VelocityNode::ChildType::LocalBasis::Traits::RangeFieldType;
            using WorldVector = Dune::FieldVector<FieldType, CG::dow>;
            std::vector<WorldVector> gradientsVS, gradientsVSAux;

            auto vSDim = basisDim(vSNode.child(0), Dune::PriorityTag<2>{});

            auto element = contextGeo.element();

            auto intersection = element.template subEntity<1>(face);

            constexpr auto intersectionDim = std::decay_t<decltype(intersection)>::mydimension;
            const auto &quadShell = Dune::QuadratureRules<double, intersectionDim>::rule(intersection.type(), 2);

            auto gfLambda2 = makeGridFunction(valueOf(*lambda2_), lambda2_->basis().gridView());
            auto lambda2 = localFunction(gfLambda2);
            lambda2.bind(contextGeo.element());

            // apply the operator only to intersections on the shell
            auto geo = intersection.geometry();

            if (part == inside_) {
                for (auto const &qp: quadShell) {
                    auto localGeo = referenceElement(element).template geometry<1>(face);

                    // Position of the current quadrature point in the reference element
                    auto &&local = localGeo.global(qp.position());
                    auto &&global = contextGeo.elementGeometry().global(local);

                    // The multiplicative factor in the integral transformation formula
                    const auto factor = geo.integrationElement(qp.position()) * qp.weight();
                    const auto distances = localFct(local);

                    // The transposed inverse Jacobian of the map from the reference facet to the real facet
                    const auto jacobian = geo.jacobianInverseTransposed(qp.position());

                    // The transposed inverse Jacobian of the map from the reference facet to the reference element
                    const auto jacobianRef = localGeo.jacobianTransposed(qp.position());

                    // The gradients of the shape functions on the reference element
                    auto const &shapeGradientsVS = vSNode.child(0).localBasisJacobiansAt(local);

                    gradientsVS.resize(shapeGradientsVS.size());
                    gradientsVSAux.resize(shapeGradientsVS.size());

                    // transform local gradients to reference intersection if necessary and the compute
                    // gradients on the real intersection
                    int dow = CG::dow;
                    if (vSDim == dow) {
                        for (std::size_t i = 0; i < gradientsVS.size(); ++i)
                            jacobianRef.mv(shapeGradientsVS[i][0], gradientsVSAux[i]);

                        for (std::size_t i = 0; i < gradientsVS.size(); ++i)
                            jacobian.mv(gradientsVSAux[i], gradientsVS[i]);
                    } else {
                        for (std::size_t i = 0; i < gradientsVS.size(); ++i)
                            jacobian.mv(shapeGradientsVS[i][0], gradientsVS[i]);
                    }

                    auto const &lambda1Values = lambda1Node.localBasisValuesAt(local);
                    auto const &lambda2Values = lambda2Node.localBasisValuesAt(local);
                    auto const &vSValues = vSNode.child(0).localBasisValuesAt(local);

                    for (std::size_t i = 0; i < lambda1Size; ++i) {
                        const auto value = factor * lambda1Values[i];
                        const auto value2 = factor * lambda2Values[i];

                        // lambda1 = ..., lambda2 = ...
                        const auto local_i = lambda1Node.localIndex(i);
                        const auto local_ii = lambda2Node.localIndex(i);
                        elementMatrix[local_i][local_i] += value * lambda1Values[i];
                        elementMatrix[local_ii][local_ii] += value2 * lambda2Values[i];

                        for (std::size_t j = i + 1; j < lambda1Size; ++j) {
                            const auto local_j = lambda1Node.localIndex(j);
                            const auto local_jj = lambda2Node.localIndex(j);

                            elementMatrix[local_i][local_j] += value * lambda1Values[j];
                            elementMatrix[local_j][local_i] += value * lambda1Values[j];
                            elementMatrix[local_ii][local_jj] += value2 * lambda2Values[j];
                            elementMatrix[local_jj][local_ii] += value2 * lambda2Values[j];
                        }

                        // = tau * lambda1 * div v + axisymmetric terms
                        const auto valueV = -tau_ * distances * (factor * lambda1Values[i]);
                        auto valueV2 = valueV;
                        if (axi_ && Ks_) valueV2 = -tau_ * lambda2(local) / Dune::at(global,1) * (factor * lambda2Values[i]);
                        for (std::size_t j = 0; j < vSSize; ++j) {
                            for (std::size_t k = 0; k < CG::dow; ++k) {
                                const auto local_kj = vSNode.child(k).localIndex(j);
                                elementMatrix[local_i][local_kj] += valueV * gradientsVS[j][k];
                                if (axi_ && !Ks_ && k == 1)
                                    elementMatrix[local_i][local_kj] += valueV / Dune::at(global,1) * vSValues[j];
                                if (axi_ && Ks_ && k == 1)
                                    elementMatrix[local_ii][local_kj] += valueV2 * vSValues[j];
                            }
                        }
                    }
                }
            }
            lambda2.unbind();
        }

        // assemble method for RHS
        template<class CG, class Node, class Quad, class LocalFct, class Vec>
        void assemble(CG const &contextGeo, Node const &tree, Quad const &quad,
                      LocalFct const &localFct, Vec &elementVector) const {
            static_assert(static_size_v<typename LocalFct::Range> == 1,
                          "Expression must be of scalar type.");

            using Lambda1Node = Dune::TypeTree::Child<Node, _lambda1.value>;
            using Lambda2Node = Dune::TypeTree::Child<Node, _lambda2.value>;

            static_assert(Lambda1Node::isLeaf && Lambda2Node::isLeaf,
                          "Nodes must have correct dimensions.");

            auto const &lambda1Node = tree.child(_lambda1);
            auto const &lambda2Node = tree.child(_lambda2);

            if (lambda1Node.size() == 0 || lambda2Node.size() == 0)
                return;

            auto lambda1Size = lambda1Node.size();
            auto lambda2Size = lambda2Node.size();

            // compute facet number and partition number from surface basis node (if no surface basis node: return -2)
            auto part = partition(lambda1Node, Dune::PriorityTag<2>{});
            auto face = facet(lambda1Node, Dune::PriorityTag<2>{});

            if (face == -2) {
                auto index = idx(indexSet_, contextGeo.context(), Dune::PriorityTag<2>{});
                test_exit(!(partitions_.size() == 0 || facets_.size() == 0 || index == -1),
                          "please provide partitions, facets and indexSet of the leafGridView when test and trial "
                          "functions are not from the surfaceLagrange basis!");
                face = facets_[index];
                part = partitions_[index];
            }

            if (face == -1 || part != inside_) // do nothing if the element is on the wrong partition
                return;

            auto element = contextGeo.context();

            auto intersection = element.template subEntity<1>(face);

            constexpr auto intersectionDim = std::decay_t<decltype(intersection)>::mydimension;
            const auto &quadShell = Dune::QuadratureRules<double, intersectionDim>::rule(intersection.type(), 2);

            auto gfLambda2 = makeGridFunction(valueOf(*lambda2_), lambda2_->basis().gridView());
            auto lambda2 = localFunction(gfLambda2);
            lambda2.bind(contextGeo.element());

            // apply the operator only to intersections on the shell
            if (part == inside_) {
                auto geo = intersection.geometry();

                for (auto const &qp: quadShell) {
                    auto localGeo = referenceElement(element).template geometry<1>(face);

                    // Position of the current quadrature point in the reference element
                    auto &&local = localGeo.global(qp.position());
                    auto &&global = contextGeo.elementGeometry().global(local);

                    // The multiplicative factor in the integral transformation formula
                    const auto factor = geo.integrationElement(qp.position()) * qp.weight();
                    const auto distances = localFct(local);
                    const auto lambda1Factor = (distances - 1.0) * factor;
                    auto lambda2Factor = lambda1Factor;
                    if (axi_ && Ks_) lambda2Factor = (lambda2(local) - 1.0) * factor;
                    auto const &lambda1Values = lambda1Node.localBasisValuesAt(local);
                    auto const &lambda2Values = lambda2Node.localBasisValuesAt(local);
                    for (std::size_t i = 0; i < lambda1Size; ++i) {
                        const auto local_i = lambda1Node.localIndex(i);
                        elementVector[local_i] += lambda1Factor * lambda1Values[i];
                        if (axi_ && Ks_) {
                            const auto local_ii = lambda2Node.localIndex(i);
                            elementVector[local_ii] += lambda2Factor * lambda2Values[i];
                        }
                    }
                }
            }
            lambda2.unbind();
        }
    };

    template<class LC, class DV, class IS>
    struct GridFunctionOperatorRegistry<tag::surface_lambda1_lambda2<DV,IS>, LC> {
        static constexpr int degree = 0;
        using type = Surface_lambda1_lambda2<DV,IS>;
    };
    /** @} **/

} // end namespace AMDiS
