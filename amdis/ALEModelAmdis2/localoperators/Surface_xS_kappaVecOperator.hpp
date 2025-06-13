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
        template<class NV, class IS = int>
        struct surface_xS_kappaVec {
            unsigned int axi = 0; // axisymmetric or not?
            double const& tau; // the time step size
            int normalMovement = 0; // move grid points with velocity (0) or with normal velocity (1)?
            double perm = 0; // membrane permeability
            double p0 = 0; // pressure difference of a closed shell in its equilibrium shape
            double kappa0 = 0; // spontaneous curvature
            NV const& normalVec = {}; // unit outer normal to the membrane
            std::size_t inside = 0; // integrate over the shell from inside (0) or outside (1)
            IS const &indexSet = {}; //the indexSet of the leafGridView of the grid!
            std::vector<int> const &partitions = {}; //partition numbers of each element
            std::vector<int> const &facets = {}; //facet number of the surface facet (-1 if element i not touching the surface)
        };        
    }


    /// An operator for the computation of the curvature vector \mathbf{\kappa} = \kappa\mathbf{n} using Dzuik's formula.
    /// The operator computes the coordinates of the grid points, when the grid is moved with fluid velocity and then
    /// computes \mathbf{\kappa} with the Laplace-Beltrami of the new coordinates
    /// The operator includes the possibility to include membrane permeability to the model
    template<class NV, class IS = int>
    class Surface_xS_kappaVec {
        unsigned int axi_;
        double const& tau_;
        int normalMovement_ = 0;
        double perm_ = 0;
        double p0_ = 0;
        double kappa0_ = 0;
        NV const& normalVec_;
        IS const &indexSet_ = {};
        std::vector<int> const &partitions_;
        std::vector<int> const &facets_;
        std::size_t inside_;
    public:
        Surface_xS_kappaVec(tag::surface_xS_kappaVec<NV,IS> t)
                : axi_(t.axi), tau_(t.tau), normalMovement_(t.normalMovement), 
                  perm_(t.perm),  p0_(t.p0), kappa0_(t.kappa0), normalVec_(FWD(t.normalVec)), indexSet_(t.indexSet), partitions_(t.partitions), facets_(t.facets), inside_(t.inside) {}

        template<class CG, class Node, class Quad, class LocalFct, class Mat>
        void assemble(CG const &contextGeo, Node const &tree, Node const /*colNode*/,
                      Quad const &quad, LocalFct const &localFct, Mat &elementMatrix) const {
            static_assert(static_size_v<typename LocalFct::Range> == 1,
                          "Expression must be of scalar type.");

            using VelocityNode = Dune::TypeTree::Child<Node, _vS.value>;
            using PNode = Dune::TypeTree::Child<Node, _p.value>;
            using XSNode = Dune::TypeTree::Child<Node, _xS.value>;
            using KappaVecNode = Dune::TypeTree::Child<Node, _kappaVec.value>;


            static_assert(VelocityNode::isPower && PNode::isLeaf && XSNode::isPower && KappaVecNode::isPower,
                          "Nodes must have correct dimensions.");

            auto const &vSNode = tree.child(_vS);
            auto const &pNode = tree.child(_p);
            auto const &xSNode = tree.child(_xS);
            auto const &kappaVecNode = tree.child(_kappaVec);

            if (pNode.size() == 0 || xSNode.child(0).size() == 0 || vSNode.child(0).size() == 0)
                return;

            // compute facet number and partition number from surface basis node (if no surface basis node: return -2)
            auto part = std::max(partition(xSNode.child(0), Dune::PriorityTag<2>{}),
                                 partition(vSNode.child(0), Dune::PriorityTag<2>{}));
            auto face = std::max(facet(xSNode.child(0), Dune::PriorityTag<2>{}),
                                 facet(vSNode.child(0), Dune::PriorityTag<2>{}));

            if (face == -2) {
                auto index = idx(indexSet_, contextGeo.element(), Dune::PriorityTag<2>{});
                test_exit(!(partitions_.size() == 0 || facets_.size() == 0 || index == -1),
                          "please provide partitions, facets and indexSet of the leafGridView when test and trial "
                          "functions are not from the surfaceLagrange basis!");
                face = facets_[index];
                part = partitions_[index];
            }

            if (face == -1 || (part != inside_ && std::abs(perm_) < 1.e-10)) // do nothing if element has no facet on the shell
                return;

            std::size_t xSSize = xSNode.child(0).size();
            std::size_t kappaVecSize = kappaVecNode.child(0).size();
            std::size_t vSSize = vSNode.child(0).size();
            std::size_t pSize = pNode.size();

            using RangeFieldType = typename KappaVecNode::ChildType::LocalBasis::Traits::RangeFieldType;
            using WorldVector = Dune::FieldVector<RangeFieldType, CG::dow>;
            std::vector<WorldVector> gradientsKappaVec, gradientsKappaVecAux;

            using XSFieldType = typename XSNode::ChildType::LocalBasis::Traits::RangeFieldType;
            using XSWorldVector = Dune::FieldVector<XSFieldType, CG::dow>;
            std::vector<XSWorldVector> gradientsXS, gradientsXSAux;

            auto kappaVecDim = basisDim(kappaVecNode.child(0), Dune::PriorityTag<2>{}); // = dow-1 if surfaceLagrange node, = dow else
            auto xSDim = basisDim(xSNode.child(0), Dune::PriorityTag<2>{});

            auto element = contextGeo.element();

            auto intersection = element.template subEntity<1>(face);

            constexpr auto intersectionDim = std::decay_t<decltype(intersection)>::mydimension;
            const auto &quadShell = Dune::QuadratureRules<double, intersectionDim>::rule(intersection.type(), 2);

            // apply the operator only to intersections on the shell
            auto geo = intersection.geometry();

            auto normalX = localFunction(valueOf(normalVec_[0]));
            auto normalY = localFunction(valueOf(normalVec_[1]));
            normalX.bind(contextGeo.element());
            normalY.bind(contextGeo.element());

            //todo: make 3D ready
            auto n = geo.corner(0);

            for (auto const &qp: quadShell) {
                auto localGeo = referenceElement(element).template geometry<1>(face);

                // Position of the current quadrature point in the reference element
                auto &&local = localGeo.global(qp.position());
                auto &&global = contextGeo.elementGeometry().global(local);

                // The multiplicative factor in the integral transformation formula
                const auto factor = geo.integrationElement(qp.position()) * qp.weight();
                const auto permeabilityFactor = perm_ * localFct(local);

                n[0] = normalX(local);
                n[1] = normalY(local);

                // The transposed inverse Jacobian of the map from the reference facet to the real facet
                const auto jacobian = geo.jacobianInverseTransposed(qp.position());

                // The transposed inverse Jacobian of the map from the reference facet to the reference element
                const auto jacobianRef = localGeo.jacobianTransposed(qp.position());

                // The gradients of the shape functions on the reference element
                auto const &shapeGradientsKappaVec = kappaVecNode.child(0).localBasisJacobiansAt(local);
                auto const &shapeGradientsXS = xSNode.child(0).localBasisJacobiansAt(local);

                gradientsKappaVec.resize(shapeGradientsKappaVec.size());
                gradientsXS.resize(shapeGradientsXS.size());
                gradientsKappaVecAux.resize(shapeGradientsKappaVec.size());
                gradientsXSAux.resize(shapeGradientsXS.size());

                // transform local gradients to reference intersection if necessary and the compute
                // gradients on the real intersection
                int dow = CG::dow;
                if (kappaVecDim == dow) {
                    for (std::size_t i = 0; i < gradientsKappaVec.size(); ++i)
                        jacobianRef.mv(shapeGradientsKappaVec[i][0], gradientsKappaVecAux[i]);

                    for (std::size_t i = 0; i < gradientsKappaVec.size(); ++i)
                        jacobian.mv(gradientsKappaVecAux[i], gradientsKappaVec[i]);
                } else {
                    for (std::size_t i = 0; i < gradientsKappaVec.size(); ++i)
                        jacobian.mv(shapeGradientsKappaVec[i][0], gradientsKappaVec[i]);
                }
                if (xSDim == dow) {
                    for (std::size_t i = 0; i < gradientsXS.size(); ++i)
                        jacobianRef.mv(shapeGradientsXS[i][0], gradientsXSAux[i]);

                    for (std::size_t i = 0; i < gradientsXS.size(); ++i)
                        jacobian.mv(gradientsXSAux[i], gradientsXS[i]);
                } else {
                    for (std::size_t i = 0; i < gradientsXS.size(); ++i)
                        jacobian.mv(shapeGradientsXS[i][0], gradientsXS[i]);
                }
                
                auto const &xSValues = xSNode.child(0).localBasisValuesAt(local);
                auto const &kappaVecValues = kappaVecNode.child(0).localBasisValuesAt(local);
                auto const &vSValues = vSNode.child(0).localBasisValuesAt(local);
                auto const &pValues = pNode.localBasisValuesAt(local);
                auto axiValue = (axi_ * Dune::at(global,1) + 1.0 - axi_);

                if (part == inside_) {
                    // xS = ...
                    for (std::size_t i = 0; i < xSSize; ++i) {
                        const auto value = factor * xSValues[i];

                        for (std::size_t k = 0; k < CG::dow; ++k) {
                            const auto local_ki = xSNode.child(k).localIndex(i);
                            elementMatrix[local_ki][local_ki] += value * xSValues[i];
                        }

                        for (std::size_t j = i + 1; j < xSSize; ++j) {
                            const auto value0 = value * xSValues[j];

                            for (std::size_t k = 0; k < xSNode.degree(); ++k) {
                                const auto local_ki = xSNode.child(k).localIndex(i);
                                const auto local_kj = xSNode.child(k).localIndex(j);

                                elementMatrix[local_ki][local_kj] += value0;
                                elementMatrix[local_kj][local_ki] += value0;
                            }
                        }
                    }

                    // = tau * v or = tau*(v*n)*n
                    for (std::size_t i = 0; i < vSSize; ++i) {
                        for (std::size_t l = 0; l < CG::dow; ++l) {
                            const auto local_li = vSNode.child(l).localIndex(i);
                            const auto value = -tau_ * (factor * vSValues[i]);
                            for (std::size_t j = 0; j < xSSize; ++j) {
                                if (normalMovement_) {
                                    for (std::size_t k = 0; k < CG::dow; ++k) {
                                        const auto local_kj = xSNode.child(k).localIndex(j);
                                        elementMatrix[local_kj][local_li] += (value * Dune::at(n,l)) * Dune::at(n,k) * xSValues[j];
                                    }
                                } else {
                                    const auto local_lj = xSNode.child(l).localIndex(j);
                                    elementMatrix[local_lj][local_li] += value * xSValues[j];

                                }
                            }
                        }
                    }

                    if (std::abs(perm_) > 1.e-10) {
                        // permeability term from outside
                        for (std::size_t i = 0; i < xSSize; ++i) {
                            const auto value = -tau_ * permeabilityFactor * factor * xSValues[i];
                            for (std::size_t k = 0; k < CG::dow; ++k) {
                                const auto local_ki = xSNode.child(k).localIndex(i);

                                for (std::size_t j = 0; j < pSize; ++j) {
                                    const auto local_j = pNode.localIndex(j);
                                    elementMatrix[local_ki][local_j] += value * Dune::at(n,k) * pValues[j];
                                }
                            }
                        }
                    }

                    // kappaVec = \Delta_\Gamma xS
                    for (std::size_t i = 0; i < kappaVecSize; ++i) {
                        const auto value = axiValue * factor * kappaVecValues[i];
                        const auto valueAxiTerm = 1.0 / Dune::at(global,1) * factor * kappaVecValues[i];

                        for (std::size_t k = 0; k < CG::dow; ++k) {
                            const auto local_ki = kappaVecNode.child(k).localIndex(i);
                            elementMatrix[local_ki][local_ki] += value * kappaVecValues[i];

                            // <gradtest, gradtrial> with kappaVec and xS
                            for (std::size_t j = 0; j < xSSize; ++j) {
                                const auto local_kj = xSNode.child(k).localIndex(j);
                                //special implementation for axisymmetry at the symmetry axis
                                if (axi_ && Dune::at(geo.corner(i),1) < 1.e-8) {
                                    elementMatrix[local_ki][local_kj] += eval(Dune::at(global,1), 4.0/3.0*factor,
                                                                              gradientsKappaVec[i], gradientsXS[j]);
                                }
                                else {
                                    elementMatrix[local_ki][local_kj] += eval(axiValue, factor, gradientsKappaVec[i],
                                                                              gradientsXS[j]);
                                }
                                if (k == 1 && axi_) {
                                    elementMatrix[local_ki][local_kj] += valueAxiTerm * xSValues[j];
                                }
                            }
                        }

                        for (std::size_t j = i + 1; j < kappaVecSize; ++j) {
                            const auto value0 = value * kappaVecValues[j];

                            for (std::size_t k = 0; k < CG::dow; ++k) {
                                const auto local_ki = kappaVecNode.child(k).localIndex(i);
                                const auto local_kj = kappaVecNode.child(k).localIndex(j);

                                elementMatrix[local_ki][local_kj] += value0;
                                elementMatrix[local_kj][local_ki] += value0;
                            }
                        }
                    }
                } else {
                    // permeability term from inside
                    for (std::size_t i = 0; i < xSSize; ++i) {
                        const auto value = tau_ * permeabilityFactor * factor * xSValues[i];
                        for (std::size_t k = 0; k < CG::dow; ++k) {
                            const auto local_ki = xSNode.child(k).localIndex(i);

                            for (std::size_t j = 0; j < pSize; ++j) {
                                const auto local_j = pNode.localIndex(j);
                                elementMatrix[local_ki][local_j] += value * Dune::at(n,k) * pValues[j];
                            }
                        }
                    }
                }
            }
            normalX.unbind();
            normalY.unbind();
        }
        
        //assemble method for RHS
        template <class CG, class Node, class Quad, class LocalFct, class Vec>
        void assemble(CG const& contextGeo, Node const& tree, Quad const& quad,
                      LocalFct const& localFct, Vec& elementVector) const
        {
            static_assert(static_size_v<typename LocalFct::Range> == 1,
                          "Expression must be of scalar type." );

            using XSNode = Dune::TypeTree::Child<Node,_xS.value>;
            using KappaVecNode = Dune::TypeTree::Child<Node,_kappaVec.value>;

            static_assert(XSNode::isPower && KappaVecNode::isPower,
                          "Operator can only be applied to power nodes");

            auto const& xSNode = tree.child(_xS);
            auto const& kappaVecNode = tree.child(_kappaVec);

            std::size_t xSSize = xSNode.child(0).size();
            std::size_t kappaVecSize = kappaVecNode.child(0).size();

            if (xSNode.child(0).size()==0)
                return;

            // compute facet number and partition number from surface basis node (if no surface basis node: return -2)
            auto part = partition(xSNode.child(0), Dune::PriorityTag<2>{});
            auto face = facet(xSNode.child(0), Dune::PriorityTag<2>{});

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
            const auto& quadShell = Dune::QuadratureRules<double,intersectionDim>::rule(intersection.type(), 2);

            // apply the operator only to intersections on the shell
            if (part == inside_) {
                auto geo = intersection.geometry();

                auto normalX = localFunction(valueOf(normalVec_[0]));
                auto normalY = localFunction(valueOf(normalVec_[1]));
                normalX.bind(contextGeo.element());
                normalY.bind(contextGeo.element());

                //todo: make 3D ready
                auto n = geo.corner(0);

                for (auto const& qp : quadShell) {
                    auto localGeo = referenceElement(element).template geometry<1>(face);

                    // Position of the current quadrature point in the reference element
                    auto&& local = localGeo.global(qp.position());
                    auto &&global = contextGeo.elementGeometry().global(local);

                    // The multiplicative factor in the integral transformation formula
                    const auto factor = geo.integrationElement(qp.position()) * qp.weight();
                    const auto exprValue = -tau_ * p0_ * perm_ * localFct(local);

                    n[0] = normalX(local);
                    n[1] = normalY(local);

                    auto const& xSValues = xSNode.child(0).localBasisValuesAt(local);
                    auto const& kappaVecValues = kappaVecNode.child(0).localBasisValuesAt(local);

                    for (std::size_t i = 0; i < xSSize; ++i) { //assuming xSSize = kappaVecSize
                        const auto value = factor * xSValues[i];
                        const auto valueKappa = kappa0_ * (axi_ * Dune::at(global,1) + 1.0 - axi_) * factor * kappaVecValues[i];
                        for (std::size_t k = 0; k < CG::dow; ++k) {
                            const auto local_ki = xSNode.child(k).localIndex(i);
                            elementVector[local_ki] += value * (exprValue * Dune::at(n,k) + Dune::at(global,k));
                            if (xSSize == kappaVecSize && kappa0_) {
                                const auto local_kii = kappaVecNode.child(k).localIndex(i);
                                elementVector[local_kii] += -valueKappa * Dune::at(n,k);
                            }
                        }
                    }
                    if (xSSize != kappaVecSize && kappa0_) {
                        for (std::size_t i = 0; i < kappaVecSize; ++i) { //assuming xSSize = kappaVecSize
                            const auto valueKappa = kappa0_ * (axi_ * Dune::at(global,1) + 1.0 - axi_) * factor * kappaVecValues[i];
                            for (std::size_t k = 0; k < CG::dow; ++k) {
                                const auto local_kii = kappaVecNode.child(k).localIndex(i);
                                elementVector[local_kii] += -valueKappa * Dune::at(n,k);
                            }
                        }
                    }
                }
                normalX.unbind();
                normalY.unbind();
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

    template<class LC, class NV, class IS>
    struct GridFunctionOperatorRegistry<tag::surface_xS_kappaVec<NV,IS>, LC> {
        static constexpr int degree = 0;
        using type = Surface_xS_kappaVec<NV,IS>;
    };
    
    /** @} **/

} // end namespace AMDiS
