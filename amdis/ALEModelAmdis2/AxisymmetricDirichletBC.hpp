#pragma once

#include <functional>
#include <utility>
#include <vector>

#include <dune/common/concept.hh>
#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

#include <amdis/Boundary.hpp>
#include <amdis/BoundaryCondition.hpp>
#include <amdis/BoundarySubset.hpp>
#include <amdis/common/Concepts.hpp>
#include <amdis/common/TypeTraits.hpp>
#include <amdis/typetree/RangeType.hpp>
#include <amdis/typetree/TreePath.hpp>

namespace AMDiS
{
    /// Implements a boundary condition of Dirichlet-type for inner intersections
    template <class SubBasis, class RowPath, class ColPath, class ValueGridFct>
    class AxisymmetricDirichletBC
    {
        using GridView = typename SubBasis::GridView;
        using LevelGridView = typename GridView::Grid::LevelGridView;
        using Intersection = typename GridView::Intersection;

        using Domain = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;
        using Range = TYPEOF(std::declval<ValueGridFct>()(std::declval<Domain>()));

    public:
        /// Make a DirichletBC from a global basis
        template <class Basis>
        AxisymmetricDirichletBC(Basis&& basis,
                       RowPath const& row,
                       ColPath const& col,
                       ValueGridFct values,
                       int lineOnly = 0)
                : basis_{wrap_or_share(FWD(basis))}
                , row_{row}
                , col_{col}
                , values_(std::move(values))
                , lineOnly_(lineOnly)
        {}

        // Make a DirichletBC from a global basis
        template <class B>
        AxisymmetricDirichletBC(B&& basis,
                       ValueGridFct values)
                : AxisymmetricDirichletBC(FWD(basis),makeTreePath(),makeTreePath(),std::move(values))
        {}




        /// Fill \ref dirichletNodes_ with 1 or 0 whether DOF corresponds to the
        /// inner intersection or not.
        void init();

        /// \brief Apply dirichlet BC to matrix and vector
        template <class Mat, class Sol, class Rhs>
        void apply(Mat& matrix, Sol& solution, Rhs& rhs);

    private:
        std::shared_ptr<SubBasis const> basis_;
        RowPath row_;
        ColPath col_;
        BoundarySubset<Intersection> boundarySubset_;
        ValueGridFct values_;
        std::vector<typename SubBasis::MultiIndex> rowIndices_;
        std::vector<typename SubBasis::MultiIndex> colIndices_;
        int lineOnly_;
    };

    // Make a DirichletBC from a basis with treepath arguments
    template <class B, class Row, class Col, class Values, class Basis = Underlying_t<B>,
            REQUIRES(Concepts::GlobalBasis<Basis>)>
    AxisymmetricDirichletBC(B const&, Row const&, Col const&, Values const&)
    -> AxisymmetricDirichletBC<Basis, Row, Col, Underlying_t<Values>>;

    // Make a DirichletBC from a basis with treepath arguments
    template <class B, class Row, class Col, class Values, class Basis = Underlying_t<B>,
            REQUIRES(Concepts::GlobalBasis<Basis>)>
    AxisymmetricDirichletBC(B const&, Row const&, Col const&, Values const&, int)
    -> AxisymmetricDirichletBC<Basis, Row, Col, Underlying_t<Values>>;


    // Make a DirichletBC from a global basis
    template <class B, class Values, class Basis = Underlying_t<B>,
            REQUIRES(Concepts::GlobalBasis<Basis>)>
    AxisymmetricDirichletBC(B const&, Values const&)
    -> AxisymmetricDirichletBC<Basis,TYPEOF(makeTreePath()),TYPEOF(makeTreePath()), Underlying_t<Values>>;

} // end namespace AMDiS

#include "AxisymmetricDirichletBC.inc.hpp"
