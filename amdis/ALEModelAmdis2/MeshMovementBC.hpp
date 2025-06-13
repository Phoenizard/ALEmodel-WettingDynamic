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
    class MeshMovementBC
    {
        using GridView = typename SubBasis::GridView;
        using LevelGridView = typename GridView::Grid::LevelGridView;
        using Intersection = typename GridView::Intersection;

        using Domain = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;
        using Range = TYPEOF(std::declval<ValueGridFct>()(std::declval<Domain>()));

    public:
        /// Make a DirichletBC from a global basis
        template <class Basis>
        MeshMovementBC(Basis&& basis,
                       RowPath const& row,
                       ColPath const& col,
                       std::vector<int> const& partitions,
                       ValueGridFct values,
                       int part = -1)
                : basis_{wrap_or_share(FWD(basis))}
                , row_{row}
                , col_{col}
                , partitions_(partitions)
                , values_(std::move(values))
                , level0GridView_(basis.gridView().grid().levelGridView(0))
                , part_(part)
        {}

        // Make a DirichletBC from a global basis
        template <class B>
        MeshMovementBC(B&& basis,
                       std::vector<int> const& partitions,
                       ValueGridFct values)
                : MeshMovementBC(FWD(basis),makeTreePath(),makeTreePath(),partitions,std::move(values))
        {}




        /// Fill \ref dirichletNodes_ with 1 or 0 whether DOF corresponds to the
        /// inner intersection or not.
        void init();

        /// \brief Apply dirichlet BC to matrix and vector
        template <class Mat, class Sol, class Rhs>
        void apply(Mat& matrix, Sol& solution, Rhs& rhs);

    private:
        std::shared_ptr<SubBasis const> basis_;
        LevelGridView level0GridView_;
        std::vector<int> const& partitions_;
        RowPath row_;
        ColPath col_;
        BoundarySubset<Intersection> boundarySubset_;
        ValueGridFct values_;
        std::vector<typename SubBasis::MultiIndex> rowIndices_;
        std::vector<typename SubBasis::MultiIndex> colIndices_;
        int part_; // partition number for BC to be set (-1 if both)
    };

    // Make a DirichletBC from a basis with treepath arguments
    template <class B, class Row, class Col, class Values, class Basis = Underlying_t<B>,
            REQUIRES(Concepts::GlobalBasis<Basis>)>
    MeshMovementBC(B const&, Row const&, Col const&, std::vector<int> const&, Values const&)
    -> MeshMovementBC<Basis, Row, Col, Underlying_t<Values>>;

    template <class B, class Row, class Col, class Values, class Basis = Underlying_t<B>,
            REQUIRES(Concepts::GlobalBasis<Basis>)>
    MeshMovementBC(B const&, Row const&, Col const&, std::vector<int> const&, Values const&, int)
    -> MeshMovementBC<Basis, Row, Col, Underlying_t<Values>>;

    // Make a DirichletBC from a global basis
    template <class B, class Values, class Basis = Underlying_t<B>,
            REQUIRES(Concepts::GlobalBasis<Basis>)>
    MeshMovementBC(B const&, std::vector<int> const&, Values const&)
    -> MeshMovementBC<Basis,TYPEOF(makeTreePath()),TYPEOF(makeTreePath()), Underlying_t<Values>>;

} // end namespace AMDiS

#include "MeshMovementBC.inc.hpp"
