#pragma once

#include <type_traits>

#include <dune/functions/functionspacebases/subentitydofs.hh>

#include <amdis/LinearAlgebra.hpp>

namespace AMDiS {

    template <class B, class Row, class Col, class V>
    void AxisymmetricDirichletBC<B, Row, Col, V>::init() {
       if (row_ == col_) {
            auto rowBasis = Dune::Functions::subspaceBasis(*basis_, col_);

            std::set<typename B::MultiIndex> rowSet;

            auto rowLocalView = rowBasis.localView();
            auto rowDOFs = Dune::Functions::subEntityDOFs(rowBasis);
            auto const &gridView = basis_->gridView();
            auto const &indexSet = gridView.indexSet();

            auto facets = basis_->preBasis().subPreBasis(_7).facets();

           auto line = Parameters::get<int>("line mesh").value_or(0);

           double ymin = 100000, ymax = -100000;
           if (line) {
               //find corner coordinates of the grid
               for (auto const &el: elements(gridView)) {
                   for (auto const &is: intersections(gridView, el)) {
                       if (is.neighbor()) //only consider boundary intersections
                           continue;
                       for (int i = 0; i < is.geometry().corners(); i++) {
                           ymin = std::min(ymin, is.geometry().corner(i)[1]);
                           ymax = std::max(ymax, is.geometry().corner(i)[1]);
                       }
                   }
               }
           }

            for (auto const &element: elements(gridView)) {
                int face = facets[indexSet.index(element)];

                if (face >= 0) {
                    auto myIntersection = element.template subEntity<1>(face);
                    rowLocalView.bind(element);
                    auto node = rowLocalView.tree();
                    for (int i = 0; i < myIntersection.geometry().corners(); ++i) {
                        //only add points on the shell that touch the symmetry axis
                        const auto localIndex = node.localIndex(i);
                        if (line ? (myIntersection.geometry().corner(i)[1] >= ymax - 1.e-7 || (!lineOnly_ && myIntersection.geometry().corner(i)[1] <= ymin + 1.e-7))
                                    : myIntersection.geometry().corner(i)[1] < 1.0e-10) {
                            rowSet.insert(rowLocalView.index(localIndex));
                        }
                    }
                }
            }

            rowIndices_.clear();
            rowIndices_.reserve(rowSet.size());
            rowIndices_.insert(rowIndices_.begin(), rowSet.begin(), rowSet.end());
            colIndices_ = rowIndices_;
        } else {
            auto rowBasis = Dune::Functions::subspaceBasis(*basis_, row_);
            auto colBasis = Dune::Functions::subspaceBasis(*basis_, col_);

            std::set<typename B::MultiIndex> rowSet, colSet;

            auto rowLocalView = rowBasis.localView();
            auto colLocalView = colBasis.localView();
            auto rowDOFs = Dune::Functions::subEntityDOFs(rowBasis);
            auto colDOFs = Dune::Functions::subEntityDOFs(colBasis);
            auto const &gridView = basis_->gridView();
            auto const &indexSet = gridView.indexSet();

            auto facets = basis_->preBasis().subPreBasis(_7).facets();

            for (auto const &element: entitySet(*basis_)) {
                int face = facets[indexSet.index(element)];

                if (face >=0) {
                    auto myIntersection = element.template subEntity<1>(face);
                    rowLocalView.bind(element);
                    colLocalView.bind(element);
                    auto rowNode = rowLocalView.tree();
                    auto colNode = colLocalView.tree();
                    for (int i = 0; i < myIntersection.geometry().corners(); ++i) {
                        const auto rowLocalIndex = rowNode.localIndex(i);
                        const auto colLocalIndex = colNode.localIndex(i);
                        if (myIntersection.geometry().corner(i)[1] < 1.0e-10) {
                            rowSet.insert(rowLocalView.index(rowLocalIndex));
                            colSet.insert(rowLocalView.index(colLocalIndex));
                        }
                    }
                }
            }

            rowIndices_.clear();
            rowIndices_.reserve(rowSet.size());
            rowIndices_.insert(rowIndices_.begin(), rowSet.begin(), rowSet.end());

            colIndices_.clear();
            colIndices_.reserve(colSet.size());
            colIndices_.insert(colIndices_.begin(), colSet.begin(), colSet.end());
        }

    }


    template <class B, class Row, class Col, class V>
    template <class Mat, class Sol, class Rhs>
    void AxisymmetricDirichletBC<B, Row, Col, V>::apply(Mat& matrix, Sol& solution, Rhs& rhs)
    {
        std::vector<typename Sol::value_type> solutionValues;
        solutionValues.reserve(colIndices_.size());
        {
            // create a copy the solution, because the valueGridFct_ might involve the solution vector.
            Sol boundarySolution{solution};
            valueOf(boundarySolution, col_).interpolate_noalias(values_, tag::assign{});
            boundarySolution.gather(colIndices_,solutionValues);
        }

        // copy boundary solution dirichlet data to solution vector
        solution.scatter(colIndices_,solutionValues);

        matrix.zeroRows(rowIndices_, true);

        // copy boundary solution dirichlet data to rhs vector
        rhs.scatter(rowIndices_,solutionValues);
    }

} // end namespace AMDiS
