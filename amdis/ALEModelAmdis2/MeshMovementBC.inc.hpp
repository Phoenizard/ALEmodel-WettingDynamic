#pragma once

#include <type_traits>

#include <dune/functions/functionspacebases/subentitydofs.hh>

#include <amdis/LinearAlgebra.hpp>

namespace AMDiS {

    template <class B, class Row, class Col, class V>
    void MeshMovementBC<B, Row, Col, V>::init() {
       if (row_ == col_) {
            auto rowBasis = Dune::Functions::subspaceBasis(*basis_, col_);

            std::set<typename B::MultiIndex> rowSet;

            auto rowLocalView = rowBasis.localView();
            auto rowDOFs = Dune::Functions::subEntityDOFs(rowBasis);
            auto const &gridView = basis_->gridView();
            auto const &indexSet = level0GridView_.indexSet();
            for (auto const &element: entitySet(*basis_)) {
                auto father = element;
                while (father.hasFather())
                    father = father.father();
                int p = partitions_[indexSet.index(father)];

                rowLocalView.bind(element);
                for (auto const &intersection: intersections(gridView, element)) {
                    if (!intersection.neighbor())
                        continue;
                    auto fatherOutside = intersection.outside();
                    while (fatherOutside.hasFather())
                        fatherOutside = fatherOutside.father();
                    int q = partitions_[indexSet.index(fatherOutside)];
                    if (p < 0 || q < 0)
                        continue;
                    if ((part_ == -1 && p != q) || (part_ == 1 && p == 1 && q == 0)  || (part_ == 0 && p == 0 && q == 1)) {
                            for (auto localIndex: rowDOFs.bind(rowLocalView, intersection))
                                rowSet.insert(rowLocalView.index(localIndex));
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
            auto const &indexSet = level0GridView_.indexSet();
            for (auto const &element: entitySet(*basis_)) {
                auto father = element;
                while (father.hasFather())
                    father = father.father();
                int p = partitions_[indexSet.index(father)];
                rowLocalView.bind(element);
                colLocalView.bind(element);
                for (auto const &intersection: intersections(gridView, element)) {
                    if (!intersection.neighbor())
                        continue;
                    auto fatherOutside = intersection.outside();
                    while (fatherOutside.hasFather())
                        fatherOutside = fatherOutside.father();
                    int q = partitions_[indexSet.index(fatherOutside)];
                    if (p != q && p >= 0 && q >= 0) {
                        for (auto localIndex: rowDOFs.bind(rowLocalView, intersection))
                            rowSet.insert(rowLocalView.index(localIndex));
                        for (auto localIndex: colDOFs.bind(colLocalView, intersection))
                            colSet.insert(colLocalView.index(localIndex));
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
    void MeshMovementBC<B, Row, Col, V>::apply(Mat& matrix, Sol& solution, Rhs& rhs)
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
