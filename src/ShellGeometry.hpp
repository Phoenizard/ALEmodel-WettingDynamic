//
// Created by Marcel Mokbel on 30.11.2020.
//

#pragma once

#include <algorithm>

#include <amdis/AMDiS.hpp>

#include <dune/grid/common/entity.hh>

using namespace AMDiS;

template<int dow>
struct CompareFieldVector {
    template<class T>
    static bool almost_equal(T x, T y, int ulp = 2) {
        return std::abs(x - y) <= std::numeric_limits<T>::epsilon() * std::abs(x + y) * ulp
               || std::abs(x - y) < std::numeric_limits<T>::min();
    }

    bool operator()(FieldVector<double, dow> const &lhs, FieldVector<double, dow> const &rhs) const {
        for (int i = 0; i < dow; i++) {
            if (almost_equal(lhs[i], rhs[i]))
                continue;
            return lhs[i] < rhs[i];
        }
        return false;
    }
};

/** @brief shell geometry class taking a previously extracted surface grid, which will be wrapped
 *         into a moving surface grid, representing the geometry of the surface,
 *         i.e. the coordinates of the surface DOFs and the surface normals and curvatures
 */
template<class SPB, class GV>
class ShellGeometry {
    static const int dow = GV::Grid::dimensionworld;
    using SurfaceGrid = typename SPB::SurfaceGrid;
    using GridView = typename SurfaceGrid::LeafGridView;
    using CoordinatesGF = Dune::DiscreteGridViewFunction<GridView, dow>;
    using Seed = typename GV::Grid::template Codim<0>::EntitySeed;
    using SeedSurface = typename SurfaceGrid::template Codim<0>::EntitySeed;


public:
    /// @brief constructor creating the shell geometry by using the previously extracted surface mesh
    ShellGeometry(SPB const &surfacePreBasis,
                  GV const& bulkGridView)
            : surfacePreBasis_(surfacePreBasis)
            , surfaceToFluidMapIdx_(surfacePreBasis.surfaceToFluidMapIdx())
            , surfaceToFluidMap_(surfacePreBasis.surfaceToFluidMap())
            , coordinates_(surfacePreBasis.surfaceGrid()->leafGridView(), 1)
            , surfaceGrid_(surfacePreBasis.surfaceGrid())
            , bulkGridView_(bulkGridView)
            {

        //fill coordinates_
        updateGridFunctions();
    }

    /// @brief Return coordinate grid-function
    auto &coordinates() {
        return coordinates_;
    }

    /// @brief update DiscreteGridViewFunctions, to be called if the grid has changed
    void updateGridFunctions() {
        coordinates_.update(surfaceGrid_->leafGridView());
        Dune::Functions::interpolate(coordinates_.basis(), coordinates_.coefficients(),
                                     [](auto const &x) { return x; });
        // update neighbor indices
        nIdx_ = findNeighborIndices();
    }

    /** @brief interpolation of values from surface to bulk
     * bulkVec = normalsVec * (factor * kappaVec)
     */
    template<class BulkVec, class NormalsVec, class CurvatureVec>
    void
    interpolateToBulk(BulkVec &bulkVec, NormalsVec const &normalsVec, double factor, CurvatureVec const &kappaVec) {
        for (auto it : surfaceToFluidMapIdx_) {
            bulkVec[it.second] = normalsVec[it.first] * (factor * kappaVec[it.first]);
        }
    }

    /** @brief interpolation of values from surface to bulk
     * bulkVec = shellVec
     */
    template<class BulkVec, class ShellVec>
    void interpolateToBulk(BulkVec &bulkVec, ShellVec const &shellVec) {
        for (auto it : surfaceToFluidMapIdx_) {
            bulkVec[it.second] = shellVec[it.first];
        }
    }

    /** @brief interpolation of values from bulk to surface
     */
    template<class BulkFct, class ShellFct>
    void interpolateToSurface(BulkFct const &bulkFct, ShellFct &shellFct) {
        for (auto it : surfaceToFluidMapIdx_) {
            shellFct[it.first] = bulkFct[it.second];
        }
    }

    /// @brief return the surface grid
    auto const &surfaceGrid() const {
        return surfaceGrid_;
    }

    /// @brief return the mapping from surface indices to bulk indices
    auto const &surfaceToFluidMapIdx() const {
        return surfaceToFluidMapIdx_;
    }

    /// @brief return the mapping from surface elements to bulk elements
    auto const &surfaceToFluidMap() const {
        return surfaceToFluidMap_;
    }
    /// @brief Return the gridView of the basis
    GridView gridView() const { return surfaceGrid_->leafGridView(); }

    /// @brief return the index map containing neighbor indices for all vertices of the grid
    auto const& neighborIndices() {
        return nIdx_;
    }

    /// @brief compute the normal direction to a given point on the surface
    template <class SEL>
    auto normalsDirection(SEL const& surfaceElement) {
        auto const& is = bulkGridView_.indexSet();
        auto facets = surfacePreBasis_.facets();
        for (const auto &map: surfaceToFluidMap_) {
            auto const &sEl = surfaceGrid_->entity(map.first);
            auto geoAux = surfaceElement.geometry();
            auto geo = sEl.geometry();
            auto x0 = geo.corner(0);
            auto x1 = geo.corner(1);
            auto x00 = geoAux.corner(0);
            auto x11 = geoAux.corner(1);

            if (((x0-x00).two_norm() < 1.e-10 && (x1-x11).two_norm() < 1.e-10)
                || ((x1-x00).two_norm() < 1.e-10 && (x0-x11).two_norm() < 1.e-10)) {
                auto const &bEl = bulkGridView_.grid().entity(map.second);
                auto face = facets[is.index(bEl)];

                // define geometry of another facet
                auto geo2 = bEl.template subEntity<1>((face + 1) % bEl.subEntities(1)).geometry();
                // create vector representing the facet
                auto x02 = geo2.corner(0);
                auto x12 = geo2.corner(1);
                auto normalDirection = x02 - x12;
                //change direction of normalDirection if x02 is on the surface
                if ((x02 - geo.corner(0)).two_norm() < 1.e-10 || (x02 - geo.corner(1)).two_norm() < 1.e-10)
                    normalDirection = -normalDirection;

                return normalDirection;
            }
        }
    }

    /// @brief find neighbor indices for all vertices stored in a vector of vectors
    /// for each index, the result has one or two indices referring to the neighboring DOFs of the vertex
    auto findNeighborIndices() {
        auto const& indexSet = gridView().indexSet();
        std::vector<std::vector<int>> neighborIdx(coordinates_.coefficients().size());

        // go through all vertices
        for (const auto& v : vertices(gridView())) {
            FieldVector<double,dow> normalDirection(0.0);
            auto midpoint = v.geometry().corner(0);

            // find all elements containing the vertex
            for (const auto& el : elements(gridView())) {
                for (int i = 0; i < el.geometry().corners(); i++) {
                    if ((el.geometry().corner(i) - midpoint).two_norm() < 1.e-5) {
                        //compute element normal
                        auto t = el.geometry().corner(1) - el.geometry().corner(0);
                        auto n = t;
                        n[0] = -t[1];
                        n[1] = t[0];
                        if (n.dot(normalsDirection(el)) < 0) n=-n;
                        // store point on the outside grid to find if the indices have to be switched
                        normalDirection += n;

                        // go through all vertices to find the respective neighbor index
                        for (const auto& v2 : vertices(gridView())) {
                            if ((v2.geometry().corner(i) - el.geometry().corner((i+1)%2)).two_norm() < 1.e-5) {
                                neighborIdx[indexSet.index(v)].push_back(indexSet.index(v2));
                            }
                        }
                    }
                }
            }
            switchNeighbors(normalDirection,neighborIdx[indexSet.index(v)]);
        }
        return neighborIdx;
    }

    /// @brief switch the neighbors indices to correct the direction, the normals point into and to correct
    /// the computation of all surface properties, where derivatives have to be computed
    void switchNeighbors(FieldVector<double,2> normalDirection,
                         std::vector<int>& neighborIndices) {
        auto const& indexSet = gridView().indexSet();
        // go through all vertices
        for (const auto& v : vertices(gridView())) {
            std::vector<FieldVector<double,dow>> neighbors(2);
            if (neighborIndices.size() == 2) {
                neighbors[0] = coordinates_.coefficients()[neighborIndices[0]];
                neighbors[1] = coordinates_.coefficients()[neighborIndices[1]];
                auto midpoint = v.geometry().corner(0);

                // compute point distances
                double dist1, dist2;
                computeDist1Dist2(neighbors,midpoint,dist1,dist2);

                // compute normal
                double Zs = (neighbors[1][0] - neighbors[0][0]) / (dist1 + dist2);
                double Rs = (neighbors[1][1] - neighbors[0][1]) / (dist1 + dist2);
                double Xs = sqrt(Zs * Zs + Rs * Rs);
                auto auxNormal = FieldVector<double,2>{-Rs/Xs,Zs/Xs};

                // see, if direction of normal has to be switched and hence, neighborIndices have to be switched
                if (auxNormal.dot(normalDirection) < 0) {
                    auto aux = neighborIndices[1];
                    neighborIndices[1] = neighborIndices[0];
                    neighborIndices[0] = aux;
                }
            }
        }
    }

    /// @brief compute distances of a given midpoint to its neigboring point(s)
    void computeDist1Dist2(std::vector<Dune::FieldVector<double,2>> neighbors,
                           Dune::FieldVector<double,2> midpoint,
                           double &dist1,
                           double &dist2){
        FieldVector<double,2> d1;
        auto d2 = neighbors[0] - midpoint;
        (neighbors.size() == 2) ? d1 = neighbors[1] - midpoint : d1 = d2;
        dist1 = d1.two_norm();
        dist2 = d2.two_norm();
    }

public:
    std::map<int, int> const &surfaceToFluidMapIdx_;
    std::vector<std::pair<SeedSurface, Seed>> const& surfaceToFluidMap_;   //map to create entity seeds of surface elements
    CoordinatesGF coordinates_;
    std::shared_ptr<SurfaceGrid> const& surfaceGrid_;
    std::vector<std::vector<int>> nIdx_;                            // indices of neigbor DOFs for each vertex
    GV const& bulkGridView_;
    SPB const& surfacePreBasis_;
};
