//
// Created by Marcel Mokbel on 22.06.2021.
//

#ifndef ALEMODELAMDIS2_MESHGENERATOR_H
#define ALEMODELAMDIS2_MESHGENERATOR_H
#include <set>
#include <unordered_set>
#include <cmath>
#include <gmsh.h>
#include <iostream>
#include <iomanip>


using Coord = FieldVector<double, 2>; // only 2D for now, 3D to be implemented
using Edge = std::pair<Coord, Coord>;
using PhiPair = std::pair<double, double>;
using BoundaryEdge = std::tuple<Edge, int, int, PhiPair>; //contains Edge, boundaryId, partition number and phi values in corners

constexpr double eps = 1e-12;

bool equalCoord(const Coord& a, const Coord& b) {
    return std::abs(a[0] - b[0]) < eps && std::abs(a[1] - b[1]) < eps;
}

struct CoordHash {
    std::size_t operator()(const Coord& c) const {
        std::size_t h1 = std::hash<int>()(static_cast<int>(c[0] / eps));
        std::size_t h2 = std::hash<int>()(static_cast<int>(c[1] / eps));
        return h1 ^ (h2 << 1);
    }
};

struct CoordEqual {
    bool operator()(const Coord& a, const Coord& b) const {
        return equalCoord(a, b);
    }
};

using NeighborMap = std::unordered_map<Coord, std::vector<std::pair<Coord, int>>, CoordHash, CoordEqual>;
using BCIDMap = std::unordered_map<Coord, std::pair<std::set<int>,double>, CoordHash, CoordEqual>;
using LevelMap = std::unordered_map<Coord, std::set<int>, CoordHash, CoordEqual>;
using PointSet = std::unordered_set<Coord,CoordHash,CoordEqual>;
using PointMapGmsh = std::unordered_map<Coord, int, CoordHash, CoordEqual>;
using LineMapGmsh = std::map<std::pair<int, int>, int>;
using MyGraph = std::unordered_map<int, std::vector<std::pair<int, std::pair<int, int>>>>;

template <class GV>
double maximumBoundaryFacetLength(GV const& gridView) {
    double distMax = 0;
    for (auto const &el: elements(gridView)) {
        for (auto const &is: intersections(gridView, el)) {
            if (is.neighbor()) //only consider boundary intersections
                continue;

            auto dist = two_norm(is.geometry().corner(1) - is.geometry().corner(0));
            if (dist > distMax) distMax = dist;
        }
    }
    return distMax;
}

/** extracts all edges on the boundary of the grid and on the membrane
 * stored in a vector containing
 *         1. the grid points of the facet,
 *         2. the boundary number of the facet (-1 if facet on the membrane)
 *         3. the partition number of the element containing this facet (-1 if element on the membrane)
 *         4. the values of phi in the grid points
 * currently only implemented in 2D
 */
template <class NschProb>
std::vector<BoundaryEdge> extractBoundaryEdges(AdaptInfo& adaptInfo, NschProb& nschProb)
{
    std::vector<BoundaryEdge> boundaryEdges;
    std::unordered_set<int> visitedSegments;
    auto const& gridView =  nschProb.gridView();
    auto& indexSet = gridView.indexSet();
    auto& indexSet0 = gridView.grid().levelGridView(0).indexSet();

    auto lfPhi = localFunction(nschProb.getPhase()); // the local values of the phase field inside the element

    for (const auto& element : elements(gridView)) {
        lfPhi.bind(element);
        for (auto const& is : intersections(gridView,element)) {
            auto father = is.inside();
            while (father.hasFather())
                father = father.father();

            // is boundary on inside or outside partition of the grid? -> extract partitions number
            int partition = nschProb.partitions()[indexSet0.index(father)];

            if (is.neighbor()) {
                auto father = is.outside();
                while (father.hasFather())
                    father = father.father();
                int q = nschProb.partitions()[indexSet0.index(father)];
                if (!(partition == 0 && q == 1))
                    continue;
                else
                    partition = -1;
            } // now the intersection is either on the boundary (is.neighbor() = false), or on the membrane

            // find boundary ID (on the membrane, set ID to -1)
            int bcId = -1;
            if (is.boundary()) bcId = nschProb.boundaryManager()->boundaryId(is);

            Coord p1 = is.geometry().corner(0);
            Coord p2 = is.geometry().corner(1);

            auto local1 = element.geometry().local(p1);
            auto local2 = element.geometry().local(p2);

            auto valuePhi1 = lfPhi(local1);
            auto valuePhi2 = lfPhi(local2);

            // only extract level 0 facets on outer boundaries
            if (bcId > -1 && !visitedSegments.count(is.boundarySegmentIndex())) {
                for (auto const &isFather: intersections(gridView.grid().levelGridView(0), father)) {
                    if (isFather.neighbor()) continue;

                    p1 = isFather.geometry().corner(0);
                    p2 = isFather.geometry().corner(1);

                    visitedSegments.insert(is.boundarySegmentIndex());
                    boundaryEdges.emplace_back(Edge{p1, p2}, bcId, partition, PhiPair{valuePhi1,valuePhi2});
                }
                continue;
            } else if (bcId > -1 && visitedSegments.count(is.boundarySegmentIndex())) {
                //already visited - do nothing
            } else {
                boundaryEdges.emplace_back(Edge{p1, p2}, bcId, partition, PhiPair{valuePhi1, valuePhi2});
            }
        }
        lfPhi.unbind();
    }

    return boundaryEdges;
}

/**
 * Compute the mean length of the boundary edges on the interface (0.05<phi<0.95) and
 * for the membrane in the bulk (phi < 0.025 || phi > 0.975)
 * @param boundaryEdges the vector storing the boundary edge information
 * @param distMeanIface the mean edge length on the interface to be computed
 * @param distMeanBulk the mean edge length in the bulk on the membrane to be computed
 */
void meanEdgeLengths(std::vector<BoundaryEdge> const& boundaryEdges,
                     double& distMeanIface,
                     double& distMeanBulk) {
    int count = 0, count2 = 0;
    for (const auto& [edge, bcId, partition,  phiValues] : boundaryEdges) {
        const Coord& p1 = edge.first;
        const Coord& p2 = edge.second;

        auto phaseValue = phiValues.first;
        auto phaseValue2 = phiValues.second;
        if ((phaseValue < 0.95 && phaseValue > 0.05) && (phaseValue2 < 0.95 && phaseValue2 > 0.05)) {
            distMeanIface += (p1 - p2).two_norm();
            count++;
        }

        if (bcId == -1) {
            if ((phaseValue > 0.975 || phaseValue < 0.025) && (phaseValue2 > 0.975 || phaseValue2 < 0.025)) {
                distMeanBulk += (p1 - p2).two_norm();
                count2++;
            }
        }
    }
    distMeanIface /= count; // mean length of a interface facet
    distMeanBulk /= count2; // mean length of a bulk facet
}

/**
 * Compute the locations on the boundary, where the interface is (phi = 0.5)
 * @param iFaceLocations the mean edge length in the bulk on the membrane to be computed
 */
std::vector<Coord> iFaceLocations(std::vector<BoundaryEdge> const& boundaryEdges) {
    std::vector<Coord> iFaceLocations;
    for (const auto& [edge, bcId, partition, phiValues] : boundaryEdges) {
        const Coord& p1 = edge.first;
        const Coord& p2 = edge.second;

        auto phaseValue = phiValues.first;
        auto phaseValue2 = phiValues.second;
        if ((phaseValue < 0.5 && phaseValue2 >= 0.5) || (phaseValue2 < 0.5 && phaseValue >= 0.5)) {
            iFaceLocations.push_back(0.5 * (p1 + p2));
        }
    }
    return iFaceLocations;
}

void buildNeighborhoods(std::vector<BoundaryEdge> const& boundaryEdges,
                        NeighborMap& pointToNeighbors,
                        BCIDMap& pointToBcIds)
{
    pointToNeighbors.clear();
    pointToBcIds.clear();

    for (const auto& [edge, bcId, partition,  phiValues] : boundaryEdges) {
        const Coord& p1 = edge.first;
        const Coord& p2 = edge.second;

        pointToNeighbors[p1].emplace_back(p2, bcId);
        pointToNeighbors[p2].emplace_back(p1, bcId);

        pointToBcIds[p1].first.insert(bcId);
        pointToBcIds[p2].first.insert(bcId);
        pointToBcIds[p1].second = phiValues.first;
        pointToBcIds[p2].second = phiValues.second;
    }
}

PointSet markPointsForRemoval(NeighborMap const& pointToNeighbors,
                          BCIDMap& pointToBcIds,
                          std::vector<Coord> iFace,
                          double distMeanBulk,
                          double distMeanIface) {
    PointSet pointsToRemove;
    for (const auto &[p, neighbors]: pointToNeighbors) {
        if (neighbors.size() != 2) continue; // only remove points with exacttly two neighbors

        if (pointToBcIds[p].first.size() > 1) continue; // multiple bcIDs -> do not remove

        const Coord &p1 = neighbors[0].first;
        const Coord &p2 = neighbors[1].first;


        double d1 = (p - p1).two_norm();
        double d2 = (p - p2).two_norm();
        double phiVal = pointToBcIds[p].second;

        auto maxDist = Parameters::get<double>("distance for removal of interface points").value_or(0);
        maxDist *= distMeanBulk;

        // distance of p to the interface of phi
        double distIface = 100000000000;
        for (int k = 0; k < iFace.size(); k++) {
            if (two_norm(p - iFace[k]) < distIface) distIface = two_norm(p - iFace[k]);
        }

        auto bcId = *pointToBcIds[p].first.begin();

        // booleans
        bool inInterface = (phiVal >= 0.05 && phiVal <= 0.95);
        bool farFromInterface = (distIface > 1.5 * distMeanBulk) && (phiVal < 0.05 || phiVal > 0.95);
        bool isCloseInterface = (d1 < distMeanIface / 3.0 || d2 < distMeanIface / 3.0);
        bool isClose = ((d1 + d2) * 0.5 < maxDist);

        // remove points if 1: far from interface and a) close enough or b) points are on refined elements
        // or 2: in interface and close enough
        // on external boundary and points are on refined elements
        if ((isClose && farFromInterface) || (isCloseInterface && inInterface)) {
            pointsToRemove.insert(p);
        }
    }
    return pointsToRemove;
}

bool removePoint(Coord const& p,
                 NeighborMap& pointToNeighbors,
                 BCIDMap& pointToBcIds,
                 PointSet& remainingToRemove,
                 PointSet& alreadyRemoved,
                 std::vector<BoundaryEdge>& updatedEdges,
                 int bcId,
                 int partition) {
    const auto &neighbors = pointToNeighbors[p];

    if (neighbors.size() == 2) {
        const Coord &n1 = neighbors[0].first;
        const Coord &n2 = neighbors[1].first;

        // only remove if at least one neighbour is not on the list of removal and if both neighbors have not yet
        // been removed
        if ((!remainingToRemove.count(n1) || !remainingToRemove.count(n2)) &&
            !(alreadyRemoved.count(n1) || alreadyRemoved.count(n2))) {
            auto phi1 = pointToBcIds[n1].second;
            auto phi2 = pointToBcIds[n2].second;

            updatedEdges.emplace_back(std::make_pair(n1, n2), bcId, partition, std::make_pair(phi1, phi2));

            remainingToRemove.erase(p);
            alreadyRemoved.insert(p);

            std::cout << "REMOVED membrane point at " << p
                      << " between neighbors " << n1 << " and " << n2 << "\n"
                      << ", phi = " << 0.5*(phi2 - phi1) << "\n";
            return true;
        }
    }
    return false;
}

/**
 * Remove too close points. Idea is, that we first mark all points that need to be removed, then we go through the
 * boundaryEdges vector and find points to be removed, where maximum 1 neighbor should also be removed, then remove them
 * BUT: if one neighbor has already been removed, skip this point for the next iteration! THis way, if a chain of points
 * is marked for removal, only the outermost points will be removed. The others will successively be removed in later
 * iterations
 */
void cleanMembraneEdges(std::vector<BoundaryEdge>& boundaryEdges,
                        double& distMeanIface,
                        double& distMeanBulk) {
    // map that contains all neighboring points of each boundary point
    NeighborMap pointToNeighbors;

    // mapping from point to respective bcIDs
    BCIDMap pointToBcIds;

    // mapping from point to respective minimum refinement level of the adjacent elements
    LevelMap pointToLevels;

    // build point neighborhoods and neighboring bcIDs
    buildNeighborhoods(boundaryEdges,pointToNeighbors,pointToBcIds);

    // find locations, where the interface is
    std::vector<Coord> iFace = iFaceLocations(boundaryEdges);

    auto doRemove = Parameters::get<int>("remove interface points").value_or(0);
    if (doRemove) {
        msg("clean boundary edges if too close points exist ...");

        // mark points for removal
        auto pointsToRemove = markPointsForRemoval(pointToNeighbors,pointToBcIds,iFace,distMeanBulk,distMeanIface);

        // updated list of boundary edges
        std::vector<BoundaryEdge> updatedEdges;

        PointSet remainingToRemove = pointsToRemove;
        PointSet alreadyRemoved;

        std::vector<BoundaryEdge> currentEdges = boundaryEdges;
        auto sizeBefore = boundaryEdges.size();
        boundaryEdges.clear();

        while (!remainingToRemove.empty()) {
            updatedEdges.clear();

            // loop over current vector of all edges
            for (const auto &[edge, bcId, partition, phiValues]: currentEdges) {
                const Coord &p1 = edge.first;
                const Coord &p2 = edge.second;

                // do not add anything to the new edge vector, if one of the points on the edge have been removed already
                if (alreadyRemoved.count(p1) || alreadyRemoved.count(p2)) continue;

                // first case: nothing to remove, just add the edge to the updatedEdges
                if (remainingToRemove.count(p1) == 0 && remainingToRemove.count(p2) == 0) {
                    updatedEdges.emplace_back(edge, bcId, partition, phiValues);
                    continue;
                }

                // second case: p1 should be removed
                if (remainingToRemove.count(p1)) {
                    if (removePoint(p1,
                                    pointToNeighbors,
                                    pointToBcIds,
                                    remainingToRemove,
                                    alreadyRemoved,
                                    updatedEdges,
                                    bcId,
                                    partition))
                        continue;
                }

                // third case: p2 should be removed
                if (remainingToRemove.count(p2)) {
                    if (removePoint(p2,
                                    pointToNeighbors,
                                    pointToBcIds,
                                    remainingToRemove,
                                    alreadyRemoved,
                                    updatedEdges,
                                    bcId,
                                    partition))
                        continue;
                }

                // fallback case if none of the above cases fits: just add the facet to updatedEdges
                updatedEdges.emplace_back(edge, bcId, partition, phiValues);
            }

            // we are finished as soon as there is no point remaining to be removed
            if (remainingToRemove.empty()) {
               // msg("all points removed");
                break;
            }

            // if there are still points to remove at this point, clear maps, rebuild neighborhoods, mark points for
            // removal again and start over
            alreadyRemoved.clear();
            buildNeighborhoods(updatedEdges,pointToNeighbors,pointToBcIds);
            remainingToRemove.clear();
            remainingToRemove = markPointsForRemoval(pointToNeighbors,pointToBcIds,iFace,distMeanBulk,distMeanIface);

            currentEdges = updatedEdges;
        }

        msg("number of suface facets before removal = {}, and after removal = {}", sizeBefore, currentEdges.size());
        boundaryEdges = std::move(currentEdges);
    }
}

void addPoint(Coord const& p,
              PointMapGmsh& pointMap,
              int& nextPointTag,
              int& pTag) {
    auto it = pointMap.find(p);
    if (it == pointMap.end()) {
        pTag = nextPointTag++;
        gmsh::model::geo::addPoint(p[0], p[1], 0.0, 0.0, pTag);
        pointMap[p] = pTag;
    } else {
        pTag = it->second;
    }
}

int findStartFacet(MyGraph const& graph, int part, std::set<int>& visitedLines) {
    for (const auto& [p1, neighbors] : graph) {
        for (const auto& [p2, lineTagPartition] : neighbors) {
            int partition = lineTagPartition.second;
            auto lineTag = lineTagPartition.first;

            // only accept start point on line which has not yet been visited
            if (visitedLines.count(std::abs(lineTag)))
                continue;

            if (partition == part) {
                return p1;
            }
        }
    }
    return -1;
}

std::vector<int> buildCurveLoop(MyGraph& graph,
                    int start,
                    int part,
                    bool& surfaceInside,
                    std::set<int>& visitedLines) {
    auto current = start;
    std::vector<int> curveLoop;
    while (true) {
        bool found = false;
        for (const auto& [neighbor, lineTagPartition] : graph[current]) {
            auto lineTag = lineTagPartition.first;
            auto partition = lineTagPartition.second;

            if (partition == 1 - part || part == -1 && partition != -1) continue;

            if (visitedLines.count(std::abs(lineTag)))
                continue;

            if (partition == part || partition == -1) {
                curveLoop.push_back(lineTag);
                if (partition == -1) surfaceInside = false;
            }

            visitedLines.insert(std::abs(lineTag));
            current = neighbor;
            found = true;
            break;
        }

        if (!found || current == start)
            break;
    }
    return curveLoop;
}

std::vector<std::vector<int>> fillCurveLoops(MyGraph& graph,
               int part,
               std::set<int>& visitedLines) {
    std::vector<std::vector<int>> curveLoops;
    int start = findStartFacet(graph,part,visitedLines);
    int i = 0;
    while (start != -1) {
        bool placeHolder = true;
        curveLoops.push_back(buildCurveLoop(graph,start,part,placeHolder,visitedLines));

        if (curveLoops[i].size() < 3) {
            throw std::runtime_error("Could not build closed loop 1 from boundary edges.");
        }

        i++;
        start = findStartFacet(graph,part,visitedLines);
    }
    return curveLoops;
}

void createGmshGeometryFromBoundaryEdges(const std::vector<BoundaryEdge>& boundaryEdges,
                                         PointMapGmsh& pointMap,
                                         LineMapGmsh &lineMap,
                                         std::map<int, std::vector<int>> &physicalLines) {

    // a graph that will be filled with lines and line tags for building curve loops that are needed to
    // create (Physical) Surfaces
    MyGraph graph;

    int nextPointTag = 1;
    int nextLineTag = 1;

    // add points, lines to the gmsh::model::geo, and fill the graph for curve loops
    for (const auto& [edge, bcId, partition, phiValues] : boundaryEdges) {
        Coord p1 = edge.first;
        Coord p2 = edge.second;

        int p1Tag, p2Tag;

        // add points to gmsh
        addPoint(p1,pointMap,nextPointTag,p1Tag);
        addPoint(p2,pointMap,nextPointTag,p2Tag);

        // store point ids and reverse ids (for curve loops later)
        std::pair<int, int> key = {p1Tag, p2Tag};
        std::pair<int, int> revKey = {p2Tag, p1Tag};

        // don't add a line if it has been added already
        if (lineMap.count(key) || lineMap.count(revKey))
            continue;

        // add the line to gmsh
        gmsh::model::geo::addLine(p1Tag, p2Tag, nextLineTag);
        lineMap[key] = nextLineTag;

        // add the physical line for boundary values
        physicalLines[bcId].push_back(nextLineTag);

        //set partition of line s.t. it is 0 in surrounding, 1 in droplet, -1 on membrane
        int newPart = partition;
        if (bcId == -1) newPart = -1;

        // add the tags to the graph for building curve loops that are needed to create (Physical) Surfaces
        graph[p1Tag].emplace_back(p2Tag, std::make_pair(nextLineTag,newPart));
        graph[p2Tag].emplace_back(p1Tag, std::make_pair(-nextLineTag,newPart));

        ++nextLineTag;
    }

    // now build curve loops
    std::set<int> visitedLines;
    std::vector<int> curveLoop0;
    std::vector<std::vector<int>> curveLoop1, curveLoopMembrane;

    int start = -1;

    // find facet with partition == 0
    start = findStartFacet(graph,0,visitedLines);

    if (start == -1)
        throw std::runtime_error("could not find suitable starting point for external part of domain.");

    bool surfaceInside = true; //is surface completely inside the domain (true) or touching the boundary (true)

    // build curve loop for outer part with membrane
    curveLoop0 = buildCurveLoop(graph,start,0,surfaceInside,visitedLines);

    if (curveLoop0.size() < 3) {
        throw std::runtime_error("Could not build closed loop 0 from boundary edges.");
    }

    visitedLines.clear();

    // build curve loop for inner part with membrane
    curveLoop1 = fillCurveLoops(graph,1,visitedLines);

    // build curve loop for pure closed inner membrane
    curveLoopMembrane = fillCurveLoops(graph,-1,visitedLines);

    int curveTag1 = 2;
    std::vector<int> curveTagsInside, curveTagsAll;
    for (int i = 0; i < curveLoop1.size(); ++i) {
        gmsh::model::geo::addCurveLoop(curveLoop1[i], curveTag1);
        gmsh::model::geo::addPlaneSurface({curveTag1}, curveTag1);
        curveTagsAll.push_back(curveTag1);
        curveTag1++;
    }

    curveTagsInside.push_back(1); // 1 for external part, others for holes
    for (int i = 0; i < curveLoopMembrane.size(); ++i) {
        gmsh::model::geo::addCurveLoop(curveLoopMembrane[i], curveTag1);
        gmsh::model::geo::addPlaneSurface({curveTag1}, curveTag1);

        curveTagsAll.push_back(curveTag1);
        curveTagsInside.push_back(curveTag1);
        curveTag1++;
    }

    gmsh::model::geo::addCurveLoop(curveLoop0, 1);
    if (surfaceInside)
        gmsh::model::geo::addPlaneSurface(curveTagsInside, 1); // outer part with inner part as a hole
    else
        gmsh::model::geo::addPlaneSurface({1}, 1);

    gmsh::model::geo::addPhysicalGroup(2, {1}, 0); //external
    for (int i = 0; i < curveTagsAll.size(); ++i) {
        gmsh::model::geo::addPhysicalGroup(2, {curveTagsAll[i]}, i+1); //internal
    }
    for (const auto& [bcId, lines] : physicalLines) {
        if (bcId > -1) gmsh::model::geo::addPhysicalGroup(1, lines, bcId);
    }
}

void fillMeshSizeFields(std::vector<BoundaryEdge>& boundaryEdges,
                        PointMapGmsh const& pointMap,
                        LineMapGmsh& lineMap,
                        double distMeanIface,
                        double distMeanBulk,
                        double distMax) {
    std::vector<double> surfaceLabelsPhi1, surfaceLabelsPhi2, surfaceLabelsPhi0, surfaceLabelsDist;

    for (const auto& [edge, bcId, partition, phiValues] : boundaryEdges) {
        //if (partition != -1) continue; //only work on membrane

        const Coord& p1 = edge.first;
        const Coord& p2 = edge.second;
        double phi1 = phiValues.first;
        double phi2 = phiValues.second;

        auto it1 = pointMap.find(p1);
        auto it2 = pointMap.find(p2);
        if (it1 == pointMap.end() || it2 == pointMap.end()) continue;

        int tag1 = it1->second;
        int tag2 = it2->second;

        std::pair<int, int> key1 = {tag1, tag2};
        std::pair<int, int> key2 = {tag2, tag1};

        int lineTag = -1;
        if (lineMap.count(key1)) lineTag = lineMap[key1];
        else if (lineMap.count(key2)) lineTag = lineMap[key2];
        else continue; // Linie nicht gefunden

        // Abstand der Punkte
        double dist = two_norm(p2 - p1);

        // Mittelwert als Kriterium
        double phaseValue = 0.5 * (phi1 + phi2);

        if (phaseValue > 0.025 && phaseValue < 0.975) {
            surfaceLabelsPhi1.push_back(lineTag);
        }

        if (partition == -1) {
            if ((phaseValue < 0.025 || phaseValue > 0.975) && dist < 1.25 * distMeanIface) {
                surfaceLabelsDist.push_back(lineTag);
            }

            if (phaseValue > 0.975 && dist >= 1.25 * distMeanIface) {
                surfaceLabelsPhi2.push_back(lineTag);
            }

            if (phaseValue < 0.025 && dist >= 1.25 * distMeanIface) {
                surfaceLabelsPhi0.push_back(lineTag);
            }
        }
    }
    // parameter default values are chosen s.t. the resulting grid looks good in most cases, however, the user can
    // choose parameters for his/her own liking
    double fieldMin = Parameters::get<double>("field distMin").value_or(distMeanIface*1.5);
    double fieldMin0 = Parameters::get<double>("field distMin shell").value_or(distMeanBulk*1.5);
    double fieldMax = Parameters::get<double>("field distMax").value_or(distMax);
    double hPhase = Parameters::get<double>("mesh resolution at phase field interface").value_or(distMeanIface*1.5);
    double hPhase1 = Parameters::get<double>("mesh resolution at shell phi 1").value_or(distMeanBulk*1.5);
    double hPhase0 = Parameters::get<double>("mesh resolution at shell phi 0").value_or(distMeanBulk*1.5);
    double hOutField = Parameters::get<double>("mesh resolution at boundary field").value_or(distMax);

    gmsh::model::mesh::field::add("Distance",1);
    gmsh::model::mesh::field::setNumbers(1, "CurvesList", surfaceLabelsPhi2); //fine grid where droplet is
    gmsh::model::mesh::field::setNumber(1, "Sampling", 100);

    gmsh::model::mesh::field::add("Distance",3);
    gmsh::model::mesh::field::setNumbers(3, "CurvesList", surfaceLabelsPhi1); //fine grid where droplet is
    gmsh::model::mesh::field::setNumber(3, "Sampling", 100);

    gmsh::model::mesh::field::add("Distance",5);
    gmsh::model::mesh::field::setNumbers(5, "CurvesList", surfaceLabelsPhi0); //fine grid where droplet is
    gmsh::model::mesh::field::setNumber(5, "Sampling", 100);

    gmsh::model::mesh::field::add("Threshold", 2);
    gmsh::model::mesh::field::setNumber(2, "InField", 1);
    gmsh::model::mesh::field::setNumber(2, "SizeMin", hPhase1);
    gmsh::model::mesh::field::setNumber(2, "SizeMax", hOutField);
    gmsh::model::mesh::field::setNumber(2, "DistMin", fieldMin);
    gmsh::model::mesh::field::setNumber(2, "DistMax", fieldMax);

    gmsh::model::mesh::field::add("Threshold", 4);
    gmsh::model::mesh::field::setNumber(4, "InField", 3);
    gmsh::model::mesh::field::setNumber(4, "SizeMin", hPhase);
    gmsh::model::mesh::field::setNumber(4, "SizeMax", hOutField);
    gmsh::model::mesh::field::setNumber(4, "DistMin", fieldMin);
    gmsh::model::mesh::field::setNumber(4, "DistMax", fieldMax);

    gmsh::model::mesh::field::add("Threshold", 6);
    gmsh::model::mesh::field::setNumber(6, "InField", 5);
    gmsh::model::mesh::field::setNumber(6, "SizeMin", hPhase0);
    gmsh::model::mesh::field::setNumber(6, "SizeMax", hOutField);
    gmsh::model::mesh::field::setNumber(6, "DistMin", fieldMin0);
    gmsh::model::mesh::field::setNumber(6, "DistMax", fieldMax);

    if (surfaceLabelsDist.size()>0) {
        gmsh::model::mesh::field::add("Distance",7);
        gmsh::model::mesh::field::setNumbers(7, "CurvesList", surfaceLabelsDist); //fine grid where droplet is
        gmsh::model::mesh::field::setNumber(7, "Sampling", 100);

        gmsh::model::mesh::field::add("Threshold", 8);
        gmsh::model::mesh::field::setNumber(8, "InField", 7);
        gmsh::model::mesh::field::setNumber(8, "SizeMin", (distMeanBulk + distMeanIface)*0.25);
        gmsh::model::mesh::field::setNumber(8, "SizeMax", hOutField);
        gmsh::model::mesh::field::setNumber(8, "DistMin", fieldMin);
        gmsh::model::mesh::field::setNumber(8, "DistMax", fieldMax);
    }

    std::vector<double> fields = {2,4,6};
    if (surfaceLabelsDist.size()>0) fields.push_back(8);
    gmsh::model::mesh::field::add("Min", 700000);
    gmsh::model::mesh::field::setNumbers(700000,"FieldsList", fields);
    gmsh::model::mesh::field::setAsBackgroundMesh(700000);

    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);
}

template <class NschProb>
void generateMesh(std::string meshName,
                  std::string meshDir,
                  NschProb& nschProb,
                  AdaptInfo& adaptInfo) {
    auto gridView = nschProb.grid()->levelGridView(0);
    double time = adaptInfo.time();

    //find maximum length of all facets at the boundary
    double distMax = maximumBoundaryFacetLength(gridView);

    // extract all codim 1 facets at the boundary and the membrane
    auto boundaryEdges = extractBoundaryEdges(adaptInfo,nschProb);

    // compute 1. the mean edge length on the phase field interface,
    //         2. the mean edge length away from the interface
    double distMeanIface = 0, distMeanBulk = 0;
    meanEdgeLengths(boundaryEdges,distMeanIface,distMeanBulk);

    // and remove unnecessaryly close points on the boundaries (if wanted)
    cleanMembraneEdges(boundaryEdges,distMeanIface,distMeanBulk);

    msg("create new mesh ...");
    PointMapGmsh pointMap;
    LineMapGmsh lineMap;
    std::map<int, std::vector<int>> physicalLines;

    gmsh::initialize();

    gmsh::option::setNumber("General.Terminal", 1); //set to 0 for no output
    gmsh::model::add(meshName);

    // fill in all points, lines, surfaces and physical entities to the gmsh model
    createGmshGeometryFromBoundaryEdges(boundaryEdges,pointMap,lineMap,physicalLines);
    gmsh::model::geo::synchronize();

    // add refinement fields for mesh size information if desired
    int field = Parameters::get<int>("size by field").value_or(0);
    if (field) {
        fillMeshSizeFields(boundaryEdges,pointMap,lineMap,distMeanIface,distMeanBulk,distMax);
    }

    //do the rest of the necessary stuff and write out files
    gmsh::option::setNumber("Mesh.Smoothing",10);

    gmsh::write(meshDir + "/" + meshName + ".geo_unrolled");
    gmsh::model::mesh::generate(2);

    gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
    gmsh::write(meshDir + "/" + meshName + ".msh");

    gmsh::write(meshDir + "/backup/" + meshName + std::to_string(time) + ".geo_unrolled");
    gmsh::write(meshDir + "/backup/" + meshName + std::to_string(time) + ".msh");

    gmsh::finalize();
}

#endif //ALEMODELAMDIS2_MESHGENERATOR_H