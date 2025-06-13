//
// Created by mmokbel on 21.08.24.
//

#ifndef ALEMODELAMDIS2_CORRECTINTERFACEMEMBRANECOORDS_HPP
#define ALEMODELAMDIS2_CORRECTINTERFACEMEMBRANECOORDS_HPP

template <class GEO>
std::list<int> sortShell(GEO &shellGeo) {
    std::list<int> indicesSorted;
    auto shellCoords = shellGeo.coordinates().coefficients();
    auto nIdx = shellGeo.neighborIndices();
    int startIdx = Parameters::get<int>("sortShell start index").value_or(0);
    int idx = startIdx;
    int idxRight = 0;
    int axi = Parameters::get<int>("axisymmetric").value_or(0);
    int line = Parameters::get<int>("line mesh").value_or(0);
    double ymin = 100000;
    // axisymmetric - find starting point: right point touching symmetry axis
    if (axi) {
        std::vector<FieldVector<double,2>> pointsOnSyxmmetryAxis;
        std::vector<int> idxOnSymmetryAxis;
        for (int i = 0; i < shellCoords.size(); i++) {
            if (nIdx[i].size() < 2) { //if there is only one neighbor ->shell points touching symmetry axis
                pointsOnSyxmmetryAxis.push_back(shellCoords[i]);
                idxOnSymmetryAxis.push_back(i);
            }
        }
        if (pointsOnSyxmmetryAxis[0][0] < pointsOnSyxmmetryAxis[1][0]) {
            idxRight = idxOnSymmetryAxis[1];
        } else {
            idxRight = idxOnSymmetryAxis[0];
        }
    }
    if (line) {
        for (int i = 0; i < shellCoords.size(); i++) {
            if (nIdx[i].size() < 2) { //if there is only one neighbor ->shell points touching boundary of mesh
                if (shellCoords[i][1] < ymin) { //find lower point
                    idxRight = i;
                    ymin = shellCoords[i][1];
                }
            }
        }
    }

    if (axi || line) idx = idxRight;
    indicesSorted.push_back(idx);
    for (int i = 0; i < shellCoords.size() - 1; i++) {
        idx = nIdx[idx][0]; //right neighbor index of i
        indicesSorted.push_back(idx);
    }
    return indicesSorted;
}

/// if membrane points on the interface have moved tangentially, correct them
template <class DGVF, class Geo, class DF>
void correctInterfaceMembraneCoords(DGVF& result, Geo& shellGeo, DF const& phi) {
    int count = 0;
    std::list<int> indices = sortShell(shellGeo);
    Dune::DiscreteGridViewFunction coords{shellGeo.gridView(),1};
    auto x = coords.coefficients();
    x = shellGeo.coordinates().coefficients();
    auto ind = indices.begin();
    for (int i = 0; i < indices.size(); i++) {
        Dune::FieldVector<double,2> x0 = x[*std::prev(ind)];
        Dune::FieldVector<double,2> x1 = x[*ind];
        Dune::FieldVector<double,2> x2 = x[*std::next(ind)];
        auto dist1 = two_norm(x2 - x1);
        auto dist2 = two_norm(x1 - x0);

        // if the distances are s.t. one of them is less than 45% of the other, move points
        if (dist1/dist2 < 0.85 || dist2/dist1 < 0.85) {
            auto phaseValue = phi.coefficients()[*ind];
            if (phaseValue >= 0.1 && phaseValue <= 0.9) { // only move interface points
                auto midpoint = 0.5 * (x0 + x2); // the middle of the line between x0 and x2
                auto m = two_norm(
                        midpoint - x0); // the length of the line from x0 to midpoint (or equally x1 to midpoint)
                auto dist = 0.5 * (dist1 + dist2); // the new dist1 and dist2 (both should be equal to dist)
                auto h = std::sqrt(dist * dist - m * m); // the height of the triangle (x0,M,x_new)
                auto vec = midpoint - x0;
                auto vecT = vec;
                vecT[1] = -vec[0]/two_norm(vec);
                vecT[0] = vec[1]/two_norm(vec);
                auto x_new = h * vecT + midpoint;
                auto x_new2 = -h * vecT + midpoint;
                if (two_norm(x_new - x1) < two_norm(x_new2 - x1)) {
                    x[*ind] = x_new;
                } else {
                    x[*ind] = x_new2;
                }
                std::cout << "correct coordinate x = " << x1 << " to x = " << x[*ind] << "\n";
                count++;
            }
        }
        ++ind;
        ++ind; // twice to move every second point
    }
    std::cout << "number of corrected coordinates: " << count << "\n";
    result.coefficients() = x;
}

template <class Prob, class SG, class DV, class DV2, class DF>
void correctSurfaceInterfacePointCoords(Prob& nschProb,
                                        SG& shellGeometry,
                                        DV& coords,
                                        DV2& coordinates,
                                        DF& u,
                                        double tau) {
    // define coords and phi_S as DiscreteGridViewFunction
    auto coordsDiscrete = vertexCoordinates(nschProb.gridView());
    auto const& surfaceGrid = nschProb.globalBasis()->preBasis().subPreBasis(_phiGamma).surfaceGrid();

    Dune::DiscreteGridViewFunction<decltype(surfaceGrid->leafGridView()), 1> phiDGVF{surfaceGrid->leafGridView(),1};
    phiSurface(phiDGVF, nschProb);

    // compute correction on the shel where the droplet interface is located
    correctInterfaceMembraneCoords(shellGeometry.coordinates(), shellGeometry, phiDGVF);

    // compute corrected coordinates and set displacement u accordingly
    auto displShell = coordsDiscrete;
    shellGeometry.interpolateToBulk(coordsDiscrete.coefficients(),
                                    shellGeometry.coordinates().coefficients());
    valueOf(coords) << coordsDiscrete;
    interpolateHostGrid(valueOf(*coordinates), valueOf(coords));
    u += coordsDiscrete;
    u -= displShell;
    u << u / tau;
}



#endif //ALEMODELAMDIS2_CORRECTINTERFACEMEMBRANECOORDS_HPP
