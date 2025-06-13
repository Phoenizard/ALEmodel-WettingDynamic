//
// Created by mmokbel on 21.08.24.
//

#ifndef ALEMODELAMDIS2_JUMPACROSSGAMMA_HPP
#define ALEMODELAMDIS2_JUMPACROSSGAMMA_HPP

/// compute the jump of p and phi across \Gamma (needed for permeability)
template<class Prob, class GF, class DF>
void jump(Prob& prob,
          GF&& gridFunction,
          DF& result,
          int cut = 0,
          double V0 = 0) {
    auto const& gridView = prob.gridView();
    auto const& partitions = prob.partitions();
    auto const& spb = prob.globalBasis()->preBasis().subPreBasis(_phiGamma);

    result.resizeZero();
    auto shellGridView = spb.surfaceGridView();
    auto const &indexSet = gridView.grid().levelGridView(0).indexSet();
    auto localViewP = gridFunction.basis().localView(); //local view for pressure/phi basis
    auto localViewV = result.basis().localView(); //local view for continuous across \Gamma basis
    std::vector<unsigned int> idxes;

    std::vector<FieldVector<double,2>> interface;
    //first, compute coordinates, where the interface is on \Gamma
    //REMARK: does not work for inside = outside = 1 !!!
    if (cut) {
        Dune::DiscreteGridViewFunction<decltype(shellGridView), 1> phiDGVF{shellGridView, 1};
        phiSurface(phiDGVF,prob);

        auto localViewS = phiDGVF.basis().localView(); //local view for continuous across \Gamma basis
        for (auto const &e: elements(shellGridView)) {
            localViewS.bind(e);

            auto const &node = localViewS.tree();
            auto idx0 = localViewS.index(node.localIndex(0))[0];
            auto idx1 = localViewS.index(node.localIndex(1))[0];
            auto v0 = phiDGVF.coefficients()[idx0];
            auto v1 = phiDGVF.coefficients()[idx1];
            if (v0 > 0.5 && v1 <= 0.5 || v0 <= 0.5 && v1 > 0.5) {
                interface.push_back(e.geometry().center());
            }
        }
    }

    // compute pressure difference in spherical state
    auto use = Parameters::get<int>("use p0 for permeability").value_or(0);
    double p0 = 0.0;
    if (!cut && use) {
        auto sigma0 = Parameters::get<double>("parameters->sigma0").value_or(0.0);
        auto Kb = Parameters::get<double>("bendingStiffness").value_or(0.0);
        auto axi = Parameters::get<int>("axisymmetric").value_or(0);
        double rCirc = std::sqrt(V0/M_PI);
        auto dow = Prob::GridView::dimension;
        if (dow == 2 && !axi)
            p0 = sigma0 / rCirc - 0.5 * Kb / std::pow(rCirc, 3);
        else if (axi || dow == 3)
            p0 = 2.0 * sigma0 / rCirc;
    }

    //now compute jump
    for (auto const &e: elements(gridView)) {
        localViewP.bind(e);
        localViewV.bind(e);

        for (auto const &is: intersections(gridView,e)) {
            if (!is.neighbor())
                continue;

            auto father = e;
            while (father.hasFather())
                father = father.father();

            auto fatherOut = is.outside();
            while (fatherOut.hasFather())
                fatherOut = fatherOut.father();

            auto p = partitions[indexSet.index(father)];
            auto q = partitions[indexSet.index(fatherOut)];

            auto const &nodeP = localViewP.tree();
            auto const &nodeV = localViewV.tree();
            if (p == 1 && q == 0) {
                // result += gridFunction...
                for (int i = 0; i < e.geometry().corners(); ++i) { //assumes that local index == corner index...
                    unsigned int bulk_idx = localViewV.index(nodeV.localIndex(i))[0];
                    unsigned int bulk_idxP = localViewP.index(nodeP.localIndex(i))[0];
                    //check if bulk_idx was visited before, if yes, do nothing
                    if (!(std::find(idxes.begin(), idxes.end(), bulk_idxP) != idxes.end())) {
                        idxes.push_back(bulk_idxP);
                        auto idx = localViewP.index(nodeP.localIndex(i));
                        auto idxV = localViewV.index(nodeV.localIndex(i));
                        if (!cut) {
                            result.coefficients().add(idxV,p0 + gridFunction.coefficients().get(idx));
                        } else {
                            auto inside = Parameters::get<int>("phase field inside").value_or(0);
                            auto threshold = Parameters::get<double>("phi jump threshold").value_or(0);
                            if (inside) {
                                auto tru = 1;
                                for (int j = 0; j < interface.size(); ++j) {
                                    tru *= ((e.geometry().corner(i) - interface[j]).two_norm() >
                                            threshold);
                                }
                                result.coefficients().add(idxV,tru * (gridFunction.coefficients().get(idx) < 0.5));
                            }
                        }

                    }
                }

                // result -= gridFunction... on outside element
                localViewP.bind(is.outside());
                localViewV.bind(is.outside());
                for (int i = 0; i < is.outside().geometry().corners(); ++i) { //assumes that local index == corner index...
                    unsigned int bulk_idx = localViewV.index(nodeV.localIndex(i))[0];
                    unsigned int bulk_idxP = localViewP.index(nodeP.localIndex(i))[0];
                    //check if bulk_idx was visited before, if yes, do nothing
                    if (!(std::find(idxes.begin(), idxes.end(), bulk_idxP) != idxes.end())) {
                        idxes.push_back(bulk_idxP);
                        auto idx = localViewP.index(nodeP.localIndex(i));
                        auto idxV = localViewV.index(nodeV.localIndex(i));
                        if (!cut) {
                            result.coefficients().add(idxV, -gridFunction.coefficients().get(idx));
                        } else {
                            auto outside = Parameters::get<int>("phase field outside").value_or(0);
                            auto threshold = Parameters::get<double>("phi jump threshold").value_or(0);
                            if (outside) {
                                auto tru = 1;
                                for (int j = 0; j < interface.size(); ++j) {
                                    tru *= ((e.geometry().corner(i) - interface[j]).two_norm() >
                                            threshold);
                                }
                                result.coefficients().add(idxV, tru * (gridFunction.coefficients().get(idx) < 0.5));
                            }
                        }
                    }
                }
            }
        }
        localViewV.unbind();
        localViewP.unbind();
    }
    result.finish();
}

#endif //ALEMODELAMDIS2_JUMPACROSSGAMMA_HPP
