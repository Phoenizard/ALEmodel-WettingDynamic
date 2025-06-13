//
// Created by mmokbel on 21.08.24.
//

#ifndef ALEMODELAMDIS2_UPDATEFUNCTIONS_HPP
#define ALEMODELAMDIS2_UPDATEFUNCTIONS_HPP

/// contains all functions that are needed to update data if the grid has
/// changed due to apdatation or grid movement

/// update of distances ElementVector and DOFVector, both used for the stretching force
/// distances contains lambda1*lambda2 in axisymmetric case and lambda1 in general 2D case
template <class P>
void updateDistances(P& prob) {
    Dune::Timer t;
    auto& distances = prob.distances();
    auto& distancesDOF = prob.distancesDOF();
    auto& oldLambda1 = prob.oldLambda1();
    auto& partitions = prob.partitionsRef();
    distances.resizeZero();
    auto localView = oldLambda1->basis().localView();
    auto axi = Parameters::get<int>("axisymmetric").value_or(0);
    auto Ks = Parameters::get<double>("areaShear").value_or(0);


    auto const &indexSet = prob.gridView().indexSet();
    auto const &indexSet0 = prob.gridView().grid().levelGridView(0).indexSet();
    for (auto const &e: elements(prob.gridView())) {
        for (const auto &is: intersections(prob.gridView(), e)) {
            if (!is.neighbor())
                continue;

            auto father = e;
            while (father.hasFather())
                father = father.father();

            auto fatherOut = is.outside();
            while (fatherOut.hasFather())
                fatherOut = fatherOut.father();

            auto p = partitions[indexSet0.index(father)];
            auto q = partitions[indexSet0.index(fatherOut)];
            if (p != q) {
                localView.bind(e);
                auto const &node = localView.tree();
                int subEntityCodim = 1;
                int subEntityIndex = is.indexInInside(); // Local index of codim 1 entity in the inside() entity where intersection is contained in

                auto const &lfe = node.finiteElement();
                auto const &lc = lfe.localCoefficients();
                auto re = Dune::referenceElement(e);

                // get the values of oldLambda1 on the corners
                std::vector<double> vals;
                for (std::size_t i = 0, k = 0; i < lc.size() && k < 2; ++i) {
                    auto localKey = lc.localKey(i);
                    if (re.subEntities(subEntityIndex, subEntityCodim, localKey.codim()).contains(
                            localKey.subEntity())) {
                        vals.push_back(valueOf(*oldLambda1).coefficients().get(localView.index(node.localIndex(i))));
                    }
                }
                // compute mean value of the two values on the corners
                auto mean = (vals[1] + vals[0]) * 0.5;

                // multiply surface segment volume with mean
                distances.data()[indexSet.index(e)] = is.geometry().volume() * mean;
                if (axi && !Ks) distances.data()[indexSet.index(e)] *= is.geometry().center()[1];
            }
        }
    }

    for (auto const &e: elements(prob.grid()->hostGrid()->hostGrid().hostGrid()->leafGridView())) {
        auto const &indexSet = prob.grid()->hostGrid()->hostGrid().hostGrid()->leafGridView().indexSet();
        auto const &indexSet0 = prob.grid()->hostGrid()->hostGrid().hostGrid()->levelGridView(0).indexSet();
        for (const auto &is: intersections(prob.grid()->hostGrid()->hostGrid().hostGrid()->leafGridView(), e)) {
            if (!is.neighbor())
                continue;

            auto father = e;
            while (father.hasFather())
                father = father.father();

            auto fatherOut = is.outside();
            while (fatherOut.hasFather())
                fatherOut = fatherOut.father();

            auto p = partitions[indexSet0.index(father)];
            auto q = partitions[indexSet0.index(fatherOut)];
            if (p != q) {
                distances.data()[indexSet.index(e)] /= is.geometry().volume();
                if (axi && !Ks) distances.data()[indexSet.index(e)] /= is.geometry().center()[1];
            }
        }
    }

    // now interpolate the solution to a DOFVector
    distancesDOF->resizeZero();

    // first, interpolate everything to a DOFVector (subtracted by 1) on the surfaceLagrange Basis
    DOFVector distancesDOFS(prob.gridView(),
                            power<1>(Dune::Surfacebasis::surfaceLagrange<1>(partitions), flatInterleaved()));
    valueOf(distancesDOFS).interpolate(valueOf(prob.distances()) - 1.0, tag::average{});

    // then, interpolate to the actual DOFVector on the bulk basis (and add 1 again)
    surfaceValueOf(*distancesDOF,partitions) << valueOf(distancesDOFS,_0);
    valueOf(*distancesDOF) += 1.0;
    AMDiS::info(2, "updateDistances() needed {} seconds", t.elapsed());
}

/// update the second principal stretch lambda2 in axisymmetric case (nor used in general 2D case)
template <class P, class DV>
void updateLambda2(P& prob,
                   DV const& y0) {
    auto axi = Parameters::get<int>("axisymmetric").value_or(0);
    if (axi) {
        valueOf(*prob.lambda2()) << valueOf(*prob.oldLambda2()) * invokeAtQP([](double y, double y0) {
                                                                  return std::abs(y - y0) < 1.e-8 ? 1.0 : y / y0;
                                                              },
                                                              X(1), valueOf(y0));

        auto const &indexSet = prob.gridView().indexSet();
        auto const &indexSet0 = prob.gridView().grid().levelGridView(0).indexSet();
        auto localView = prob.lambda2()->basis().localView();

        for (auto const &e: elements(prob.gridView())) {
            for (const auto &is: intersections(prob.gridView(), e)) {
                if (!is.neighbor())
                    continue;

                auto father = e;
                while (father.hasFather())
                    father = father.father();

                auto fatherOut = is.outside();
                while (fatherOut.hasFather())
                    fatherOut = fatherOut.father();

                auto p = prob.partitions()[indexSet0.index(father)];
                auto q = prob.partitions()[indexSet0.index(fatherOut)];
                if (p != q) {
                    localView.bind(e);
                    auto const &node = localView.tree();

                    for (std::size_t j = 0; j < e.geometry().corners(); ++j) {
                        auto x = e.geometry().corner(j);
                        if (x[1] < 1.e-10) {
                            auto val = prob.distances().data()[indexSet.index(e)];
                            valueOf(*prob.lambda2()).coefficients().set(localView.index(node.localIndex(j)), val);
                        }
                    }
                }
            }
        }
        prob.lambda2()->finish();
    }
}

/// update the DOFVectors containing the jump of phi and/or p across \Gamma (needed for permeability)
template <class Prob>
void updateJumps(Prob& nschProb, double V0) {
    auto perm = Parameters::get<double>("parameters->permeability").value_or(0);

    if (perm) {
        Dune::Timer t;
        auto const &gridView = nschProb.gridView();
        auto const &surfacePreBasis = nschProb.globalBasis()->preBasis().subPreBasis(_phiGamma);
        auto const &partitions = nschProb.partitions();

        auto p = nschProb.solution(_p);
        auto dpCorrected = makeDOFVector(nschProb.nuBasis()); //DOFVector on P1 basis
        auto phi = makeDOFVector(nschProb.nuBasis()); //DOFVector on P1 basis
        valueOf(phi) << nschProb.getOldPhase();
        std::cout << "compute jump of phi\n";
        jump(nschProb, valueOf(phi), *nschProb.deltaPhi(), 1);

        auto sigma = Parameters::get<double>("parameters->sigma").value_or(0);
        sigma = 6 * sqrt(2) * sigma; //rescale physically correct sigma_ to model sigma
        auto eps = Parameters::get<double>("parameters->eps").value_or(0);
        valueOf(dpCorrected) << p - (sigma * eps * 0.5 * pow(two_norm(gradientOf(phi)), 2) +
                                     sigma / eps * 0.25 * (pow(valueOf(phi), 2) * pow(1.0 - valueOf(phi), 2)));
        std::cout << "compute jump of p\n";
        jump(nschProb, valueOf(dpCorrected), *nschProb.deltaP(), 0, V0);
        AMDiS::info(2, "updateJumps() needed {} seconds", t.elapsed());
    }
}


/// update DiscreteGridViewFunctions to new grid and set them to zero
template <class DGVF>
void updateZero(DGVF& dgvf) {
    Dune::Functions::interpolate(dgvf.basis(),
                                 dgvf.coefficients(),
                                 [](auto const &x) { return 0.0; });
}

template <class DGVF, class... Args>
void updateZero(DGVF& dgvf, Args... args) {
    Dune::Functions::interpolate(dgvf.basis(),
                                 dgvf.coefficients(),
                                 [](auto const &x) { return 0.0; });
    updateZero(args...);
}

/// update DiscreteGridViewFunctions to new grid and set them to zero
template <class GV, class DGVF>
void updateZero(GV const& gridView, DGVF& dgvf) {
    dgvf.update(gridView);
    Dune::Functions::interpolate(dgvf.basis(),
                                 dgvf.coefficients(),
                                 [](auto const &x) { return 0.0; });
}

template <class GV, class DGVF, class... Args>
void updateZero(GV const& gridView, DGVF& dgvf, Args... args) {
    dgvf.update(gridView);
    Dune::Functions::interpolate(dgvf.basis(),
                                 dgvf.coefficients(),
                                 [](auto const &x) { return 0.0; });
    updateZero(gridView, args...);
}

/// perform a manual update of the surface grid stored in each instance of the surfaceLagrangeBasis
/// the surfaceLagrangeBasis is in general updated automatically, but to save time for large number of instances of
/// the surfaceLagrangeBasis, it makes sense to only to the automatic update once and then store the computed
/// results of the surfaceGrid, surfaceToFluidMap and surfaceToFluidMapIdx in the other instances
template <int i, class P, class SG, class Map, class TP>
void updateSurfaceBasisManually(P& prob,
                                SG const& sg,
                                Map const& stf,
                                std::map<int,int> const& stfI,
                                TP tp) {
    if constexpr (i == 1) {
        prob.globalBasis()->preBasis().subPreBasis(tp).update(sg,stf,stfI);
    }
    else if constexpr (i == 2) {
        prob.globalBasis()->preBasis().subPreBasis(tp).subPreBasis().update(sg,stf,stfI);
    }
}

template <int i, class P, class SG, class Map, class TP, class... Types>
void updateSurfaceBasisManually(P& prob,
                                SG const& sg,
                                Map const& stf,
                                std::map<int,int> const& stfI,
                                TP tp,
                                Types... args) {
    if constexpr (i == 1) {
        prob.globalBasis()->preBasis().subPreBasis(tp).update(sg,stf,stfI);
    }
    else if constexpr (i == 2) {
        prob.globalBasis()->preBasis().subPreBasis(tp).subPreBasis().update(sg,stf,stfI);
    }
    updateSurfaceBasisManually<i>(prob, sg, stf, stfI, args...);
}

#endif //ALEMODELAMDIS2_UPDATEFUNCTIONS_HPP
