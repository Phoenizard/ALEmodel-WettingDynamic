#pragma once

template<class NSCHProb>
class Activity : public ExtensionInterface {
public:
    using GridView = typename NSCHProb::GridView;
    using PBF1 = decltype(surface(lagrange<2>(),std::declval<std::vector<int>>()));
    using PhiBasis = decltype(GlobalBasis{std::declval<GridView>(),std::declval<PBF1>()});
    using PBF0 = decltype(lagrange<1>());
    using Lambda2Basis = decltype(GlobalBasis{std::declval<GridView>(),std::declval<PBF0>()});

    enum { dow = GridView::Grid::dimensionworld };
    enum { dim = GridView::Grid::dimension };

    // Konstruktor
    Activity(NSCHProb& nschProb_) : nschProb_(nschProb_) {
        Parameters::get("axisymmetric", axi);
        Parameters::get("parameters->eps", eps_);
        Parameters::get("molecular bulk volume", molecularBulkVol_);
        Parameters::get("molecular surface area", molecularSurfaceArea_);
        Parameters::get("k0", k0_);
        Parameters::get("kOn", kOn_);
        Parameters::get("kOff", kOff_);
        k0_ *= eps_/std::pow(molecularBulkVol_,1.0/3.0);
        double sigma =  Parameters::get<double>("parameters->sigma").value_or(0);
        sigmaGamma_ = Parameters::get<double>("parameters->sigmaGamma").value_or(sigma); //suggestion by Susanne
        sigmaGamma_ *= 6.0*sqrt(2)*molecularBulkVol_/molecularSurfaceArea_;
        double mobility = Parameters::get<double>("parameters->M0").value_or(1.0);
        mobilityGamma_ = Parameters::get<double>("parameters->mobilityGamma").value_or(mobility);
        mobilityGamma_ *= molecularSurfaceArea_/molecularBulkVol_;
    }

    void initData(AdaptInfo& adaptInfo) override {
        partitions_ = nschProb_.partitions();
        // compute facets vector
        std::vector<int> facets;
        facets.resize(partitions_.size(),-1);
        facets = computeFacets(*nschProb_.grid(),partitions_);

        lambda2Basis_ = std::make_shared<Lambda2Basis>(nschProb_.gridView(),lagrange<1>());
        phiBasis_ = std::make_shared<PhiBasis>(nschProb_.gridView(),surface(lagrange<2>(),facets));
        oldPhi_.reset(new DOFVector<PhiBasis>(*phiBasis_));

        oldPhiGamma_.reset(new DOFVector<Lambda2Basis>(*lambda2Basis_));
        bindingFlux_.reset(new DOFVector<Lambda2Basis>(*lambda2Basis_));

        //set initial phiGamma:
        getSurfacePhase() << Parameters::get<double>("initial phiGamma").value_or(0.5);
        nschProb_.solution(_phiGamma) << Parameters::get<double>("initial phiGamma").value_or(0.5);
        invTau_ = 1.0 / adaptInfo.timestep();
    }

    void initTimestep(AdaptInfo& adaptInfo) override {
        surfaceValueOf(*oldPhiGamma_,partitions_) << nschProb_.problem().solution(_phiGamma);
        surfaceValueOf(*bindingFlux_,partitions_,0) << k0_ * nschProb_.solution(_mu) + kOn_ * nschProb_.solution(_phiS)
                                                       - k0_ * molecularSurfaceArea_/molecularBulkVol_ * nschProb_.solution(_muGamma)
                                                       - kOff_ * nschProb_.solution(_phiGamma);
        valueOf(*oldPhi_) << nschProb_.problem().solution(_phi);

        auto newS = Parameters::get<int>("use new binding flux").value_or(0);
        if (newS) {
            surfaceValueOf(*bindingFlux_, partitions_, 0)
                    << k0_ * nschProb_.solution(_mu) + kOn_ * (1.0 - nschProb_.solution(_phiGamma))
                       - k0_ * molecularSurfaceArea_ / molecularBulkVol_ * nschProb_.solution(_muGamma)
                       - kOff_ * nschProb_.solution(_phiGamma) * nschProb_.solution(_phiS);
        }
        //valueOf(*oldMu_) << problem().solution(_mu);
    }

    void closeTimestep(AdaptInfo& adaptInfo) override {
        auto const& gridView = nschProb_.gridView();

        // compute data for post processing
        auto every = Parameters::get<int>("every output file").value_or(10000);
        if (adaptInfo.timestepNumber()%every == 0 || adaptInfo.time() == adaptInfo.timestep()) {
            auto const &spb = nschProb_.solution(_kappaVec).basis().preBasis().subPreBasis(_kappaVec).subPreBasis();

            // compute dipole moment over droplet Radius
            auto R = Parameters::get<double>("parameters->radius1_0").value_or(10000);
            DOFVector muS(gridView, Dune::Surfacebasis::surfaceLagrange<1>(partitions_));
            DOFVector muS2(gridView, Dune::Surfacebasis::surfaceLagrange<1>(partitions_));
            surfaceValueOf(muS,partitions_,0) << X(1)*(k0_* nschProb_.solution(_mu) + kOn_ * nschProb_.solution(_phiS)
                                                       - k0_ * molecularSurfaceArea_/molecularBulkVol_ * nschProb_.solution(_muGamma)
                                                       - kOff_ * nschProb_.solution(_phiGamma));

            surfaceValueOf(muS2,partitions_,0) << X(1)*X(1)*(k0_* nschProb_.solution(_mu) + kOn_ * nschProb_.solution(_phiS)
                                                             - k0_ * molecularSurfaceArea_/molecularBulkVol_ * nschProb_.solution(_muGamma)
                                                             - kOff_ * nschProb_.solution(_phiGamma));
            auto newS = Parameters::get<int>("use new binding flux").value_or(0);
            if (newS) {
                surfaceValueOf(muS, partitions_, 0)
                        << X(1) * (k0_ * nschProb_.solution(_mu) + kOn_ * (1.0 - nschProb_.solution(_phiGamma))
                                   -
                                   k0_ * molecularSurfaceArea_ / molecularBulkVol_ * nschProb_.solution(_muGamma)
                                   - kOff_ * nschProb_.solution(_phiGamma) * nschProb_.solution(_phiS));

                surfaceValueOf(muS2, partitions_, 0)
                        << X(1) * X(1) * (k0_ * nschProb_.solution(_mu) + kOn_ * (1.0 - nschProb_.solution(_phiGamma))
                                          -
                                          k0_ * molecularSurfaceArea_ / molecularBulkVol_ * nschProb_.solution(_muGamma)
                                          - kOff_ * nschProb_.solution(_phiGamma) * nschProb_.solution(_phiS));
            }

            auto mus = surfaceGridFct(valueOf(muS),spb);
            auto mus2 = surfaceGridFct(valueOf(muS2),spb);
            auto dipoleMoment = 1.0/R*integrate(mus2,spb.surfaceGrid()->leafGridView(),4);
            auto bindingFlux = 2.0*M_PI*integrate(mus,spb.surfaceGrid()->leafGridView(),4);

            DOFVector muS_2(gridView, Dune::Surfacebasis::surfaceLagrange<1>(partitions_));
            DOFVector muS2_2(gridView, Dune::Surfacebasis::surfaceLagrange<1>(partitions_));
            surfaceValueOf(muS_2,partitions_,0) << X(1)*(k0_* nschProb_.solution(_mu) + kOn_ * clamp(nschProb_.solution(_phiS),-0.4632,1.6606)
                                                         - k0_ * molecularSurfaceArea_/molecularBulkVol_ * nschProb_.solution(_muGamma)
                                                         - kOff_ * nschProb_.solution(_phiGamma));

            surfaceValueOf(muS2_2,partitions_,0) << X(1)*X(1)*(k0_* nschProb_.solution(_mu) + kOn_ * clamp(nschProb_.solution(_phiS),-0.4632,1.6606)
                                                               - k0_ * molecularSurfaceArea_/molecularBulkVol_ * nschProb_.solution(_muGamma)
                                                               - kOff_ * nschProb_.solution(_phiGamma));

            if (newS) {
                surfaceValueOf(muS_2, partitions_, 0)
                        << X(1) * (k0_ * nschProb_.solution(_mu) + kOn_ * (1.0 - nschProb_.solution(_phiGamma))
                                   -
                                   k0_ * molecularSurfaceArea_ / molecularBulkVol_ * nschProb_.solution(_muGamma)
                                   - kOff_ * nschProb_.solution(_phiGamma) * clamp(nschProb_.solution(_phiS),-0.4632,1.6606));

                surfaceValueOf(muS2_2, partitions_, 0)
                        << X(1) * X(1) * (k0_ * nschProb_.solution(_mu) + kOn_ * (1.0 - nschProb_.solution(_phiGamma))
                                          -
                                          k0_ * molecularSurfaceArea_ / molecularBulkVol_ * nschProb_.solution(_muGamma)
                                          - kOff_ * nschProb_.solution(_phiGamma) * clamp(nschProb_.solution(_phiS),-0.4632,1.6606));
            }


            auto mus_2 = surfaceGridFct(valueOf(muS_2),spb);
            auto mus2_2 = surfaceGridFct(valueOf(muS2_2),spb);
            auto dipoleMoment2 = 1.0/R*integrate(mus2_2,spb.surfaceGrid()->leafGridView(),4);
            auto bindingFlux2 = 2.0*M_PI*integrate(mus_2,spb.surfaceGrid()->leafGridView(),4);

            double phiInt = nschProb_.integrateShell(nschProb_.getOldPhase(),gridView, partitions_,4);

            std::string path;
            Parameters::get("output directory", path);
            std::string output_filename;
            Parameters::get("output_filename", output_filename);
            output_filename = path + "/" + output_filename + "_activity";
            FILE *file;
            if (adaptInfo.time() == adaptInfo.timestep()) {
                file = fopen(output_filename.c_str(), "w");
                fprintf(file, "time[s], phi, dipoleMoment, dipoleMoment2, bindingFlux, bindingFlux2\n");
                fclose(file);

            } else {
                file = fopen(output_filename.c_str(), "a");
                fprintf(file, "%10.6f, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e,\n", adaptInfo.time(),
                        phiInt, dipoleMoment, dipoleMoment2, bindingFlux, bindingFlux2);
                fclose(file);
            }
        }
    }

    void fillOperators(AdaptInfo& adaptInfo) override {
        auto const &is = nschProb_.gridView().indexSet();
        using IS = decltype(is);
        auto const &partitionsLeaf = nschProb_.problem().solution(_xS).basis().preBasis().subPreBasis(
                _xS).subPreBasis().partitionsLeaf();
        auto const &facetsVec = nschProb_.problem().solution(_xS).basis().preBasis().subPreBasis(_xS).subPreBasis().facets();


        // bindingFlux explicit
        auto explicitS = Parameters::get<int>("use explicit binding flux").value_or(0);
        if (explicitS) {
            auto opBF = makeOperator(tag::surface_test{}, valueOf(*bindingFlux_), 2);
            nschProb_.problem().addVectorOperator(opBF, _phiGamma);

            auto opBF2 = makeOperator(tag::surface_test<IS>{0, is, partitionsLeaf, facetsVec},
                                      -molecularBulkVol_/molecularSurfaceArea_*valueOf(*bindingFlux_), 2);
            nschProb_.problem().addVectorOperator(opBF2, _phi);
        }

        auto oldTermSCH = Parameters::get<int>("useOldTermSCH").value_or(0);
        auto phiGamma = getSurfacePhase();

        /// complete Surface Cahn-Hilliard equation LHS for W=0.25*(phi-0.5)^4
        /// and constant mobility, including axisymmetric terms and coupling term to velocity
        if (!oldTermSCH) {
            auto explicitS = Parameters::get<int>("use explicit binding flux").value_or(0);
            double k0 = 0, kon = 0, koff = 0;
            if (!explicitS) {
                k0 = k0_;
                kon = kOn_;
                koff = kOff_;
            }

            auto invTau = std::ref(invTau_);
            using Phi = decltype(oldPhi_);
            auto chS = makeOperator(tag::surface_cahn_hilliard<Phi,IS>{axi, sigmaGamma_, eps_, mobilityGamma_, oldPhi_,invTau,
                                                                       k0, kon, koff,
                                                                       molecularBulkVol_, molecularSurfaceArea_},
                                    phiGamma, 2);
            nschProb_.problem().addMatrixOperator(chS, makeTreePath(), makeTreePath());
            auto chSRHS = makeOperator(tag::surface_cahn_hilliard<Phi,IS>{axi, sigmaGamma_, eps_, mobilityGamma_, oldPhi_,invTau,
                                                                          k0, kon, koff,
                                                                          molecularBulkVol_, molecularSurfaceArea_},
                                       phiGamma, 4);
            nschProb_.problem().addVectorOperator(chSRHS, makeTreePath());
        }

        auto ignoreLaplace = Parameters::get<int>("ignore laplace problem").value_or(0);
        auto normalVel = Parameters::get<int>("move grid points with normal velocity").value_or(0);
        if (ignoreLaplace || normalVel) {
            // material derivative convection vanishes only if grid movement with velocity is imposed
            // <grad(phi)*u^new, psi>
            auto phiGamma = getSurfacePhase();
            auto opAdvect2 = makeOperator(tag::surface_test_trialvec{}, gradientOf(phiGamma), 2);
            nschProb_.problem().addMatrixOperator(opAdvect2, _phiGamma, _vS);

            auto opAdvect2b = makeOperator(tag::surface_test_gradtrial{}, -valueOf(*nschProb_.laplaceSolution()), 2);
            nschProb_.problem().addMatrixOperator(opAdvect2b, _phiGamma, _phiGamma);
        }
    }

    void fillBoundaryConditions(AdaptInfo &adaptInfo) override {

    }

    void writeFiles() override {
        DOFVector muS2(nschProb_.gridView(), Dune::Surfacebasis::surfaceLagrange<1>(nschProb_.partitionsRef()));
        DOFVector muS3(nschProb_.gridView(), Dune::Surfacebasis::surfaceLagrange<1>(nschProb_.partitionsRef()));

        auto k0_ = Parameters::get<double>("k0").value_or(0);
        auto kOn_ = Parameters::get<double>("kOn").value_or(0);
        auto kOff_ = Parameters::get<double>("kOff").value_or(0);
        auto molecularBulkVol_ = Parameters::get<double>("molecular bulk volume").value_or(0);
        auto molecularSurfaceArea_ = Parameters::get<double>("molecular surface area").value_or(0);

        surfaceValueOf(muS2,nschProb_.partitionsRef(),0) << k0_* nschProb_.solution(_mu) + kOn_ * nschProb_.solution(_phiS)
                                                           - k0_ * molecularSurfaceArea_/molecularBulkVol_ * nschProb_.solution(_muGamma)
                                                           - kOff_ * nschProb_.solution(_phiGamma);

        surfaceValueOf(muS3,nschProb_.partitionsRef(),0) << k0_* nschProb_.solution(_mu) + kOn_ * clamp(nschProb_.solution(_phiS),-0.4632,1.6606)
                                                           - k0_ * molecularSurfaceArea_/molecularBulkVol_ * nschProb_.solution(_muGamma)
                                                           - kOff_ * nschProb_.solution(_phiGamma);

        auto newS = Parameters::get<int>("use new binding flux").value_or(0);
        if (newS) {
            surfaceValueOf(muS2, nschProb_.partitionsRef(), 0)
                    << (k0_ * nschProb_.solution(_mu) + kOn_ * (1.0 - nschProb_.solution(_phiGamma))
                        -
                        k0_ * molecularSurfaceArea_ / molecularBulkVol_ * nschProb_.solution(_muGamma)
                        - kOff_ * nschProb_.solution(_phiGamma) * nschProb_.solution(_phiS) * nschProb_.solution(_phiS));

            surfaceValueOf(muS3, nschProb_.partitionsRef(), 0)
                    <<  (k0_ * nschProb_.solution(_mu) + kOn_ * (1.0 - nschProb_.solution(_phiGamma))
                         -
                         k0_ * molecularSurfaceArea_ / molecularBulkVol_ * nschProb_.solution(_muGamma)
                         - kOff_ * nschProb_.solution(_phiGamma) * clamp(nschProb_.solution(_phiS),-0.4632,1.6606));
        }

        auto const &spb = nschProb_.globalBasis()->preBasis().subPreBasis(_phiGamma);

        auto mus2 = surfaceGridFct(valueOf(muS2), spb);
        auto mus3 = surfaceGridFct(valueOf(muS3), spb);
        auto phiGamma = surfaceGridFct(nschProb_.solution(_phiGamma), spb, _phiGamma);
        auto muGamma = surfaceGridFct(nschProb_.solution(_muGamma), spb, _muGamma);

        using SurfaceGridFct = TYPEOF(mus2);
        using Output1D = OutputSurfaceField<SurfaceGridFct, 1>;
        this->outputSurfaceFields_["bindingFlux"] =
                std::make_shared<Output1D>(std::make_shared<SurfaceGridFct>(mus2));
        this->outputSurfaceFields_["bindingFluxClamped"] =
                std::make_shared<Output1D>(std::make_shared<SurfaceGridFct>(mus3));
        this->outputSurfaceFields_["phiGamma"] =
                std::make_shared<Output1D>(std::make_shared<SurfaceGridFct>(phiGamma));
        this->outputSurfaceFields_["muGamma"] =
                std::make_shared<Output1D>(std::make_shared<SurfaceGridFct>(muGamma));
    }

    /// surface phase field component
    auto getSurfacePhase(int = 0)       { return valueOf(*oldPhiGamma_);  }
    auto getSurfacePhase(int = 0) const { return valueOf(*oldPhiGamma_);  }

private:
    NSCHProb& nschProb_;
    std::shared_ptr<DOFVector<Lambda2Basis>> bindingFlux_;      // DOFVectors for net binding flux s(phiGamma,phi)
    std::shared_ptr<DOFVector<Lambda2Basis>> oldPhiGamma_;      // DOFVector for old phiGamma
    std::shared_ptr<DOFVector<PhiBasis>> oldPhi_;               // DOFVector for phase
    std::shared_ptr<PhiBasis> phiBasis_;                        // basis for phi
    std::shared_ptr<Lambda2Basis> lambda2Basis_;                // basis for lambda2
    double molecularBulkVol_=1.0, molecularSurfaceArea_=1.0;    // surface modelExtensions parameters
    double k0_, kOn_, kOff_;                                    // surface modelExtensions parameters
    unsigned int axi;                                           // boolean if axisymmetry is used
    double sigmaGamma_;                                         // surface tension of surface CH
    double mobilityGamma_;                                      // Cahn-Hilliard mobility (in bulk and on surface)
    double eps_;                                                // epsilon of Cahn-Hilliard equation
    double invTau_;
    std::vector<int> boundaryIDs_, partitions_;                 // information gathered from the mesh: IDs of boundary
};

// registry of the extension for usage in the code
using NSCHProb = NavierStokesCahnHilliard<CurvedGrid,PBF>;
REGISTER_EXTENSION(NSCHProb, "activity", Activity);
