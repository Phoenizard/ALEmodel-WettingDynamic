#pragma once

template<class NSCHProb>
class Binodal : public ExtensionInterface {
public:
    using WorldVector = typename NSCHProb::WorldVector;
    using GridView = typename NSCHProb::GridView;
    using PBF1 = decltype(surface(lagrange<2>(),std::declval<std::vector<int>>()));
    using PhiBasis = decltype(GlobalBasis{std::declval<GridView>(),std::declval<PBF1>()});
    using PBF0 = decltype(lagrange<1>());
    using Lambda2Basis = decltype(GlobalBasis{std::declval<GridView>(),std::declval<PBF0>()});
    using PBFRandoms = decltype(power<Grid::dimensionworld>(surface(lagrange<1>(),std::declval<std::vector<int>>())));
    using RandomVecBasis = decltype(GlobalBasis{std::declval<GridView>(),std::declval<PBFRandoms>()});

    enum { dow = GridView::Grid::dimensionworld };
    enum { dim = GridView::Grid::dimension };

    // Konstruktor
    Binodal(NSCHProb& nschProb_) : nschProb_(nschProb_) {
        Parameters::get("axisymmetric", axi);
        Parameters::get("parameters->eps", eps_);
        mobility_ = Parameters::get<double>("parameters->M0").value_or(1.0);
        sigma_ = 6*sqrt(2)*Parameters::get<double>("parameters->sigma").value_or(1.0);
    }

    void initData(AdaptInfo& adaptInfo) override {
        partitions_ = nschProb_.partitions();
        // compute facets vector
        std::vector<int> facets;
        facets.resize(partitions_.size(),-1);
        facets = computeFacets(*nschProb_.grid(),partitions_);

        lambda2Basis_ = std::make_shared<Lambda2Basis>(nschProb_.gridView(),lagrange<1>());
        randomBasis_ = std::make_shared<RandomVecBasis>(nschProb_.gridView(),power<Grid::dimensionworld>(surface(lagrange<1>(),facets)));
        phiBasis_ = std::make_shared<PhiBasis>(nschProb_.gridView(),surface(lagrange<2>(),facets));

        neckForce_.reset(new DOFVector<Lambda2Basis>(*lambda2Basis_));
        randoms_.reset(new DOFVector<RandomVecBasis>(*randomBasis_));
        randoms00_.reset(new DOFVector<PhiBasis>(*phiBasis_));
        randoms01_.reset(new DOFVector<PhiBasis>(*phiBasis_));
        randoms11_.reset(new DOFVector<PhiBasis>(*phiBasis_));

        invTau_ = 1.0 / adaptInfo.timestep();
    }

    void initTimestep(AdaptInfo& adaptInfo) override {
        // fill the DOFVectors for the random noise values for thermodynamic and hydrodynamic fluctuations
        nschProb_.nuPhase();
        nschProb_.rhoPhase();
        nschProb_.updateElementVectors();
        fillRandomNoisePhi();
        fillRandomNoiseV();

        invTau_ = 1.0 / adaptInfo.timestep();

        // neck force
        auto neckForce = Parameters::get<double>("binodalBC neckForce").value_or(0);
        auto radius = Parameters::get<double>("parameters->radius1_0").value_or(0);
        auto center = Parameters::get<WorldVector>("parameters->center0").value_or(WorldVector{0});

        msg("neckForce = {}, radius = {}, center = {}",neckForce, radius, center[1]);
        if (neckForce) {
            valueOf(*neckForce_) << invokeAtQP([neckForce, radius, center](FieldVector<double, dow> x) {
                if (std::abs(x[0]) < 0.5 * radius)
                    return -neckForce / std::abs(x[1] - center[1]);
                return 0.0;
            }, X());
        }


        auto binodalBC = Parameters::get<int>("binodal BCs for phi").value_or(0);
        auto binodalBCPhi = Parameters::get<double>("phase field initial value").value_or(0);
        auto reductionTime = Parameters::get<double>("time for reduction of binodalBCMu").value_or(1e15);
        if (adaptInfo.time() > reductionTime) {
            //find corner coordinates of the grid
            double xmax=-10000070;
            for (auto const& el : elements(nschProb_.gridView())) {
                for (auto const& is : intersections(nschProb_.gridView(),el)) {
                    if (is.neighbor()) //only consider boundary intersections
                        continue;
                    for (int i = 0; i<is.geometry().corners(); i++) {
                        xmax = std::max(xmax,is.geometry().corner(i)[0]);
                    }
                }
            }
            auto right = [xmax](auto const& x){ return x[0] >= xmax - 1e-7; }; // define boundary
            deltaMu_ = Parameters::get<double>("binodal deltaMu").value_or(0);

            if (!binodalBCMu_) binodalBCMu_ = integrate(nschProb_.problem().solution(_mu)*nschProb_.getPhase(),nschProb_.gridView(),4) /
                                              integrate(nschProb_.getPhase(),nschProb_.gridView(),4) + deltaMu_;
            binodalShrinkageRate_ = -Parameters::get<double>("binodal shrinkage rate").value_or(0);

            int binodalShrinkageOff =  Parameters::get<int>("no binodal shrinkage").value_or(0);
            nschProb_.dropletVolumeChange() = 0.0; //do not keep droplet volume constant any longer
            auto binodalBCMuRef = std::ref(binodalBCMu_);

            if (!binodalShrinkageRate_ && !binodalShrinkageOff) {
                nschProb_.problem().addDirichletBC(right, _phi, _mu, binodalBCMuRef);
            }
        }
    }

    void closeTimestep(AdaptInfo& adaptInfo) override {

    }

    void fillOperators(AdaptInfo& adaptInfo) override {
        auto randomNoise = Parameters::get<double>("phase field random noise").value_or(0);
        if (randomNoise > 0) {
            auto opNoise = makeOperator(tag::gradtest{}, -valueOf(*randoms_),2);
            nschProb_.problem().addVectorOperator(opNoise, _phi);
        }
        auto randomNoiseV = Parameters::get<double>("velocity random noise").value_or(0);
        if (randomNoiseV > 0) {
            auto opNoise00 = makeOperator(tag::partialtest{0}, valueOf(*randoms00_), 2);
            nschProb_.problem().addVectorOperator(opNoise00, makeTreePath(_v, 0));

            auto opNoise01 = makeOperator(tag::partialtest{0}, valueOf(*randoms01_), 2);
            nschProb_.problem().addVectorOperator(opNoise01, makeTreePath(_v, 1));

            auto opNoise10 = makeOperator(tag::partialtest{1}, valueOf(*randoms01_), 2);
            nschProb_.problem().addVectorOperator(opNoise10, makeTreePath(_v, 0));

            auto opNoise11 = makeOperator(tag::partialtest{1}, valueOf(*randoms11_), 2);
            nschProb_.problem().addVectorOperator(opNoise11, makeTreePath(_v, 1));
        }

        //binodal droplet shrinkage
        auto shrink = std::ref(binodalShrinkageRate_);
        auto deltaMu = std::ref(deltaMu_);

        auto vCons = zot(shrink * ((get(gradientOf(getPhase()), 0) * get(gradientOf(getPhase()), 0)
                                    + get(gradientOf(getPhase()), 1) * get(gradientOf(getPhase()), 1)) * deltaMu), 2);
        nschProb_.problem().addVectorOperator(vCons, _phi);

        auto vCons2 = zot(shrink * (-((get(gradientOf(getPhase()), 0) * get(gradientOf(getPhase()), 0)
                                    + get(gradientOf(getPhase()), 1) * get(gradientOf(getPhase()), 1)))), 2);
        nschProb_.problem().addMatrixOperator(vCons2, _phi,_mu);

        auto fNeck = makeOperator(tag::surface_testvec_trialvec{},valueOf(*neckForce_), 4);
        nschProb_.problem().addMatrixOperator(fNeck, makeTreePath(_fShell), makeTreePath(_kappaVec));
    }

    void fillBoundaryConditions(AdaptInfo &adaptInfo) override {
        //find corner coordinates of the grid
        double xmin = 100000, ymin = 100000, xmax = -100000, ymax = -100000;
        for (auto const &el: elements(nschProb_.gridView())) {
            for (auto const &is: intersections(nschProb_.gridView(), el)) {
                if (is.neighbor()) //only consider boundary intersections
                    continue;
                for (int i = 0; i < is.geometry().corners(); i++) {
                    xmin = std::min(xmin, is.geometry().corner(i)[0]);
                    ymin = std::min(ymin, is.geometry().corner(i)[1]);
                    xmax = std::max(xmax, is.geometry().corner(i)[0]);
                    ymax = std::max(ymax, is.geometry().corner(i)[1]);
                }
            }
        }
        FieldVector<double, 2> zero(0);

        auto bottom = [ymin](auto const &x) { return x[1] <= ymin + 1e-7; }; // define boundary
        auto left = [xmin](auto const &x) { return x[0] <= xmin + 1e-7; }; // define boundary
        auto right = [xmax](auto const &x) { return x[0] >= xmax - 1e-7; }; // define boundary
        auto top = [ymax](auto const &x) { return x[1] >= ymax - 1e-7; }; // define boundary

        // set phi on boundaries for binodal simulations
        auto binodalBC = Parameters::get<int>("binodal BCs for phi").value_or(0);
        auto binodalBCPhi = Parameters::get<double>("phase field initial value").value_or(0);
        auto neumann = Parameters::get<double>("Neumann BC value for mu binodal").value_or(0);

        if (binodalBC) {
            nschProb_.problem().addDirichletBC(top, _mu, _phi, 0.0);
            nschProb_.problem().addDirichletBC(bottom, _mu, _phi, 0.0);
            nschProb_.problem().addDirichletBC(right, _mu, _phi, binodalBCPhi);
            binodalBCMu_ = sigma_/(2.0*eps_)*(binodalBCPhi*(1.0 - binodalBCPhi)*(1.0-2.0*binodalBCPhi));
            auto bBC = std::ref(binodalBCMu_);
            if (!neumann) {
                nschProb_.problem().addDirichletBC(right, _phi, _mu, bBC);
            } else {
                auto angle = makeOperator(tag::test{}, neumann, 2);
                nschProb_.problem().addVectorOperator(3, angle, _phi);
            }
        }

    }

    void writeFiles() override {
        this->outputBulkFields_["randoms"] =
                std::make_shared<OutputField<TYPEOF(*randoms_), Grid::dimensionworld>>(randoms_);
    }

private:
    void fillRandomNoisePhi() {
        auto randomNoise = Parameters::get<double>("phase field random noise").value_or(0.0);
        randomNoise *= sqrt(mobility_);
        if (randomNoise > 0) {
            // compute the average volume of the phase field outside
            fillRandomDOFVec(valueOf(*randoms_,_0),randomNoise,0,1);
            fillRandomDOFVec(valueOf(*randoms_,_1),randomNoise,0,1);
        }
    }

    void fillRandomNoiseV() {
        auto randomNoise = Parameters::get<double>("velocity random noise").value_or(0.0);
        auto viscosity = Parameters::get<double>("stokes->viscosityOutside").value_or(0.0);
        randomNoise *= sqrt(viscosity);

        if (randomNoise > 0) {
            fillRandomDOFVec(valueOf(*randoms00_),randomNoise,1);
            fillRandomDOFVec(valueOf(*randoms01_),randomNoise);
            fillRandomDOFVec(valueOf(*randoms11_),randomNoise, 1);
        }
    }

    /// fill randomNoise vectors for thermodynamic and hydrodynamic fluctuations
    template <class DF>
    void fillRandomDOFVec(DF&& input, double noise, int diag = 0, int phiMean = 0) {
        double mean = Parameters::get<double>("normal distribution mean").value_or(0.0);
        double dev = Parameters::get<double>("normal distribution standard dev").value_or(0.2);
        double nuIn = Parameters::get<double>("stokes->viscosityInside").value_or(0.0);
        double L = Parameters::get<double>("L").value_or(0.0);
        double T = Parameters::get<double>("T").value_or(0.0);
        double kbT = 4.e-21/std::pow(L,5)*std::pow(T,2);

        input << invokeAtQP([this,&noise,&mean,&dev,&diag,nuIn,kbT,L,phiMean](FieldVector<double,2> x, double phi, double size, double visc)
                            {
                                double radius = Parameters::get<double>("initial radius").value_or(0.15);
                                double radiusIn = Parameters::get<double>("initial radius inside").value_or(0.15);
                                double use = Parameters::get<int>("use phase dependent noise").value_or(0);
                                dev = (diag ? 2.0 : sqrt(2.0)) * sqrt(kbT*invTau_/(size*eps_)) * noise;
                                std::normal_distribution<double> dist(mean, dev); // Gaussian Normal distribution with mean 0 und standard deviation 0.2
                                return ((std::pow(x[0], 2) + std::pow(x[1], 2) <= std::pow(radius, 2) &&
                                         std::pow(x[0], 2) + std::pow(x[1], 2) >= std::pow(radiusIn, 2)) ?
                                        dist(getGenerator()) : 0.0) * (use && phiMean ? 2.0*phi*(1.0-phi) : 1.0) * ((visc > nuIn) ? 1.0 : 0.0);
                            }, X(), clamp(nschProb_.getPhase(),0.0,1.0), valueOf(nschProb_.getElementSizes()), valueOf(nschProb_.getViscosityDOF()));
        makeBoundariesZero(input);
    }

    template <class DF>
    void makeBoundariesZero(DF&& input) {
        auto localView = input.basis().localView();
        auto const& gridView = nschProb_.gridView();

        auto const &level0GridView = gridView.grid().levelGridView(0);
        auto const &indexSet = level0GridView.indexSet();

        int surfaceElementIndex = 0;
        // for all elements, find surface intersections (using original grid)
        for (auto const &e: elements(gridView)) {
            auto father = e;
            while (father.hasFather())
                father = father.father();
            int p = partitions_[indexSet.index(father)];
            int q = -1;
            localView.bind(e);
            for (auto const &is: intersections(gridView, e)) {
                if (is.neighbor()) {
                    auto fatherOutside = is.outside();
                    while (fatherOutside.hasFather())
                        fatherOutside = fatherOutside.father();
                    q = partitions_[indexSet.index(fatherOutside)];
                }
                if (p == 0 && q == 1 || p == 1 && q == 0 || is.boundary()) { // now we are on a surface intersection
                    int subEntityCodim = 1;
                    int subEntityIndex = is.indexInInside(); // Local index of codim 1 entity in the inside() entity where intersection is contained in
                    auto &&subTree = localView.tree();
                    Traversal::forEachLeafNode(subTree, [&](auto const &node, auto const &tp) {
                        auto const &lfe = node.finiteElement();
                        auto const &lc = lfe.localCoefficients();
                        auto re = Dune::referenceElement(e);

                        for (std::size_t i = 0, j = 0; i < lc.size() && j < dow; ++i) {
                            auto localKey = lc.localKey(i);
                            if (re.subEntities(subEntityIndex, subEntityCodim, localKey.codim()).contains(
                                    localKey.subEntity())) {
                                auto bulk_idx = localView.index(node.localIndex(i));
                                input.coefficients().set(bulk_idx, 0.0);
                            }
                        }
                    });
                }
            }
        }
        input.coefficients().finish();
    }

    auto const& getPhase() {
        return nschProb_.getPhase();
    }

    std::mt19937& getGenerator() {
        thread_local std::mt19937 gen(std::random_device{}());
        return gen;
    }

private:
    NSCHProb& nschProb_;
    std::shared_ptr<DOFVector<Lambda2Basis>> neckForce_;        // DOFVector for neck force computation
    std::shared_ptr<DOFVector<RandomVecBasis>> randoms_;        // DOFVector for random values
    std::shared_ptr<DOFVector<PhiBasis>> randoms00_, randoms01_, randoms11_; // DOFVectors for random values

    std::shared_ptr<RandomVecBasis> randomBasis_;               // basis for randoms vector
    std::shared_ptr<PhiBasis> phiBasis_;                        // basis for randoms vector
    std::shared_ptr<Lambda2Basis> lambda2Basis_;                // basis for lambda2
    unsigned int axi;                                           // boolean if axisymmetry is used
    double mobility_;                                           // Cahn-Hilliard mobility (in bulk and on surface)
    double sigma_;                                              // Cahn-Hilliard surface tension
    double eps_;                                                // epsilon of Cahn-Hilliard equation
    double invTau_;
    double binodalBCMu_, binodalShrinkageRate_=0, deltaMu_ = 0; // for binodal shrinkage
    std::vector<int> partitions_;                               // information gathered from the mesh: IDs of boundary
};

// registry of the extension for usage in the code
using NSCHProb = NavierStokesCahnHilliard<CurvedGrid,PBF>;
REGISTER_EXTENSION(NSCHProb, "binodal", Binodal);
