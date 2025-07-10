/**
 * Implementation of the Navier Stokes problem
 * Therefore, we use the AMDiS::BaseProblem
 * You'll need the amdis-extensions module for that
 **/

#pragma once

#include <amdis/Integrate.hpp>

namespace AMDiS {

    template <class HostGrid, class PBF>
    NavierStokesCahnHilliard<HostGrid, PBF>::NavierStokesCahnHilliard(std::string const& name
                                              , HostGrid& grid
                                              , PBF& pbf
                                              , std::vector<int> boundaryIDs
                                              , std::vector<int> partitions)
        : ProblemIterationInterface()
        , ProblemTimeInterface()
        , ProblemInstatBase(name)
        , problem_(name + "->space", grid, pbf)
        , name_(name)
        , boundaryIDs_(boundaryIDs)
        , partitions_(partitions)
        , elementSizes(grid,-1)
        , elementLevels(grid,-1)
        , distances_(grid,0.0)
    {
        // set the correct boundary IDs to the mesh
        problem().boundaryManager()->setBoundaryIds(boundaryIDs_);
        // read all needed parameters
        Parameters::get("stokes->density", rho_);
        Parameters::get("stokes->density2", rho2_);
        Parameters::get("stokes->densityOutside", rhoOut_);
        Parameters::get("stokes->viscosity", nu_);
        Parameters::get("axisymmetric", axi);
        Parameters::get("parameters->eps", eps_);
        Parameters::get("parameters->sigma", sigma_);
        Parameters::get("bendingStiffness",B_);
        Parameters::get("areaDilation", Ka_);
        Parameters::get("areaShear", Ks_);
        Parameters::get("spontaneous curvature", kappa0_);
        sigma_ = 6.0*sqrt(2)*sigma_; //rescale physically correct sigma_ to model sigma
        mobility_ = Parameters::get<double>("parameters->M0").value_or(1.0);
        permeability_ = Parameters::get<double>("parameters->permeability").value_or(0.0);
    }

    template <class HostGrid, class PBF>
    void NavierStokesCahnHilliard<HostGrid, PBF>::solveInitialProblem(AdaptInfo& adaptInfo) {
        msg("NavierStokesCahnHilliard::solveInitialProblem");

        double eps = Parameters::get<double>("parameters->eps").value_or(0.1);
        auto phase = [eps](auto s)
        {
            return 0.5*(1.0 - invokeAtQP(Operation::Tanh{}, evalAtQP(s)*(1.0/(4.0*std::sqrt(2.0)*eps))));
        };
        // indicator function to indicate, that the phase field only exists inside and/or outside
        auto elementVectorIn = ElementVector(*grid(),0);
        auto inside = Parameters::get<int>("phase field inside").value_or(0);
        auto elementVectorOut = ElementVector(*grid(),0);
        auto outside = Parameters::get<int>("phase field outside").value_or(0);
        for (auto const& e : elements(gridView())) {
            auto father = e;
            while (father.hasFather())
                father = father.father();

            if (partitions_[problem().grid()->levelGridView(0).indexSet().index(father)] == 1 && inside) {
                elementVectorIn.data()[gridView().indexSet().index(e)] = 1;
            }
            if (partitions_[problem().grid()->levelGridView(0).indexSet().index(father)] == 0 && outside) {
                elementVectorOut.data()[gridView().indexSet().index(e)] = 1;
            }
        }
        valueOf(*rhoDOF_,_mu) << valueOf(elementVectorIn) + valueOf(elementVectorOut);
        if (adaptInfo.time() < 1.e-6) {
            // initial refinement
            int ref_int = Parameters::get<int>("refinement->interface").value_or(10);
            int numCircles = Parameters::get<int>("parameters->number of circles").value_or(1);
            int initial = Parameters::get<int>("circles or initial value").value_or(0);
            /*
            if (initial) {
                double initialVal = Parameters::get<double>("phase field initial value").value_or(0);
                double randomVal = Parameters::get<double>("phase field initial value random").value_or(0);
                double radius = Parameters::get<double>("initial radius").value_or(0.15);
                double radiusIn = Parameters::get<double>("initial radius inside").value_or(0.15);
                auto randomDOF = [initialVal, randomVal, radius, radiusIn](WorldVector const &x) {
                    double r = std::pow(x[0], 2) + std::pow(x[1], 2);
                    if (r <= std::pow(radius, 2) && r >= std::pow(radiusIn, 2)) {
                        theta = std::atan2(x[1], x[0]);
                        if (theta >= - M_PI / 4 && theta <= M_PI / 4) {
                            return initialVal;
                        } else {
                            return 0.0;
                        }
                    } else {
                        return 0.0;
                    }

                    // return (std::pow(x[0], 2) + std::pow(x[1], 2) <= std::pow(radius, 2) &&
                    //         std::pow(x[0], 2) + std::pow(x[1], 2) >= std::pow(radiusIn, 2)) ?
                    //        initialVal + randomVal * ((double) rand() / (RAND_MAX)) : 0.0;
                };
                this->getPhase() << randomDOF;
            }
            */ 
            for (int i = 0; i <= ref_int; ++i) {
                if (!initial) {
                    this->getPhase() << 0.0;
                    for (int n = 0; n < numCircles; n++) {
                        //circles
                        double radius1 = Parameters::get<double>("parameters->radius1_" + std::to_string(n)).value_or(
                                0.15);
                        double radius2 = Parameters::get<double>("parameters->radius2_" + std::to_string(n)).value_or(
                                0.25);
                        WorldVector center{0.5, 0.5};
                        Parameters::get("parameters->center" + std::to_string(n), center);
                        int hiv = Parameters::get<int>("use square").value_or(0);

                        auto circleN = [radius1, radius2, &center, hiv](WorldVector const &x) {
                            if (hiv) {
                                if (x[1] < 0.0 || x[1] > 2.0*center[1]) {
                                    return 1e6;
                                }

                                double half_width = 0.5 * (radius2 + (radius1 - radius2) * (x[1] / (2.0*center[1])));

                                return 2.0*(-half_width + std::abs(x[0] - center[0]));
                            }
                            return (radius1 + radius2) * (std::sqrt(
                                    Math::sqr((x[0] - center[0]) / radius1) + Math::sqr((x[1] - center[1]) / radius2)) -
                                                          1.0);
                        };
                        this->getPhase() += clamp(phase(circleN), -0.6249, 0.3751); // phase field only inside the shell
                    }
                } else {
                    // Quarter Ring Shape without smoothing
                    double initialVal = Parameters::get<double>("phase field initial value").value_or(0);
                    double EnvVal = Parameters::get<double>("env initial phase field").value_or(0);
                    double radius = Parameters::get<double>("initial radius").value_or(0.15);
                    double radiusIn = Parameters::get<double>("initial radius inside").value_or(0.15);
                    this->getPhase() << EnvVal; // set the environment value
                    WorldVector center{0.0, 0.0};
                    auto QuarterRing = [initialVal, EnvVal, radius, radiusIn](WorldVector const &x) {
                        double r = std::pow(x[0], 2) + std::pow(x[1], 2);
                        double theta = std::atan2(x[1], x[0]);
                        if (!((r <= std::pow(radius, 2) && r >= std::pow(radiusIn, 2)) && 
                                (theta >= - M_PI / 8 && theta <= M_PI / 8))) {
                            return 0.0;
                        } else {
                            return initialVal - EnvVal; // Set the phase of the Quarter Ring to initialVal
                        }
                    };
                    this->getPhase() += QuarterRing;
                }
                // this->getPhase() << valueOf(*rhoDOF_, _mu) * this->getPhase();
                problem().markElements(adaptInfo);
                problem().adaptGrid(adaptInfo);
                updateSurfaceBasis();
                updateElementVectors();
            }
            minimumElementVolume *= 2.5;
            std::cout << "minimum element volume = " << minimumElementVolume << "\n";
        } else { // if the method is called after remeshing
            problem().markElements(adaptInfo);
            problem().adaptGrid(adaptInfo);
            updateSurfaceBasis();

            // this->getPhase() << valueOf(*rhoDOF_, _mu) * this->getPhase();
        }
        transferInitialSolution(adaptInfo);

        auto incompressibleShell = Parameters::get<int>("is shell incompressible").value_or(0);
        if (incompressibleShell) {
            // initialize concentration c with c = 1
            auto localView = problem().solution(_lambda2).basis().localView();
            for (const auto &el: elements(gridView())) {
                localView.bind(el);
                if (localView.size() == 0) continue;

                std::vector<double> localCoeff(localView.size(), 1.0);

                //fill rhsVector with calculated Values stored in intersectionVector
                problem().solution(_lambda2).coefficients().scatter(localView, localView.tree().child(_lambda2), localCoeff, Assigner::assign{});
            }
            this->solution(_lambda2).coefficients().finish();
        }
        dropletVolume0_ = integrate((2.0*M_PI*X(1)*axi + (1.0 - axi)) * clamp(getPhase(),-0.6249,0.3751) *
                invokeAtQP([](double phi) {return (phi > -0.1249) ? 1 : 0;}, clamp(getPhase(),-0.6249,0.3751)),gridView(),2);
        volume0_ = integrateShell(1.0, gridView(), partitions_,1,1);
    }

    template <class HostGrid, class PBF>
    void NavierStokesCahnHilliard<HostGrid, PBF>::initTimestep(AdaptInfo& adaptInfo)
    {
        Dune::Timer t;
        //*oldSolution_ = *problem().solutionVector();
        valueOf(*oldVelocity_) << problem().solution(_v);
        valueOf(*oldPhi_) << problem().solution(_phi);
        surfaceValueOf(*oldKappaVec_,partitions_) << problem().solution(_kappaVec);

        //valueOf(*oldMu_) << problem().solution(_mu);
        tau_ = adaptInfo.timestep();
        invTau_ = 1.0 / adaptInfo.timestep();

        nuPhase();
        rhoPhase();
        updateElementVectors();

        auto line = Parameters::get<double>("line mesh").value_or(0.0);
        if (!line) {
            auto center = Parameters::get<int>("parameters->keep membrane centered").value_or(0);

            auto vShell = integrateShell(1.0, gridView(), partitions_,1, 1);
            vGrid_[0] = -center * integrateShell(getOldVelocity(0), gridView(), partitions_, 2, 1)
                        / vShell;
            if (!axi) vGrid_[1] = -center * integrateShell(getOldVelocity(1), gridView(), partitions_, 2, 1)
                        / vShell;
            shellArea_ = integrateShell(1.0,gridView(),partitions_,1); // relaxation term
            volume_ = vShell;
        } else {
            vGrid_[0] = 0.0;
            vGrid_[1] = 0.0;
            shellArea_ = 1.0; // relaxation term
        }
        auto relaxation = Parameters::get<double>("surface concentration relaxation").value_or(0.0);
        auto incompressibleShell = Parameters::get<int>("is shell incompressible").value_or(0);

        if (incompressibleShell) {
            shellAreaChange_ = relaxation * (shellArea0_ - shellArea_) / shellArea0_;
            std::cout << "relative area change = " << shellAreaChange_ / (relaxation + 1e-10) << "\n";
        }
        dropletVolume_ = integrate((2.0*M_PI*X(1)*axi + (1.0 - axi)) * clamp(getOldPhase(),-0.6249, 0.) *
                                           invokeAtQP([](double phi) {return (phi > -0.1249) ? 1 : 0;}, clamp(getOldPhase(),-0.6249, 1.0)),gridView(),2);

        auto delta_D = Parameters::get<double>("droplet volume conservation parameter").template value_or(0);
        auto delta_V = Parameters::get<double>("volume restore constant").template value_or(0);
        dropVolChangeOverTime_ = dropletVolume0_ > 1.e-10 ? (dropletVolume0_-dropletVolume_)/dropletVolume0_*delta_D*invTau_ : 0.0;

        volChangeOverTime_ = (!line)*-delta_V * (volume_-volume0_)/volume0_;
        std::cout << "volume change = " << dropletVolume_-dropletVolume0_ << ", " << dropletVolume_ << ", " << dropletVolume0_ << "\n";

        AMDiS::info(2, "nschProb.initTimestep() needed {} seconds", t.elapsed());
    }

    template <class HostGrid, class PBF>
    void NavierStokesCahnHilliard<HostGrid, PBF>::closeTimestep(AdaptInfo& adaptInfo)
    {
        Dune::Timer t;
        writeFiles(adaptInfo);

        // clamp phase field
        auto clam = Parameters::get<int>("clamp phase field after timestep").value_or(0);
        if (clam) {
            getPhase() << clamp(getPhase(),-0.6249, 0.3751);
            solution(_phiS) << clamp(solution(_phiS),-0.6249, 0.3751);
        }

        // compute data for post processing
        auto every = Parameters::get<int>("every output file").value_or(10000);
        if (adaptInfo.timestepNumber()%every == 0 || adaptInfo.time() == adaptInfo.timestep()) {
            auto sigma0 = Parameters::get<double>("parameters->sigma0").value_or(0);
            auto sigma1 = Parameters::get<double>("parameters->sigma1").value_or(0);

            auto const &spb = solution(_kappaVec).basis().preBasis().subPreBasis(_kappaVec).subPreBasis();
            auto curvature = surfaceGridFct2D(valueOf(solution(_kappaVec)), spb, _kappaVec);

            double E_kin = integrate(
                    (2.0 * M_PI * X(1) * axi + (1.0 - axi)) * 0.5 * rhoPhase() * unary_dot(getOldVelocity()), gridView(),
                    4);
            double E_dissipate = integrate((2.0 * M_PI * X(1) * axi + (1.0 - axi)) * -0.5 * nuPhase()
                                           * 2.0 * (2.0 * get(gradientOf(getOldVelocity(0)), 0) *
                                                          get(gradientOf(getOldVelocity(0)), 0)
                                                    + 2.0 * get(gradientOf(getOldVelocity(1)), 0) *
                                                            get(gradientOf(getOldVelocity(0)), 1)
                                                    + 2.0 * get(gradientOf(getOldVelocity(1)), 1) *
                                                            get(gradientOf(getOldVelocity(1)), 1)
                                                    + get(gradientOf(getOldVelocity(0)), 1) *
                                                      get(gradientOf(getOldVelocity(0)), 1)
                                                    + get(gradientOf(getOldVelocity(1)), 0) *
                                                      get(gradientOf(getOldVelocity(1)), 0)), gridView(), 4);
            double E_bend = integrate((2.0 * M_PI * X(1) * axi + (1.0 - axi)) * 0.5 * B_ * unary_dot(makeGridFunction(curvature,spb.surfaceGridView())),
                                      solution(_kappaVec).basis().preBasis().subPreBasis(
                                              _kappaVec).subPreBasis().surfaceGridView(),
                                      4); // integrateShell(0.5 * B * unary_dot(valueOf(kappa)),gridView(),partitions_,-1);
            // double E_surfTen = integrateShell(
            //         ((sigma1 - sigma0) * (getOldPhase() * getOldPhase()) * (3.0 - 2.0 * getOldPhase()) + sigma0),
            //         gridView(), partitions_,4);

            // integrate the surface tension energy with NEW sigma function
            double E_surfTen = integrateShell((
                (sigma1 - sigma0) * (getOldPhase() * getOldPhase()) / 3 + 2 * (sigma1 - sigma0) * getOldPhase() / 3 + sigma0
            ), gridView(), partitions_,4);

            // Define E_ch with Phase
            
            double E_ch = integrate((2.0 * M_PI * X(1) * axi + (1.0 - axi)) * sigma_ *
                                    (eps_ * 0.5 * unary_dot(gradientOf(getOldPhase())) +
                                     0.25 / eps_ * pow<2>(getOldPhase()) * pow<2>(1.0 - getOldPhase())), gridView(), 4);
            double E_sum_old = E_sum;
            E_sum = E_kin + E_bend + E_ch + E_surfTen;
            double dtEsum = (E_sum - E_sum_old) / adaptInfo.timestep();

            double phiInt = integrateShell(getOldPhase(),gridView(), partitions_,4);

            std::string path;
            Parameters::get("output directory", path);
            std::string output_filename;
            Parameters::get("output_filename", output_filename);
            output_filename = path + "/" + output_filename;
            FILE *file;
            if (adaptInfo.time() == adaptInfo.timestep()) {
                file = fopen(output_filename.c_str(), "w");
                fprintf(file, "time[s], phi, dipoleMoment, dipoleMoment2, bindingFlux, bindingFlux2, E_kin, E_dissipate, E_bend, E_surfTen, E_ch, E_sum, dtE_sum\n");
                fclose(file);

            } else {
                file = fopen(output_filename.c_str(), "a");
                fprintf(file, "%10.6f, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6e, %10.6f, %10.6f, %10.6e, %10.6e\n", adaptInfo.time(),
                        phiInt, E_kin, E_dissipate, E_bend, E_surfTen, E_ch, E_sum, dtEsum);
                fclose(file);
            }
        }

        // update element size vector
        elementSizes.resizeZero();
        elementLevels.resizeZero();
        updateElementVectors(true);

        AMDiS::info(2, "nschProb.closeTimestep() needed {} seconds", t.elapsed());
    }


    template <class HostGrid, class PBF>
    void NavierStokesCahnHilliard<HostGrid, PBF>::initData(AdaptInfo& adaptInfo)
    {
        std::cout << "Navier Stokes initMarkers() ...\n";

        // compute facets vector
        std::vector<int> facets;
        facets.resize(partitions_.size(),-1);
        facets = computeFacets(*this->grid(),partitions_);

        // build basis functions
        vBasis_ = std::make_shared<VBasis>(gridView(),power<Grid::dimensionworld>(lagrange<2>(),flatInterleaved()));
        lambda2Basis_ = std::make_shared<Lambda2Basis>(gridView(),lagrange<1>());
        kappaVecBasis_ = std::make_shared<KappaVecBasis>(gridView(),power<Grid::dimensionworld>(lagrange<1>()));
        phiBasis_ = std::make_shared<PhiBasis>(gridView(),surface(lagrange<2>(),facets));
        nuBasis_ = std::make_shared<NuBasis>(gridView(),surface(lagrange<1>(),facets));

        // create DOFVectors
        oldVelocity_.reset(new DOFVector<VBasis>(*vBasis_));
        oldPhi_.reset(new DOFVector<PhiBasis>(*phiBasis_));
        oldMu_.reset(new DOFVector<PhiBasis>(*phiBasis_));
        oldSolution_.reset(new DOFVector<Basis>(*this->globalBasis()));
        rhoDOF_.reset(new DOFVector<Basis>(*this->globalBasis()));
        nuDOF_.reset(new DOFVector<NuBasis>(*nuBasis_));
        lambda2_.reset(new DOFVector<Lambda2Basis>(*lambda2Basis_));
        oldLambda1_.reset(new DOFVector<Lambda2Basis>(*lambda2Basis_));
        oldLambda2_.reset(new DOFVector<Lambda2Basis>(*lambda2Basis_));
        distancesDOF_.reset(new DOFVector<Lambda2Basis>(*lambda2Basis_));
        deltaP_.reset(new DOFVector<Lambda2Basis>(*lambda2Basis_));
        deltaPhi_.reset(new DOFVector<Lambda2Basis>(*lambda2Basis_));
        oldKappaVec_.reset(new DOFVector<KappaVecBasis>(*kappaVecBasis_));
        normalDOF_.reset(new DOFVector<KappaVecBasis>(*kappaVecBasis_));
        laplaceSolution_.reset(new DOFVector<KappaVecBasis>(*kappaVecBasis_));

        // set initial time info and area/volume of shell
        tau_ = adaptInfo.timestep();
        invTau_ = 1.0 / adaptInfo.timestep();
        int line = Parameters::get<int>("line mesh").value_or(0);
        if (!line) {
            shellArea_ = integrateShell(1.0, gridView(), partitions_,1);
            shellArea0_ = integrateShell(1.0, gridView(), partitions_,1);
            volume_ = integrateShell(1.0, gridView(), partitions_, 1, 1);
        } else {
            shellArea_ = 1;
            shellArea0_ = 1;
            volume_ = 1;
        }

        valueOf(*laplaceSolution_) << FieldVector<double,dow>(0.0);
        // initialize normalVec_
        ElementVector normalX(*grid(),0.0);
        ElementVector normalY(*grid(),0.0);
        normalVec_.push_back(normalX);
        normalVec_.push_back(normalY);
    }

    template <class HostGrid, class PBF>
    template <class Lambda>
    void NavierStokesCahnHilliard<HostGrid, PBF>::initMarkers(AdaptInfo& adaptInfo, Lambda const& lambda, double& vMin)
    {
        int ref_bulk = Parameters::get<int>("refinement->bulk").value_or(2);
        int ref_int  = Parameters::get<int>("refinement->interface").value_or(10);
        int doRefine  = Parameters::get<int>("refinement->refine with stretch").value_or(10);

        updateElementVectors();
        if (doRefine) {
            auto markerVolumes = GridFunctionMarker("a-lambda", problem().grid(),
                                                    invokeAtQP([ref_int, ref_bulk, &vMin]( double lambda, double size, double level) -> int {
                                                        return size > vMin ? std::max(
                                                                std::min((int) std::trunc(lambda * 10), ref_int), 0) :
                                                               (lambda > 1.e-6 ? level : ref_bulk);
                                                    }, lambda, valueOf(elementSizes), valueOf(elementLevels)));
            problem().addMarker(Dune::wrap_or_move(std::move(markerVolumes)));
        }
        // refine by element size: if element size at the interface is larger than vMin -> refine,
        // if smaller -> do nothing,
        // else -> coarsen
        auto marker = GridFunctionMarker("b-interface", problem().grid(),
                                         invokeAtQP([ref_int, ref_bulk, &vMin](double phi, double size, int level) -> int {
                                             return phi > 0.05 && phi < 0.95 && size > vMin ? ref_int :
                                                    (phi > 0.05 && phi < 0.95 && size <= vMin ? std::min(level,ref_int) : ref_bulk);
                                         }, this->getPhase(), valueOf(elementSizes), valueOf(elementLevels)));
        problem().addMarker(Dune::wrap_or_move(std::move(marker)));
    }


    /// the oneIteration method, including updates of kappaVecOld, distances, lambda2, and normal after grid adaptation
    template <class HostGrid, class PBF>
    template <class KV, class Y>
    void NavierStokesCahnHilliard<HostGrid, PBF>::oneIteration(AdaptInfo& adaptInfo,
                                                               KV& kappaVecOld,
                                                               Y const& y0) {
        // modification nschProb.oneIteration(adaptInfo), due to some data updates,
        // that are needed after grid refinement and before system assemble + solve
        {
            problem().beginIteration(adaptInfo);
            buildAndAdapt(adaptInfo);

            // recompute curvature vector since the grid might have changed (with finite elements)
            Dune::Timer t;
            updateNormal();
            auto const &spb = solution(_kappaVec).basis().preBasis().subPreBasis(_kappaVec).subPreBasis();
            auto const &surfaceGrid = spb.surfaceGrid();

            computeKappaVecOld(kappaVecOld, surfaceGrid, spb, normal());
            surfaceValueOf(*oldKappaVec(),partitions_) << valueOf(kappaVecOld);
            AMDiS::info(2, "computation of kappaVecOld needed {} seconds", t.elapsed());

            // update distances, lambda2 on bulk grid
            Dune::Timer t2;
            updateDistances(*this);
            if (axi && Ks_) updateLambda2(*this, y0);
            AMDiS::info(2, "update distances/lambdas needed {} seconds", t2.elapsed());

            // finally, assemble system matrix and solve the system
            assembleAndSolve(adaptInfo);
            problem().endIteration(adaptInfo);
        }
    }   

    template <class HostGrid, class PBF>
    void NavierStokesCahnHilliard<HostGrid, PBF>::fillOperators(AdaptInfo& adaptInfo)
    {
        /// Navier Stokes equation
        msg("NavierStokesCahnHilliard::fillOperators()");

        fillClassicOperators(adaptInfo);
        fillOptimizedOperators(adaptInfo);

        auto const &is = gridView().indexSet();
        using IS = decltype(is);
        auto const &partitionsLeaf = problem().solution(_xS).basis().preBasis().subPreBasis(
                _xS).subPreBasis().partitionsLeaf();
        auto const &facetsVec = problem().solution(_xS).basis().preBasis().subPreBasis(_xS).subPreBasis().facets();

        /// contact angle condition if no laplace solved
        auto ignoreLaplace = Parameters::get<int>("ignore laplace problem").value_or(0);
        if (ignoreLaplace) {
            auto outside = Parameters::get<int>("phase field outside").value_or(0);
            auto sigma0 = Parameters::get<double>("parameters->sigma0").value_or(0);
            auto sigma1 = Parameters::get<double>("parameters->sigma1").value_or(0);

            if (outside) {
                // auto angle = makeOperator(tag::surface_test<IS>{0, is, partitionsLeaf, facetsVec},
                //                           6.0 * getOldPhase() * (1.0 - getOldPhase()) * (sigma1 - sigma0), 2);
                // integrate the surface tension energy with NEW sigma function
                auto angle = makeOperator(tag::surface_test<IS>{0, is, partitionsLeaf, facetsVec},
                                              2.0 * (sigma1 - sigma0) * (getOldPhase() + 1.0) / 3.0, 2);
                problem_.addVectorOperator(angle, _mu);
            }
        }

        // gravitational force only for phase field
        WorldVector dir{0.5, 0.5};
        Parameters::get("direction", dir);
        auto gravity = Parameters::get<double>("gravity").value_or(0.0);
        if (gravity) {
            problem().addVectorOperator(zot(dir[0] * rhoPhase() * getOldPhase() * gravity), makeTreePath(_v, 0));
            problem().addVectorOperator(zot(dir[1] * rhoPhase() * getOldPhase() * gravity), makeTreePath(_v, 1));
        }

        auto incompressibleShell = Parameters::get<int>("is shell incompressible").value_or(0);
        // incompressibility operators (no optimized operators here)
        if (incompressibleShell) {
            auto const &is = gridView().indexSet();
            using IS = decltype(is);

            // div v_S =
            auto opDivVS = makeOperator(tag::surface_test_divtrialvec{}, 1.0, 2);
            problem_.addMatrixOperator(opDivVS, _lambda1, _vS);

            auto theta1 = Parameters::get<double>("surface inextensibility diffusion").value_or(0.0);
            // = \Delta_\Gamma lambda1
            auto opLaplL1 = makeOperator(tag::surface_gradtest_gradtrial{}, theta1, 2);
            problem_.addMatrixOperator(opLaplL1, _lambda1, _lambda1);

            // lagrange multiplier in NS eq
            auto opSurfGradSigma = makeOperator(tag::surface_divtestvec_trial{}, (X(1)*axi + (1.0-axi)*1.0),2);
            problem_.addMatrixOperator(opSurfGradSigma, _fShell, _lambda1);

            if (axi) {
                //  v_s_r / R =
                auto opDivVSAxi = makeOperator(tag::surface_test_trial{}, 1.0/X(1), 4);
                problem_.addMatrixOperator(opDivVSAxi, _lambda1, makeTreePath(_vS,1));

                auto laplL1Axi = makeOperator(tag::surface_test_partialtrial<IS>{1},theta1/X(1),4);
                problem().addMatrixOperator(laplL1Axi,_lambda1,_lambda1);

                auto surfGradSigmaAxi = makeOperator(tag::surface_test_trial{},1.0,4);
                problem().addMatrixOperator(surfGradSigmaAxi,makeTreePath(_fShell,1),_lambda1);

            }
            auto sAC = std::ref(shellAreaChange_);

            // div v_S = relax*(c-1)/c
            auto relax = makeOperator(tag::surface_test{}, sAC, 2);
            problem_.addVectorOperator(relax, _lambda1);
        }
    }

    template <class HostGrid, class PBF>
    void NavierStokesCahnHilliard<HostGrid, PBF>::fillClassicOperators(AdaptInfo& adaptInfo) {}

    template <class HostGrid, class PBF>
    void NavierStokesCahnHilliard<HostGrid, PBF>::fillOptimizedOperators(AdaptInfo& adaptInfo) {
        auto oldTermNS = Parameters::get<int>("useOldTermNS").value_or(0);
        auto oldTermCH = Parameters::get<int>("useOldTermCH").value_or(0);
        auto oldTermSCH = Parameters::get<int>("useOldTermSCH").value_or(0);
        auto oldTermFShell = Parameters::get<int>("useOldTermFShell").value_or(0);
        auto oldTermXS = Parameters::get<int>("useOldTermXS").value_or(0);
        auto oldTermLambdas = Parameters::get<int>("useOldTermLambdas").value_or(0);
        auto oldTermCoupling = Parameters::get<int>("useOldTermCoupling").value_or(0);

        // do nothing if all equations should be assembled with classic operators
        if (oldTermSCH && oldTermCH && oldTermNS && oldTermFShell & oldTermXS && oldTermLambdas && oldTermCoupling)
            return;

        auto invTau = std::ref(invTau_);
        auto mobility = std::ref(mobility_);
        auto vDiff = std::ref(dropVolChangeOverTime_);
        auto volDiff = std::ref(volChangeOverTime_);
        auto tau = std::ref(tau_);
        auto perm = std::ref(permeability_);
        using GV = decltype(problem().grid()->leafGridView());

        auto phi = getPhase();
        auto phiOld = getOldPhase();

        auto const &partitionsLeaf = problem().solution(_xS).basis().preBasis().subPreBasis(
                _xS).subPreBasis().partitionsLeaf();
        auto const &facetsVec = problem().solution(_xS).basis().preBasis().subPreBasis(_xS).subPreBasis().facets();
        auto const &is = gridView().indexSet();
        using IS = decltype(is);

        /// Navier-Stokes system
        if (!oldTermNS) {
            // is viscosity constant?
            auto nuIn = Parameters::get<double>("stokes->viscosityInside").value_or(1);
            auto nuOut = Parameters::get<double>("stokes->viscosityOutside").value_or(1);
            auto nuPhi = Parameters::get<double>("stokes->viscosityPhase").value_or(1);
            auto rho = Parameters::get<double>("stokes->density").value_or(1000);
            int constantViscosity = nuIn == nuOut && nuOut == nuPhi;

            using Nu = decltype(nuDOF_);
            using Phi = decltype(oldPhi_);
            using U = decltype(laplaceSolution_);
            // complete Navier-Stokes equation LHS for (nearly) constant density and constant viscosity
            // including axisymmetric terms and coupling terms to Cahn-Hilliard phase field (capillary stress only)
            auto ns = makeOperator(
                    tag::navier_stokes<Nu,Phi,U>{axi, rho, nuDOF_, invTau, volDiff,
                                                 constantViscosity, 1, oldPhi_, sigma_,
                                                 eps_, vGrid_, laplaceSolution_},
                    getOldVelocity(), 2);
            problem().addMatrixOperator(ns, makeTreePath(), makeTreePath());
            problem().addVectorOperator(ns, makeTreePath());
        }

        /// complete Cahn-Hilliard equation LHS for W=0.25\phi^2*(1-phi)^2
        /// and constant mobility, including axisymmetric terms and coupling term to velocity
        if (!oldTermCH) {
            using U = decltype(laplaceSolution_);
            auto ch = makeOperator(tag::cahn_hilliard<U>{axi, sigma_, eps_, mobility, invTau, vDiff, 1,
                                                         vGrid_, laplaceSolution_},
                                   phiOld, 2);
            problem().addMatrixOperator(ch, makeTreePath(), makeTreePath());
            auto chRHS = makeOperator(tag::cahn_hilliard<U>{axi, sigma_, eps_, mobility, invTau, vDiff, 1},
                                   phiOld, 4); // polynomial degree of RHS function in chem. potential is >2
            problem().addVectorOperator(chRHS, makeTreePath());
        }

        auto incompressibleShell = Parameters::get<int>("is shell incompressible").value_or(0);
        /// surface Force including bending stiffness, surface tension, stretching, forces due to phi and phiGamma,
        /// and axisymmetric terms
        if (!oldTermFShell) {
            /// contact angle condition
            auto inside = Parameters::get<int>("phase field inside").value_or(0);
            auto outside = Parameters::get<int>("phase field outside").value_or(0);
            auto sigma0 = Parameters::get<double>("parameters->sigma0").value_or(0);
            auto sigma1 = Parameters::get<double>("parameters->sigma1").value_or(0);
            auto lineTension = Parameters::get<double>("lineTension").value_or(0);
            auto st = Parameters::get<double>("surfaceTension").value_or(0.0);

            if (inside && outside) { //if both operators from in and out are added, sigmas have to be halved
                sigma0 *= 0.5;
                sigma1 *= 0.5;
            }

            //complete shell force
            using NType = decltype(normalVec_);
            using OType = decltype(oldKappaVec_);
            using PType = decltype(oldPhi_);
            auto shellForce = makeOperator(tag::surface_fShell<NType,OType,PType,IS>{axi, st, Ka_, Ks_, B_,
                                                                               sigma0, sigma1, kappa0_,normalVec_,
                                                                               oldKappaVec_, oldPhi_, inside, outside,
                                                                               incompressibleShell,lineTension},
                                       1.0, 4);
            auto ignoreLaplace = Parameters::get<int>("ignore laplace problem").value_or(0);

            if (!ignoreLaplace) {
                problem_.addMatrixOperator(shellForce, makeTreePath(), makeTreePath());
                problem_.addVectorOperator(shellForce, makeTreePath());
            } else {
                problem_.addMatrixOperator(makeOperator(tag::surface_testvec_trialvec{},1.0), _fShell, _fShell);
            }
        }


        // xS and kappaVec equation
        if (!oldTermXS) {
            auto normalMovement = Parameters::get<int>("move grid points with normal velocity").value_or(0);

            // prepare for permeability operators
            auto use = Parameters::get<int>("use p0 for permeability").value_or(0);
            auto sigma = Parameters::get<int>("parameters->sigma0").value_or(0);
            volume0_ = integrateShell(1.0, gridView(), partitions_, 1, 1);

            double rCirc = std::sqrt(volume0_ / M_PI);
            auto dow = GV::dimension;
            double p0 = 0.0;
            if (use) {
                if (dow == 2 && !axi)
                    p0 = sigma / rCirc - 0.5 * B_ / std::pow(rCirc, 3);
                else if (axi || dow == 3)
                    p0 = 2.0 * sigma / rCirc;
            }

            auto opXS = makeOperator(tag::surface_xS_kappaVec<decltype(normalVec_),IS>{axi, tau, normalMovement, perm, p0, kappa0_, normalVec_},
                                     abs(valueOf(*deltaPhi_)), 2);
            problem_.addMatrixOperator(opXS, makeTreePath(), makeTreePath());
            problem_.addVectorOperator(opXS, makeTreePath());
        }

        if (!oldTermLambdas && !incompressibleShell) {
            auto opLambdas = makeOperator(tag::surface_lambda1_lambda2<decltype(lambda2_),IS>{axi,tau,Ks_,lambda2_}, valueOf(distances_), 2);
            problem_.addMatrixOperator(opLambdas, makeTreePath(), makeTreePath());
            problem_.addVectorOperator(opLambdas, makeTreePath());
        }

        ///surface and bulk coupling operator (including v_S=v_\Gamma, phi_S=phi_Gamma (outside))
        if (!oldTermCoupling) {
            auto inside = Parameters::get<int>("phase field inside").value_or(0);
            auto outside = Parameters::get<int>("phase field outside").value_or(0);

            auto opCoupling = makeOperator(tag::surface_coupling<IS>{axi,inside,outside}, 1.0, 2);
            problem_.addMatrixOperator(opCoupling, makeTreePath(), makeTreePath());
        }
    }

    template <class HostGrid, class PBF>
    void NavierStokesCahnHilliard<HostGrid, PBF>::fillBoundaryConditions(AdaptInfo& adaptInfo) {
        msg("NavierStokesCahnHilliard::fillBoundaryConditions()");

        //find corner coordinates of the grid
        double xmin = 100000, ymin = 100000, xmax = -100000, ymax = -100000;
        for (auto const &el: elements(gridView())) {
            for (auto const &is: intersections(gridView(), el)) {
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

        // no-slip velocity
        int line = Parameters::get<int>("line mesh").value_or(0);
        int noMovementOnBoundariesLine = Parameters::get<int>("line mesh no movement on top and bottom").value_or(0);
        if (axi || line) {
            problem().addDirichletBC(bottom, makeTreePath(_v, 1), makeTreePath(_v, 1), 0.0);
        } else {
            problem().addDirichletBC(bottom, _v, _v, zero);
        }

        if (noMovementOnBoundariesLine && line) {
            problem().addDirichletBC(top, _v, _v, zero);
            problem().addDirichletBC(bottom, _v, _v, zero);
        }

        if (line) {
            problem().addDirichletBC(top, makeTreePath(_v, 1), makeTreePath(_v, 1), 0.0);
        } else {
            problem().addDirichletBC(top, _v, _v, zero);
        }

        auto zeroLeftRight = Parameters::get<int>("no movement on left and right").value_or(0);
        if (zeroLeftRight) {
            problem().addDirichletBC(left, _v, _v, zero);
            problem().addDirichletBC(right, _v, _v, zero);
        }

        auto zeroRight = Parameters::get<int>("no movement on right").value_or(0);
        if (zeroRight) {
            problem().addDirichletBC(right, _v, _v, zero);
        }

        // prescribe pressure
        auto pLeft = Parameters::get<int>("pressure zero left").value_or(0);
        auto pRight = Parameters::get<int>("pressure zero right").value_or(0);
        auto pTop = Parameters::get<int>("pressure zero top").value_or(0);
        auto pBottom = Parameters::get<int>("pressure zero bottom").value_or(0);
        auto p2 = Parameters::get<int>("pressure zero corner").value_or(1);
        if (p2 && !pLeft && !pRight && !pTop && !pBottom) {
            auto predicate = [xmin, xmax, ymin, ymax, line](auto const &x) {
                return (x[0] > xmax - 0.2 * (xmax - xmin) && x[1] > ymax - 1e-7)
                       || (line && (x[0] < xmin + 0.2 * (xmax - xmin) && x[1] > ymax - 1e-7));
            }; // define boundary
            problem().addDirichletBC(predicate, _p, _p, 0.0);
        } else if (pLeft) {
            problem().addDirichletBC(left, _p, _p, 0.0);
        } else if (pRight) {
            problem().addDirichletBC(right, _p, _p, 0.0);
        } else if (pTop) {
            problem().addDirichletBC(top, _p, _p, 0.0);
        } else if (pBottom) {
            problem().addDirichletBC(bottom, _p, _p, 0.0);
        }

        // set pressure to zero on internal part of membrane if no laplace problem is solved (stiff membrane)
        auto ignoreLaplace = Parameters::get<int>("ignore laplace problem").value_or(0);
        if (ignoreLaplace) {
            auto outside = Parameters::get<int>("phase field outside").value_or(0);
            if (outside) {
                problem().addConstraint(
                        MeshMovementBC{*problem().globalBasis(), _p, _p,
                                       partitions_,
                                       makeGridFunction(invokeAtQP([]() { return 0.0; }), gridView()),
                                       1
                        });
            }
        }

        //no curvature and shell force in y direction on symmetry axis
        int doConstr = 0;
        Parameters::get("axisymmetric dirichlet", doConstr);
        auto zeroFct = [](auto const &x) { return 0.0; }; // define boundary
        if ((axi || line) && doConstr) {
            problem().addConstraint(
                    AxisymmetricDirichletBC{*globalBasis(), makeTreePath(_fShell, 1), makeTreePath(_fShell, 1),
                                            zeroFct});
            problem().addConstraint(
                    AxisymmetricDirichletBC{*globalBasis(), makeTreePath(_kappaVec, 1), makeTreePath(_kappaVec, 1),
                                            zeroFct});
            problem().addConstraint(
                    AxisymmetricDirichletBC{*globalBasis(), makeTreePath(_vS, 1), makeTreePath(_vS, 1), zeroFct});
            //        problem().addConstraint(
            //                AxisymmetricDirichletBC{*globalBasis(), _lambda1,_lambda1, zeroFct});
            //        problem().addConstraint(
            //                AxisymmetricDirichletBC{*globalBasis(), _lambda2,_lambda2, zeroFct});
        }
        if (line && noMovementOnBoundariesLine) {
            auto zeroFct = [](auto const &x) { return FieldVector<double, dow>{0.0}; };  // define boundary

            problem().addConstraint(
                    AxisymmetricDirichletBC{*globalBasis(), _fShell, _fShell, zeroFct});
            problem().addConstraint(
                    AxisymmetricDirichletBC{*globalBasis(), _kappaVec, _kappaVec, zeroFct});
            problem().addConstraint(
                    AxisymmetricDirichletBC{*globalBasis(), _vS, _vS, zeroFct});
        }

        // set BCs for hivProject
        auto use_hiv = Parameters::get<int>("use hivProject").value_or(0);
        if (use_hiv) {
            auto phi1 = Parameters::get<int>("phi equals 1 inner").value_or(0);
            auto phi0 = Parameters::get<int>("phi equals zero sides").value_or(0);

            // phi = 1 on hole boundary
            if (phi1) {
                auto l = 0.5*Parameters::get<double>("nucleus diameter").value_or(0);
                auto phase = [l,this](auto const &x) {return 0.5*(1.0 + std::tanh((l - std::abs(x[0]))/(2.0*this->eps_)));};
                problem().addDirichletBC(44, _mu, _phi, phase);
            }
            // phi = 0 on nucleus boundary
            if (phi0) {
                auto l = 0.5*Parameters::get<double>("nucleus diameter").value_or(0);
                auto phase = [l,this](auto const &x) {return 0.5*(1.0 - std::tanh((std::abs(x[1]) - l)/(2.0*this->eps_)));};

                problem().addDirichletBC(33, _mu, _phi, phase);
            }

            problem().addDirichletBC(33, _v, _v, zero);
            problem().addDirichletBC(44, _v, _v, zero);

        }
    }

    template <class HostGrid, class PBF>
    void NavierStokesCahnHilliard<HostGrid, PBF>::updateNormal() {
        Dune::Timer t;
        normalVec_[0].resizeZero();
        normalVec_[1].resizeZero();
        normalDOF_->resizeZero();

        auto const &indexSet = gridView().indexSet();
        auto const &indexSet0 = gridView().grid().levelGridView(0).indexSet();
        auto localView = normalDOF_->basis().localView();

        for (auto const &e: elements(gridView())) {
            for (const auto &is: intersections(gridView(), e)) {
                if (!is.neighbor())
                    continue;

                auto father = e;
                while (father.hasFather())
                    father = father.father();

                auto fatherOut = is.outside();
                while (fatherOut.hasFather())
                    fatherOut = fatherOut.father();

                auto p = partitions_[indexSet0.index(father)];
                auto q = partitions_[indexSet0.index(fatherOut)];

                if (p == 0 && q == 1) {
                    normalVec_[0].data()[indexSet.index(e)] = is.centerUnitOuterNormal()[0];
                    normalVec_[1].data()[indexSet.index(e)] = is.centerUnitOuterNormal()[1];

                    localView.bind(e);
                    auto nodeX = localView.tree().child(0);
                    auto nodeY = localView.tree().child(1);
                    std::size_t size = nodeX.size();

                    for (std::size_t i = 0; i < size; ++i) {
                        // multiply with 0.5 to locally average values, except at the symmetry axis
                        auto axiFactor = (axi && (e.geometry().corner(i)[1] < 1.e-10)) ? 1.0 : 0.5;
                        normalDOF_->coefficients().add(localView.index(nodeX.localIndex(i)),
                                                       axiFactor * is.centerUnitOuterNormal()[0]);

                        if ((axi && (e.geometry().corner(i)[1] > 1.e-10)) || !axi) {
                            normalDOF_->coefficients().add(localView.index(nodeY.localIndex(i)),
                                                           axiFactor * is.centerUnitOuterNormal()[1]);
                        }
                    }
                    localView.unbind();
                } else if (p == 1 && q == 0) {
                    normalVec_[0].data()[indexSet.index(e)] = -is.centerUnitOuterNormal()[0];
                    normalVec_[1].data()[indexSet.index(e)] = -is.centerUnitOuterNormal()[1];
                }
            }
        }
        normalDOF_->finish();
        AMDiS::info(2, "updateNormal() needed {} seconds", t.elapsed());
    }

    template <class HostGrid, class PBF>
    template <class Expr, class GV>
    auto NavierStokesCahnHilliard<HostGrid, PBF>::integrateShell(Expr&& expr, GV const& gridView, std::vector<int> const& partitions, int order, int part)
    {
        auto axi = Parameters::get<int>("axisymmetric").value_or(0);
        auto gridFct = makeGridFunction(FWD(expr), gridView);
        auto localFct = localFunction(FWD(gridFct));
        using GridFct = TYPEOF(gridFct);
        using Range = typename GridFct::Range;
        Range result(0);

        auto const& lvgv = gridView.grid().levelGridView(0);

        //if part >= 0, integrate over subdomain with number "part", else integrate over \Gamma
        if (part >= 0) {
            for (auto const &element: elements(gridView, Dune::Partitions::interior)) {
                auto fatherEl = element;
                while (fatherEl.hasFather())
                    fatherEl = fatherEl.father();

                // integrate only if element index is different from zero
                if (part != partitions[lvgv.indexSet().index(fatherEl)])
                    continue;

                auto geometry = element.geometry();

                localFct.bind(element);
                auto const &quad = Dune::QuadratureRules<typename GV::ctype, GV::dimension>::rule(element.type(), order);

                for (auto const qp: quad) {
                    auto value = localFct(qp.position());
                    result += (axi * 2.0 * M_PI * geometry.global(qp.position())[1] + (1.0 - axi)) * value *
                              geometry.integrationElement(qp.position()) * qp.weight();
                }
                localFct.unbind();
            }
        } else {
            for (auto const &element: elements(gridView, Dune::Partitions::interior)) {
                auto fatherEl = element;
                while (fatherEl.hasFather())
                    fatherEl = fatherEl.father();

                for (auto const& is : intersections(gridView,element)) {
                    if (!is.neighbor())
                        continue;

                    auto fatherOutside = is.outside();
                    while (fatherOutside.hasFather())
                        fatherOutside = fatherOutside.father();

                    auto p = partitions[lvgv.indexSet().index(fatherEl)];
                    auto q = partitions[lvgv.indexSet().index(fatherOutside)];

                    if (p == 0 && q == 1) {
                        auto geometry = element.template subEntity<1>(is.indexInInside()).geometry();
                        auto localGeo = referenceElement(element).template geometry<1>(is.indexInInside());

                        localFct.bind(element);
                        auto const &quad = Dune::QuadratureRules<typename GV::ctype, GV::dimension - 1>::rule(is.type(), order);

                        for (auto const qp: quad) {
                            auto value = localFct(localGeo.global(qp.position()));
                                result += (axi * 2.0 * M_PI * geometry.global(qp.position())[1] + (1.0 - axi)) *
                                             value *
                                             geometry.integrationElement(qp.position()) * qp.weight();
                        }
                        localFct.unbind();
                    }

                }
            }
        }

        return gridView.comm().sum(result);
    }


} // end namespace AMDiS
