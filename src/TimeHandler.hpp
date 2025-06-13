//
// Created by mmokbel on 22.08.24.
//

#ifndef ALEMODELAMDIS2_TIMEHANDLER_HPP
#define ALEMODELAMDIS2_TIMEHANDLER_HPP

using Coord = FieldVector<double, 2>; // only 2D for now, 3D to be implemented
using Edge = std::pair<Coord, Coord>;

/** \brief
   * Class to handle the timesteps with a set of functions that are designed to make the code more readable
   */
template <class Grid>
class TimeHandler {
    using MainProb = NavierStokesCahnHilliard<CurvedGrid,PBF>;
    using LProb = Laplace<CurvedGrid>;
    using ExtensionManager = ExtensionManager<MainProb>;

public:
    explicit TimeHandler(double& time) // the current time
                :
                mainProb_(std::nullopt),
                laplaceProb_(std::nullopt),
                adaptInfo_(std::nullopt),
                time_(time),
                tau_(std::nullopt)
    {
        if (time > 0) noRemeshingHappenedBefore_ = 0;
    }

    /// \brief set the problems into the handler
    void setProblems(MainProb& prob, LProb& lProb, ExtensionManager& ex, AdaptInfo& aI) {
        mainProb_ = std::ref(prob);
        laplaceProb_ = std::ref(lProb);
        extensions_ = std::ref(ex);
        adaptInfo_ = std::ref(aI);
        tau_ = std::ref(aI.timestep());
    }

    /** \brief compute element sizes needed for refinement at the interface
     * @param gridView - the gridView of the bulk grid
     * @param partitions - a vector containing elementIds (1 inside the shell, 0 outside)
     * @param time - the current time of the simulation
     */
    template <class GV>
    void estimatedElementSizes(GV const& gridView,
                               std::vector<int> const& partitions,
                               double const& time) {
        int vAvgCount = 0;
        int surfOnly  = Parameters::get<int>("refinement vAvg only on surface").value_or(1);

        for (auto const& e : elements(gridView)) {
            for (auto const& is : intersections(gridView, e)) {
                if (!is.neighbor())
                    continue;

                auto p = partitions[gridView.indexSet().index(is.inside())];
                auto q = partitions[gridView.indexSet().index(is.outside())];
                if (p != q || !surfOnly) {
                    vAvg_ += e.geometry().volume();
                    vAvgCount++;
                    //if (vMin > e.geometry().volume()) vMin = e.geometry().volume();
                }
            }
        }
        vAvg_ /= vAvgCount;
        int ref_int  = Parameters::get<int>("refinement->interface").value_or(10);
        if (time == 0) {
            int startReduced = Parameters::get<int>("start with reduced time step").value_or(0);
            if (startReduced) ref_int = Parameters::get<int>("refinement->interface at start").value_or(10);
            vMin_ = vAvg_ * 4 * std::pow(0.5,ref_int);              //minimum element volume of the grid after the initial refinement
            std::cout << "minimum element size (computed): " << vMin_ << "\n";
        }
    }

    /// \brief initialize the main problem and the laplace problem
    void initializeData() {
        auto& mainProb = mainProb_->get();
        auto& laplaceProb = laplaceProb_->get();
        auto& extensions = extensions_->get();

        auto& adaptInfo = adaptInfo_->get();


        // initialize the data for the problem(s)
        mainProb.initialize(INIT_ALL);
        laplaceProb.initialize(INIT_ALL);

        // set the datatransfer to properly update data that lives on the surfaceLagrange basis
        mainProb.problem().solutionVector()->setDataTransfer(tag::interpolation_datatransfer_surface{});

        mainProb.initBaseProblem(adaptInfo);
        laplaceProb.initBaseProblem(adaptInfo);
        mainProb.initTimeInterface();
        laplaceProb.initTimeInterface();

        extensions.initData(adaptInfo);
    }

    /** \brief initialize marker for refinement and create the initial phase field
     * @param shellSize0 - the area (not volume!) of the shell in the initial state
     */
    void initializeMainProb(double& shellSize0) {
        auto& mainProb = mainProb_->get();
        auto& extensions = extensions_->get();
        auto& adaptInfo = adaptInfo_->get();


        mainProb.initMarkers(adaptInfo,
                             makeGridFunction(valueOf(mainProb.distances()) - 1.0, mainProb.gridView()),
                             vMin_);
        if (time_ == 0)
            shellSize0 = mainProb.integrateShell(1.0,mainProb.gridView(),mainProb.partitions());
        mainProb.setShellArea0(shellSize0);

        // perform initial grid refinement
        mainProb.solveInitialProblem(adaptInfo);

        extensions.fillOperators(adaptInfo);
        extensions.fillBoundaryConditions(adaptInfo);
    }

    /** \brief compute values for droplet and shell volume in initial state (once comuted - saved throughout the simulation)
     */
    void computeInitialVolumes() {
        auto& mainProb = mainProb_->get();
        if (time_ == 0) {
            dropletVolume0_ = mainProb.getDropletVolume0(); // intitial volume of the droplet (where phi = 1)
            volume0_ = mainProb.getVolume0(); // initial volume of the closed shell
        } else if (remeshingJustHappened_) {
            if (volume0_ == 0) volume0_ = mainProb.getVolume0(); // initial volume of the closed shell
            if (dropletVolume0_ == 0) dropletVolume0_ = mainProb.getDropletVolume0(); // initial volume of the droplet (where phi = 1)
            mainProb.setDropletVolume0(dropletVolume0_);
            mainProb.setVolume0(volume0_);
            std::cout << "minimum element volume = " << vMin_ << "\n";
        }
    }

    template <class G>
    std::vector<Edge> extractBoundaryIndices(G const& gridPtr) {
        auto lvgv = gridPtr->levelGridView(0);
        std::vector<Edge> boundaryIndices(gridPtr->numBoundarySegments());

        for (auto const &element: elements(lvgv)) {
            for (const auto &intersection: intersections(lvgv, element)) {
                if (intersection.boundary()) {
                    // Zugriff auf die boundaryId
                    auto idx = intersection.boundarySegmentIndex();
                    boundaryIndices[idx] = std::make_pair(intersection.geometry().corner(0),
                                                          intersection.geometry().corner(1));
                }
            }
        }
        return boundaryIndices;
    }

    template <class G>
    void restoreBoundaryIds(G& gridPtr,
                            std::vector<Edge> const& boundaryIndices,
                            std::vector<int>& boundaryIds) {
        auto lvgv = gridPtr->levelGridView(0);

        auto boundaryIdsNew = boundaryIds;
        for (auto const &element: elements(lvgv)) {
            for (const auto &intersection: intersections(gridPtr->levelGridView(0), element)) {
                if (intersection.boundary()) {
                    for (int i = 0; i < boundaryIndices.size(); i++) {
                        auto x1 = intersection.geometry().corner(0);
                        auto x2 = intersection.geometry().corner(1);
                        auto [x12,x22] = boundaryIndices[i];
                        auto equalX1 = (x1-x12).two_norm() < 1e-10 || (x1-x22).two_norm() < 1e-10;
                        auto equalX2 = (x2-x12).two_norm() < 1e-10 || (x2-x22).two_norm() < 1e-10;

                        auto id = boundaryIds[i];
                        if (equalX1 && equalX2) {
                            auto idx = intersection.boundarySegmentIndex();
                            boundaryIdsNew[idx] = boundaryIds[i];
                        }
                    }
                }
            }
        }
        boundaryIds = std::move(boundaryIdsNew);
    }

    /** an update function, that updates the grid to the new one, if remeshing has been done in the previous time step
     *
     * @param gridPtr - the unterlying grid to be updated
     * @param partitions - a vector containing elementIds (1 inside the shell, 0 outside)
     * @param facets  - a vector containing facet indices of the facet touching the surface (-1 if element not on the surface)
     * @param boundaryIds - a vector containing the ID of the boundary intersections of all elements that touch a boundary
     * @param count - the number of the current time step
     * @param time - the current time of the simulation
     * @param vAvg - average element volume for elements with a facet at the interface to be estimated
     * @param vMin - minimum element volume at the interface to be estimated
    */
    template <class G>
    void updateGridAfterRemeshing(G& gridPtr,
                                  std::vector<int>& partitions,
                                  std::vector<int>& facets,
                                  std::vector<int>& boundaryIds) {
        //update data if grid has been swapped
        if (time_ > 0) {
            std::string name;
            if (count_ == 1) {// this case occurs, when we manually restart a simulation from a point after remeshing
                name = "newMeshRestore";
            } else {
                name = "newMesh";
            }
            MeshCreator<Grid> meshCreatorNew(name); //unrefined mesh
            auto gridPtrNew = meshCreatorNew.create();
            auto partitionsNew = meshCreatorNew.elementIds();

            partitions = std::move(partitionsNew);
            for (int & partition : partitions) partition = (partition > 0);

            facets.clear();
            facets = computeFacets(*gridPtr, partitions);

            gridPtr.swap(gridPtrNew);

            // store boundarySegmentIndices before restore since they might change
            std::vector<Edge> boundaryIndices = extractBoundaryIndices(gridPtr);

            //restore refinement information
            auto path = Parameters::get<std::string>("output directory").value_or("");
            std::string filename(path + "/backup/gridRefined" + std::to_string(time_));
            std::filesystem::path path2 = filename;
            std::ifstream file(path2);
            gridPtr->hostGrid()->restore(file);
            gridPtr->update();

            // now restore boundaryIds with new boundarySegmentIndex values
            restoreBoundaryIds(gridPtr,boundaryIndices,boundaryIds);

            int ref_int = Parameters::get<int>("refinement->interface").value_or(10);
            vMin_ = vAvg_ * 4 * std::pow(0.5, ref_int);
        }
    }

    /** \brief after remeshing, restore the data that are needed for the next time step (backuped in remeshing routine)
     * @param facets  - a vector containing facet indices of the facet touching the surface (-1 if element not on the surface)
     */
    void restoreDataAfterRemeshing(std::vector<int> const& facets) {
        auto& mainProb = mainProb_->get();
        if (time_ > 0) {
            auto deg = Parameters::get<int>("polynomial degree ch").value_or(2);
            auto const& gv = mainProb.gridView();

            DOFVector oldPhiGamma(gv, lagrange<1>());
            DOFVector oldMuGamma(gv, lagrange<1>());
            DOFVector oldVelocity(gv, power<2>(lagrange<2>(),flatInterleaved()));
            DOFVector oldPhi(gv, surface(lagrange(deg),facets),tag::no_datatransfer{});
            DOFVector oldMu(gv, surface(lagrange(deg),facets),tag::no_datatransfer{});

            auto path = Parameters::get<std::string>("output directory").value_or("");

            oldVelocity.restore(path + "/backup/uhOld" + std::to_string(time_));
            oldPhi.restore(path + "/backup/phiOld" + std::to_string(time_));
            oldMu.restore(path + "/backup/muOld" + std::to_string(time_));
            mainProb.oldLambda1()->restore(path + "/backup/lambda1Old" + std::to_string(time_));
            mainProb.oldLambda2()->restore(path + "/backup/lambda2Old" + std::to_string(time_));
            oldPhiGamma.restore(path + "/backup/phiGammaOld" + std::to_string(time_));
            oldMuGamma.restore(path + "/backup/muGammaOld" + std::to_string(time_));
            mainProb.solution(_v) << valueOf(oldVelocity);
            mainProb.solution(_phi) << valueOf(oldPhi);
            mainProb.solution(_mu) << valueOf(oldMu);
            mainProb.solution(_phiGamma) << valueOf(oldPhiGamma);
            mainProb.solution(_muGamma) << valueOf(oldMuGamma);

            remeshingJustHappened_ = 1;       
        } else {
            valueOf(*mainProb.oldLambda1()) << 1.0;
            valueOf(*mainProb.oldLambda2()) << 1.0;
        }
    }

    /// \brief after remeshing of at the beginning of the simulation, the time steps can be smaller, this is done here
    void adaptTimestepIfNecessary() {
        auto& mainProb = mainProb_->get();
        auto& adaptInfo = adaptInfo_->get();
        auto& tau = tau_->get();

        int reduceTimestep = Parameters::get<int>("reduce time step after remeshing").value_or(0);
        int startReduced = Parameters::get<int>("start with reduced time step").value_or(1);
        double reductionInitial= Parameters::get<double>("adapt->time step reduction factor init").value_or(1);
        int numberInitial = Parameters::get<int>("number of time steps for reduction init").value_or(3);
        double reduction = Parameters::get<double>("adapt->time step reduction factor").value_or(1);
        int number = Parameters::get<int>("number of time steps for reduction").value_or(3);

        // reduce the time step size if remeshing has happened in the last time step for number time steps
        if (remeshingJustHappened_ && reduceTimestep ||
            !remeshingJustHappened_ && startReduced && !countAfterRemeshing_) {
            if (remeshingJustHappened_)
                adaptInfo.setTimestep(tau * reduction);
            else if (startReduced)
                adaptInfo.setTimestep(tau * reductionInitial);
            mainProb.setMobility(Parameters::get<double>("parameters->M").value_or(1.0));
            remeshingJustHappened_ = 0;
            if (startReduced)
                mainProb.setPermeability(0.0);
        }
        countAfterRemeshing_++;

        if (countAfterRemeshing_ == number && !noRemeshingHappenedBefore_ ) {
            double step = Parameters::get<double>("adapt->timestep").value_or(1);
            mainProb.setMobility(Parameters::get<double>("parameters->M0").value_or(1.0));
            mainProb.setPermeability(Parameters::get<double>("parameters->permeability").value_or(0.0));
            adaptInfo.setTimestep(step);
            if (count_ < numberInitial && startReduced) {
                adaptInfo.setTimestep(tau * reductionInitial);
            }
        }
        if (count_ == numberInitial && startReduced) {
            int ref_int  = Parameters::get<int>("refinement->interface").value_or(10);
            vMin_ = vAvg_ * 4 * std::pow(0.5,ref_int);
            adaptInfo.setTimestep(Parameters::get<double>("adapt->timestep").value_or(1));
        }
    }

    /** \brief this is the main iteration of one time step including the main problem iteration and laplace iteration
     * @param kappaVecOld - the curvature vector of the current surface, to be computed here
     * @param y0 - the y coordinates of each DOF in the initial state (only used in axisymmetric case)
     */
    template <class KV, class DV>
    void iteration(KV& kappaVecOld, DV const& y0) {
        auto& mainProb = mainProb_->get();
        auto& laplaceProb = laplaceProb_->get();
        auto& extensions = extensions_->get();
        auto& adaptInfo = adaptInfo_->get();

        // this is the iteration of the main problem (nschProb)
        extensions.initTimestep(adaptInfo);
        mainProb.iteration(adaptInfo,kappaVecOld,y0);
        extensions.closeTimestep(adaptInfo);

        // print grid information
        printBulkGridInfo(*mainProb.grid());

        // update deltaP and deltaPhi on bulk grid
        updateForPermeability();

        // this is the laplace iteration for computation of the grid movement velocity
        auto ignoreLaplace = Parameters::get<int>("ignore laplace problem").value_or(0);
        if(!ignoreLaplace) laplaceProb.iteration(adaptInfo);

        // update laplace solution in main problem
        valueOf(*mainProb.laplaceSolution()) << laplaceProb.solution();
    }

    /// update methods

    /// \brief increase time by time step size and print time info
    void updateTime() {
        auto& adaptInfo = adaptInfo_->get();
        auto& tau = tau_->get();


        adaptInfo.incTimestepNumber();
        adaptInfo.setTime(adaptInfo.time() + tau);
        time_ = adaptInfo.time();

        // print time information
        std::cout << "time = " << adaptInfo.time() << ", end time = " << adaptInfo.endTime() << ", timestep = "
                  << tau << ", timestep number = " << count_ << "\n";

    }

    /** \brief update data structures needed for the stretching force
     * @param y0 - the y coordinates of each DOF in the initial state (only used in axisymmetric case)
     */
    template <class DV>
    void updateForStretching(DV const& y0) {
        auto& mainProb = mainProb_->get();

        auto axi = Parameters::get<int>("axisymmetric").value_or(0);
        auto Ks = Parameters::get<double>("areaShear").value_or(0);

        updateDistances(mainProb);
        if (axi && Ks) updateLambda2(mainProb, y0);
    }

    /// \brief update data structures needed for the permeability of the membrane
    void updateForPermeability() {
        auto& mainProb = mainProb_->get();

        assert(volume0_ > 0); // V0 has to be set before calling this method
        updateJumps(mainProb,volume0_);
    }

    /// \brief update the surface basis when the grid has been moved
    void updateBasis() {
        auto& mainProb = mainProb_->get();

        mainProb.globalBasis()->update(mainProb.gridView());
        mainProb.updateSurfaceBasis(); //update manually the data stored in each instance of the surface basis (automatically updated only for _xS)
    }

    /** \brief move the grid using the solution of the laplace problem
     * @param coords - a DOFVector containing the coordinates on the grid of the main problem (AdaptiveGrid)
     * @param coordinates - a pointer to a DOFVector containing the coordinates living on the underlying CurvedGrid
     */
    template <class DV, class DV2>
    void moveGrid(DV& coords, DV2& coordinates) {
        auto& mainProb = mainProb_->get();
        auto& laplaceProb = laplaceProb_->get();
        auto& tau = tau_->get();

        auto displacement = makeGridFunction(tau * laplaceProb.solution(), mainProb.gridView());
        valueOf(coords) << X() + displacement;
        interpolateHostGrid(valueOf(*coordinates),valueOf(coords));
    }

    /// printing methods

    /// \brief print info for bulk and surface grid and number of DOFs of the main problem
    template<class G, class SG>
    void printGridInfo(G const& grid, SG const& surfaceGrid) {
        auto& mainProb = mainProb_->get();

        printBulkGridInfo(grid);
        msg("main problem number of DOFs: {}", mainProb.globalBasis()->dimension());
        printSurfaceGridInfo(surfaceGrid);
    }

    /// \brief print info for bulk grid
    template<class G>
    void printBulkGridInfo(G const& grid) {
        msg("number of grid points: {}, number of faces: {}, number of elements: {}",
            grid.size(2), grid.size(1), grid.size(0));
    }

    /// \brief print info for surface grid
    template<class SG>
    void printSurfaceGridInfo(SG const& surfaceGrid) {
        msg("surfaceGrid: vertices={}, elements={}", surfaceGrid->size(1), surfaceGrid->size(0));
    }

    void printShellInfo() {
        auto& mainProb = mainProb_->get();

        std::cout << "shell area = " << mainProb.getShellArea() << ", shell velocity = " << mainProb.vGrid()[0]
                  << ", shell volume = "
                  << mainProb.volume() << "\n";
    }

    /// \brief print info for the end of the time step and the time needed for the whole time step
    void printEndOfTimestep(Dune::Timer const& timestepTime) {
        std::cout << "\n";
        AMDiS::info(1, "time step {} needed {} seconds", count_, timestepTime.elapsed());
        AMDiS::info(1, "====================================================================");
        std::cout << "\n";
        count_++;
    }

    double& vMin() {return vMin_;}

    int& count() {return count_;}

    void countAfterRemeshing(int num) {countAfterRemeshing_ = num;}

    void noRemeshingHappenedBefore(int num) {noRemeshingHappenedBefore_ = num;}

protected:
    std::optional<std::reference_wrapper<MainProb>> mainProb_; // the main problem
    std::optional<std::reference_wrapper<LProb>> laplaceProb_; // the laplace problem to move the grid
    std::optional<std::reference_wrapper<ExtensionManager>> extensions_; // the extensions of the main problem

    std::optional<std::reference_wrapper<AdaptInfo>> adaptInfo_; // the adaptation information
    std::optional<std::reference_wrapper<const double>> tau_; // the time step size
    double& time_; // the current time of the simulation
    int count_ = 1; // the number of the current time step
    int countAfterRemeshing_ = 0; // the number of the current time step after the last remeshing
    int remeshingJustHappened_ = 0; // inficator if remeshing happened in the last time step
    int noRemeshingHappenedBefore_ = 1; // indicator if remeshing has ever happened in this simulation
    double vMin_ = 15e15, vAvg_ = 0; // mimimum and average element volumes for refinement
    double volume0_ = -1, dropletVolume0_ = 0; // initial shell and droplet volume
};

#endif //ALEMODELAMDIS2_TIMEHANDLER_HPP