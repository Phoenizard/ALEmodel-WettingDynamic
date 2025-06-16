/**
 * Navier Stokes Cahn Hilliard Problem to be solved in every time step.
 *
 * The result is a 2D/3D velocity (_v) and a 1D pressure (_p). The basis for the velocity is of type lagrange<2>()
 * and the basis for the pressure is divided into subdomains, one for the internal and one for the external fluid
 * with respect to the shell (if the shell is closed). If the shell is not closed, the basis is defined using the
 * partitions vector, which then has to be modified s.t. it is equal to -1 apart from the surface.
 * The resulting velocity is continuous where the pressure is discontinuous over the shell.
 *
 * In addition to the Navier Stokes, a Cahn-Hilliard equation is solved in the system, resulting in a phase field,
 * which is between zero and 1. The Cahn-Hilliard equations are solved on the same basis as the pressure
 *
 * created by Marcel Mokbel 2020
**/

#pragma once

#include <memory>
#include <vector>

#include <amdis/extensions/BaseProblemInterface.hpp>
using namespace MyIndices;

namespace AMDiS
{
    template <class HostGrid, class PBF>
    class NavierStokesCahnHilliard
            : public BaseProblemInterface
            , public ProblemInstatBase
    {
    public:
        using WorldVector = FieldVector<double, HostGrid::dimensionworld>;
        using AdaptiveGrid = AdaptiveGrid_t<HostGrid>;
        using GridView = typename AdaptiveGrid::LeafGridView;

        using Basis = decltype(GlobalBasis{std::declval<GridView>(),std::declval<PBF>()});
        using Problem = ProblemStat<DefaultProblemTraits<Basis>>;// TYPEOF(ProblemStat{"prob", std::declval<HostGrid&>(), std::declval<GlobalBasis>()});
        using Super = BaseProblemInterface;

        using PBFDOW = decltype(power<Grid::dimensionworld>(lagrange<2>()));
        using VBasis = decltype(GlobalBasis{std::declval<GridView>(),std::declval<PBFDOW>()});
        using PBF1 = decltype(surface(lagrange<2>(),std::declval<std::vector<int>>()));
        using PhiBasis = decltype(GlobalBasis{std::declval<GridView>(),std::declval<PBF1>()});
        using PBFNu = decltype(surface(lagrange<1>(),std::declval<std::vector<int>>()));
        using NuBasis = decltype(GlobalBasis{std::declval<GridView>(),std::declval<PBFNu>()});
        using PBF0 = decltype(lagrange<1>());
        using Lambda2Basis = decltype(GlobalBasis{std::declval<GridView>(),std::declval<PBF0>()});
        using PBF2 = decltype(power<Grid::dimensionworld>(lagrange<1>()));
        using KappaVecBasis = decltype(GlobalBasis{std::declval<GridView>(),std::declval<PBF2>()});

    public:
        enum { dow = HostGrid::dimensionworld };
        enum { dim = HostGrid::dimension };

    public:
        NavierStokesCahnHilliard(std::string const& name,
                                 HostGrid& grid,
                                 PBF& pbf,
                                 std::vector<int> boundaryIDs,
                                 std::vector<int> partitions);

        ///Overriding Methods
        /// \brief Return the name of the BaseProblem. Implements \ref ProblemIterationInterface::name()
        std::string const& name() const override { return name_; }

        /// @biref Initialisation of the problem.
        void initialize(Flag initFlag,
                        BaseProblemInterface* adoptProblem = nullptr,
                        Flag adoptFlag = INIT_NOTHING) override
        {
            if (adoptProblem != nullptr) {
                auto* baseProblem = dynamic_cast<NavierStokesCahnHilliard*>(adoptProblem);
                test_exit(baseProblem != nullptr, "Can not determine BaseProblem type in initialization.");
                problem_.initialize(initFlag, &baseProblem->problem(), adoptFlag);
            } else {
                problem_.initialize(initFlag);
            }
        }

        /// @biref Method is called at the end of \ref initBaseProblem
        void finalizeData(AdaptInfo& /*adaptInfo*/) override { /* do nothing */ }

        /// \brief Implementation of \ref ProblemTimeInterface::transferInitialSolution()
        void transferInitialSolution(AdaptInfo& adaptInfo) override
        {
            writeFiles(adaptInfo);
        }


        /// \brief Implementation of \ref ProblemTimeInterface::closeTimestep()
        void closeTimestep(AdaptInfo& adaptInfo) override;

        /// \brief prepare iteration
        void beginIteration(AdaptInfo& /*adaptInfo*/) override
        {
            msg("");
            msg("[[ <{}> iteration ]]", name_);
        }

        /// \brief compute iteration
        template <class KV, class Y>
        void oneIteration(AdaptInfo& adaptInfo, KV& kappaVecOld, Y const& y0);

        /// \brief compute iteration
        Flag oneIteration(AdaptInfo& /*adaptInfo*/, Flag toDo = FULL_ITERATION) override
        {
            return 0;
        }

        /// \brief end iteration
        void endIteration(AdaptInfo& /*adaptInfo*/) override
        {
            msg("");
            msg("[[ end of <{}> iteration ]]", name_);
        }

        /// \brief full iteration
        template <class KV, class Y>
        void iteration(AdaptInfo& adaptInfo, KV& kappaVecOld, Y const& y0) {
            Dune::Timer t;
            this->initTimestep(adaptInfo);
            this->beginIteration(adaptInfo);
            this->oneIteration(adaptInfo,kappaVecOld,y0);
            this->endIteration(adaptInfo);
            this->closeTimestep(adaptInfo);
            AMDiS::info(2, "Bulk iteration needed {} seconds", t.elapsed());
        }

        /// \brief Calls writeFiles of the problem
        void writeFiles(AdaptInfo& adaptInfo) override
        {
            problem_.writeFiles(adaptInfo, false);
        }

        /// \brief intiialize the data for the current time step
        void initTimestep(AdaptInfo& adaptInfo) override;

        /// \brief intitalize the data for the simulation
        void initData(AdaptInfo& adaptInfo) override;

        /// \brief intitalize the data for the simulation
        template <class Lambda>
        void initMarkers(AdaptInfo& adaptInfo, Lambda const& lambda, double& vMin);

        /// \brief adds the operators of the equations to the system
        void fillOperators(AdaptInfo& adaptInfo) override;

        /// \brief adds the operators of the equations to the system
        void fillClassicOperators(AdaptInfo& adaptInfo);

        /// \brief adds the operators of the equations to the system
        void fillOptimizedOperators(AdaptInfo& adaptInfo);

        /// \brief adds the boundary conditions to the system
        void fillBoundaryConditions(AdaptInfo& adaptInfo) override;

        /// \brief builds the initial phase field and performs initial refinement
        void solveInitialProblem(AdaptInfo& adaptInfo) override;

        /// \brief assemble the whole system and solve it
        void assembleAndSolve(AdaptInfo& adaptInfo) {
            Flag toDo = FULL_ITERATION;
            Flag flag = 0, markFlag = 0;

            // assemble and solve
            if (toDo.isSet(BUILD))
                problem().buildAfterAdapt(adaptInfo, markFlag, true, true);

            if (toDo.isSet(BUILD_RHS))
                problem().buildAfterAdapt(adaptInfo, markFlag, false, true);

            // solve nschProb!
            if (toDo.isSet(SOLVE))
                problem().solve(adaptInfo, true, false);

            if (toDo.isSet(SOLVE_RHS))
                problem().solve(adaptInfo, true, false);

            if (toDo.isSet(ESTIMATE))
                problem().estimate(adaptInfo);
        }

        /** \brief update the surface basis when the grid has changed, update manually
         * if there is more than one instance of the surface basis. Manual update is used to prevent
         * multiple calls of extractSurfaceGrid() - this saves time.
         */
        void updateSurfaceBasis() {
            // update surfaceBasis instances manually!
            auto const& surfacePreBasis = this->globalBasis()->preBasis().subPreBasis(_xS).subPreBasis();
            auto sg = surfacePreBasis.surfaceGrid();
            auto stf = surfacePreBasis.surfaceToFluidMap();
            auto stfI = surfacePreBasis.surfaceToFluidMapIdx();

            updateSurfaceBasisManually<1>(problem(),sg,stf,stfI,
                                          _lambda1,
                                          _lambda2,
                                          _phiS,
                                          _phiGamma,
                                          _muGamma);

            updateSurfaceBasisManually<2>(problem(),sg,stf,stfI,
                                          _kappaVec,
                                          _fShell,
                                          _vS);
        }

        void buildAndAdapt(AdaptInfo& adaptInfo) {
            Flag toDo = FULL_ITERATION;
            Flag flag = 0, markFlag = 0;
            int adapted = 0;

            ///buildAndAdapt
            if (toDo.isSet(MARK))
                markFlag = problem().markElements(adaptInfo);

            if (toDo.isSet(ADAPT) && markFlag.isSet(MESH_ADAPTED)) {
                flag |= problem().adaptGrid(adaptInfo);
            }

            updateSurfaceBasis();
        }

        /// integrate over the shell part of the domain (only for closed shells)
        template <class Expr, class GV>
        auto integrateShell(Expr&& expr,GV const& gridView, std::vector<int> const& partitions, int order = 2, int part=-1);

    public:
        /// velocity vector components
        auto getVelocity(std::size_t i)       { return this->solution(_v,i); }
        auto getVelocity(std::size_t i) const { return this->solution(_v,i); }

        /// velocity vector
        auto getVelocity()       { return this->solution(_v); }
        auto getVelocity() const { return this->solution(_v); }

        /// Velocity vector components at old timestep
        auto getOldVelocity(std::size_t i)        { return valueOf(*oldVelocity_,i); }
        auto getOldVelocity(std::size_t i) const  { return valueOf(*oldVelocity_,i); }

        /// velocity vector at old timestep
        auto getOldVelocity()        { return valueOf(*oldVelocity_); }
        auto getOldVelocity() const  { return valueOf(*oldVelocity_); }

        /// phase field component
        auto getPhase(int = 0)       { return this->solution(_phi);  }
        auto getPhase(int = 0) const { return this->solution(_phi);  }

        /// phase field component at old time step
        auto getOldPhase(int = 0)       { return valueOf(*oldPhi_); }
        auto getOldPhase(int = 0) const { return valueOf(*oldPhi_); }

        /// pressure component
        auto getPressure()       { return this->solution(_p); }
        auto getPressure() const { return this->solution(_p); }

        /// mu component
        auto getW(int = 0)       { return valueOf(*oldMu_); }
        auto getW(int = 0) const { return valueOf(*oldMu_); }

        /// Return a pointer to the grid
        std::shared_ptr<AdaptiveGrid>       grid()       { return problem_.grid(); }
        std::shared_ptr<AdaptiveGrid const> grid() const { return problem_.grid(); }

        /// Return the gridView of the leaf-level
        GridView gridView() const { return problem_.gridView(); }

        /// Return the boundary manager to identify boundary segments
        std::shared_ptr<BoundaryManager<AdaptiveGrid>>       boundaryManager()       { return problem_.boundaryManager(); }
        std::shared_ptr<BoundaryManager<AdaptiveGrid> const> boundaryManager() const { return problem_.boundaryManager(); }

        /// Return the GlobalBasis of the Problem
        std::shared_ptr<Basis>       globalBasis()       { return problem_.globalBasis(); }
        std::shared_ptr<Basis const> globalBasis() const { return problem_.globalBasis(); }

        /// Return the number of problems this BasieProblem contains. Implement \ref ProblemIterationInterface::numProblems()
        int numProblems() const override { return 1; }

        /// Return the \ref problem_
        Problem&       problem(int = 0)       { return problem_; }
        Problem const& problem(int = 0) const { return problem_; }

        Problem&       problem(std::string const&)       { return problem_; }
        Problem const& problem(std::string const&) const { return problem_; }

        /// Return a mutable/const view to a solution component
        template <class TreePath = RootTreePath>
        decltype(auto) solution(TreePath path = {}) { return problem_.solution(path); }

        template <class TreePath = RootTreePath>
        decltype(auto) solution(TreePath path = {}) const { return problem_.solution(path); }

        /// compute density on the 2/3 phases - right now return constant density!
        auto rhoPhase() {
            auto rhoInt = Parameters::get<double>("stokes->density").value_or(0.0);
            auto rhoExt = Parameters::get<double>("stokes->density2").value_or(0.0);
            auto rhoOut = Parameters::get<double>("stokes->densityOutside").value_or(0.0);

            //valueOf(*rhoDOF_,_phi) << clamp(getPhase(),0,1)*rhoInt + (1.0 - clamp(getPhase(),0,1))*rhoExt ;
            return rhoInt; //valueOf(*rhoDOF_,_phi);
        }

        /// compute viscosity on the 2/3 phases
        auto nuPhase() {
            auto nuInt = Parameters::get<double>("stokes->viscosityPhase").value_or(0.0);
            auto nuExt = Parameters::get<double>("stokes->viscosityOutside").value_or(0.0);
            auto nuInExt = Parameters::get<double>("stokes->viscosityInside").value_or(0.0);

            auto elementVectorIn = ElementVector(*grid(),0);
            auto elementVectorOut = ElementVector(*grid(),0);
            for (auto const& e : elements(gridView())) {
                auto father = e;
                while (father.hasFather())
                    father = father.father();

                if (partitions_[problem().grid()->levelGridView(0).indexSet().index(father)] == 1) {
                    elementVectorIn.data()[gridView().indexSet().index(e)] = 1;
                }
                if (partitions_[problem().grid()->levelGridView(0).indexSet().index(father)] == 0) {
                    elementVectorOut.data()[gridView().indexSet().index(e)] = 1;
                }

            }

            valueOf(*nuDOF_) << clamp(getOldPhase(),0.0,2.0)*nuInt
                                     + (1.0 - clamp(getOldPhase(),0.0,2.0))*nuExt*valueOf(elementVectorOut)
                                       + (1.0 - clamp(getOldPhase(),0.0,2.0))*nuInExt*valueOf(elementVectorIn);
            return valueOf(*nuDOF_);
        }

        /// update the Elementvectors for element size and element refinement level
        void updateElementVectors(bool all = false) {
            elementSizes.resizeZero();
            elementLevels.resizeZero();
            auto maxLevel = 0.0;
            avgElVol_ = 0;
            int count = 0;
            for (auto const &e: elements(gridView())) {
                elementSizes.data()[gridView().indexSet().index(e)] = e.geometry().volume();
                elementLevels.data()[gridView().indexSet().index(e)] = e.level();
                if (e.geometry().volume() < minimumElementVolume) minimumElementVolume = e.geometry().volume();
                avgElVol_ += e.geometry().volume();
                count++;
                if (all && e.level() > maxLevel) maxLevel = e.level();
            }
            avgElVol_ /= count;
            if (all) std::cout << "maximum refinement level on grid = " << maxLevel << "\n";
        }

        /// set the phase field mobility to mob
        void setMobility(double mob) {
            mobility_ = mob;
        }

        /// set the permeability to P
        void setPermeability(double P) {
            permeability_ = P;
        }

        ///set the initial shell area (only for closed shells)
        void setShellArea0(double area) {
            shellArea0_ = area;
        }

        ///set the initial droplet volume area
        void setDropletVolume0(double vol) {
            dropletVolume0_ = vol;
        }

        ///set the initial shell volume area
        void setVolume0(double vol) {
            volume0_ = vol;
        }

        void updateNormal();

        /// get functions

        /// get the average velocity of the grid (only for closed shells)
        auto const& vGrid() {return vGrid_;}

        /// get the shell volume (only for closed shells)
        double const& volume() {return volume_;}

        /// get the shell area (only for closed shells)
        double const& getShellArea() {return shellArea_;}

        /// get the initial droplet volume
        double const& getDropletVolume0() {return dropletVolume0_;}

        double& dropletVolumeChange() {return dropVolChangeOverTime_;}

        /// get the initial shell volume
        double const& getVolume0() {return volume0_;}

        /// get the min value for element sizes on the mesh
        double getMinimumElementVolume() {return minimumElementVolume;}

        /// get the partitions vector
        auto const& partitions() {return partitions_;}
        auto& partitionsRef() {return partitions_;}

        /// get the DOFVector for the pressure jump across Gamma
        auto& deltaP() {return deltaP_;}

        /// get the DOFVector for the phase field jump across Gamma
        auto& deltaPhi() {return deltaPhi_;}

        /// get the ElementVector for the principal stretches on Gamma (if axi, then lambda1*lambda2, else lambda1)
        auto& distances() {return distances_;}

        /// get the DOFVector for the principal stretches on Gamma (if axi, then lambda1*lambda2, else lambda1)
        auto& distancesDOF() {return distancesDOF_;}

        /// get the DOFVector for lambda1 on Gamma, which was stored after last remeshing
        auto& oldLambda1() {return oldLambda1_;}

        /// get the DOFVector for lambda2 on Gamma, which was stored after last remeshing
        auto& oldLambda2() {return oldLambda2_;}

        /// get the DOFVector for the principal stretches on Gamma (lambda2)
        auto& lambda2() {return lambda2_;}

        /// get the DOFVector for the solution of the problem at last time step
        auto& oldSolution() {return oldSolution_;}

        /// get the DOFVector for the curvature vector of the current shape of Gamma
        auto& oldKappaVec() {return oldKappaVec_;}

        /// get the DOFVector for the grid velocity
        auto& laplaceSolution() {return laplaceSolution_;}

        /// get the vector of ElementVectors for the surface normal
        auto const& normal() {return normalVec_;}

        /// get the DOFVector for the surface normal
        auto const& normalDOF() {return normalDOF_;}

        /// get the Basis for the pressure/phi space (with jump across Gamma
        auto const& nuBasis() {return nuBasis_;}

        auto const& getElementSizes() {return elementSizes;}

        auto const& getElementLevels() {return elementLevels;}

        auto const& getViscosityDOF() {return *nuDOF_;}

    protected:
        double rho_ = 1.0, rho2_ = 1.0, rhoOut_ = 1.0;              // densities
        double nu_ = 1.0;                                           // fluid visosity
        double eps_;                                                // epsilon of Cahn-Hilliard equation
        double sigma_;                                              // surface tension of CH
        double B_;                                                  // bending stiffness
        double Ks_, Ka_;                                            // area shear, area dilation
        double kappa0_;                                             // spontaneous curvature
        double mobility_;                                           // Cahn-Hilliard mobility (in bulk and on surface)
        double permeability_;                                       // the shell's permeability
        double volume0_;                                            // initial shell volume (only if closes shell is used)
        double volume_, volChangeOverTime_;                         // shell volume (only if closed shell is used)
        FieldVector<double,Grid::dimensionworld> vGrid_;            // grid velocity
        unsigned int axi;                                           // boolean if axisymmetry is used
        double E_sum;                                               // Gesamtenergie
        double shellArea_, shellArea0_, shellAreaChange_ = 0.0;
        double dropletVolume_,dropletVolume0_, dropVolChangeOverTime_, minimumElementVolume = 1e15;
        double avgElVol_;
        std::vector<int> boundaryIDs_, partitions_;                 // information gathered from the mesh: IDs of boundary
                                                                    // segments and elemtent partition IDs

        std::shared_ptr<VBasis> vBasis_;                            // basis for velocity
        std::shared_ptr<PhiBasis> phiBasis_;                        // basis for phi
        std::shared_ptr<NuBasis> nuBasis_;                          // basis for viscosity
        std::shared_ptr<Lambda2Basis> lambda2Basis_;                // basis for lambda2
        std::shared_ptr<KappaVecBasis> kappaVecBasis_;              // basis for curvature vector

        std::shared_ptr<DOFVector<Basis>> oldSolution_;             // DOFVector for oldSolution
        std::shared_ptr<DOFVector<Basis>> rhoDOF_;                  // DOFVector for density
        std::shared_ptr<DOFVector<NuBasis>> nuDOF_;                 // DOFVector for viscosity
        std::shared_ptr<DOFVector<VBasis>> oldVelocity_;            // DOFVector for velocity
        std::shared_ptr<DOFVector<PhiBasis>> oldPhi_;               // DOFVector for phase
        std::shared_ptr<DOFVector<KappaVecBasis>> oldKappaVec_;     // DOFVector for old phiGamma
        std::shared_ptr<DOFVector<KappaVecBasis>> laplaceSolution_; // DOFVector for old phiGamma
        std::shared_ptr<DOFVector<KappaVecBasis>> normalDOF_;       // DOFVector for normal
        std::shared_ptr<DOFVector<PhiBasis>> oldMu_;                // DOFVector for mu
        std::shared_ptr<DOFVector<Lambda2Basis>> lambda2_, deltaP_, deltaPhi_; // DOFVectors for lambda2, pressure and phi jump
        std::shared_ptr<DOFVector<Lambda2Basis>> distancesDOF_;     // DOFVectors for distances
        std::shared_ptr<DOFVector<Lambda2Basis>> oldLambda1_, oldLambda2_; // DOFVectors for lambdas that were stored from lasr remeshing

        Problem problem_;                                           // the underlying ProblemStat
        std::string name_;                                          // name of the problem
        double invTau_, tau_;                                       // (inverse of) the time step
        ElementVector<HostGrid,int> elementLevels;                  // element vector storing refinement levels
        ElementVector<HostGrid,double> elementSizes, distances_;    // element vectors storing element volume and principal stretches
        using EV = ElementVector<typename GridView::Grid,double>;
        std::vector<EV> normalVec_;                                 // normal vector to the surface pointing outwards
                                                                    // from the respective element touching the surface
    };

} // end namespace AMDiS

#include "NavierStokesCahnHilliard.impl.hpp"
