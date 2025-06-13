/**
 * Laplace Problem to be solved in every time step.
 * The shell surface is identified and the movement on the shell is given as a dirichlet boundary condition.
 * This condition can be explicitly set by the user using the MeshMovementBC class. It can also be given by the
 * solution of a Navier-Stokes problem (probably also with shell forces).
 *
 * The result is a displacement, which can be used to update the mesh (the mesh has to be a Dune::CurvedGrid<HostGrid>)
 *
 * created by Marcel Mokbel 2020
**/

#pragma once

#include <memory>

#include <amdis/extensions/BaseProblem.hpp>

#include <dune/grid-glue/adapter/gridgluevtkwriter.hh>
#include <dune/grid-glue/extractors/codim0extractor.hh>
#include <dune/grid-glue/gridglue.hh>
#include <dune/grid-glue/merging/overlappingmerge.hh>

namespace AMDiS
{
    template <class Grid, class PBF>
    class InterpolateGrids
            : public BaseProblemInterface
            , public ProblemInstatBase
    {
    public:
        using AdaptiveGrid = AdaptiveGrid_t<Grid>;
        using GridView = typename AdaptiveGrid::LeafGridView;

        using Basis = decltype(GlobalBasis{std::declval<GridView>(),std::declval<PBF>()});
        using Problem = ProblemStat<DefaultProblemTraits<Basis>>;
        using Super = BaseProblemInterface;

        enum { dow = Grid::dimensionworld };
        enum { dim = Grid::dimension };

        static constexpr auto _v = Dune::Indices::_0;
        static constexpr auto _lambda1 = Dune::Indices::_1; //lagrange multiplier for v
        static constexpr auto _phi = Dune::Indices::_2;
        static constexpr auto _mu = Dune::Indices::_3;

    public:
        InterpolateGrids(std::string const& name, Grid& grid, Grid& gridOld, PBF& pbf, double const& vMin);

        ///Overriding Methods
        /// @brief Return the name of the BaseProblem. Implements \ref ProblemIterationInterface::name()
        std::string const& name() const override { return name_; }

        /// @biref Initialisation of the problem.
        void initialize(Flag initFlag,
                        BaseProblemInterface* adoptProblem = nullptr,
                        Flag adoptFlag = INIT_NOTHING) override
        {
            if (adoptProblem != nullptr) {
                auto* baseProblem = dynamic_cast<InterpolateGrids*>(adoptProblem);
                test_exit(baseProblem != nullptr, "Can not determine BaseProblem type in initialization.");
                problem_.initialize(initFlag, &baseProblem->problem(), adoptFlag);
            } else {
                problem_.initialize(initFlag);
            }
        }

        /// @biref Method is called at the end of \ref initBaseProblem
        void finalizeData(AdaptInfo& /*adaptInfo*/) override { /* do nothing */ }

        /// @brief Implementation of \ref ProblemTimeInterface::transferInitialSolution()
        void transferInitialSolution(AdaptInfo& adaptInfo) override
        {
            writeFiles(adaptInfo);
        }


        /// @brief Implementation of \ref ProblemTimeInterface::closeTimestep()
        void closeTimestep(AdaptInfo& adaptInfo) override;

        /// @brief prepare iteration
        void beginIteration(AdaptInfo& /*adaptInfo*/) override
        {
            msg("");
            msg("[[ <{}> iteration ]]", name_);
        }

        /// @brief compute iteration
        Flag oneIteration(AdaptInfo& adaptInfo, Flag toDo = FULL_ITERATION) override
        {
            Flag flag;

            problem_.beginIteration(adaptInfo);
            flag |= problem_.oneIteration(adaptInfo, toDo);
            problem_.endIteration(adaptInfo);

            return flag;
        }

        /// @brief end iteration
        void endIteration(AdaptInfo& /*adaptInfo*/) override
        {
            msg("");
            msg("[[ end of <{}> iteration ]]", name_);
        }

        /// @brief Calls writeFiles of the problem
        void writeFiles(AdaptInfo& adaptInfo) override
        {
            problem_.writeFiles(adaptInfo, false);
        }

        /// @brief intiialize the data for the current time step
        void initTimestep(AdaptInfo& adaptInfo) override;

        /// @brief intitalize the data for the simulation
        void initData(AdaptInfo& adaptInfo) override;

        /// @brief adds the operators of the equations to the system
        void fillOperators(AdaptInfo& adaptInfo) override;

        /// @bried adds the boundary conditions to the system
        void fillBoundaryConditions(AdaptInfo& adaptInfo) override;

        /// @brief builds the initial phase field and performs initial refinement
        void solveInitialProblem(AdaptInfo& adaptInfo) override;

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

        /// phase field component
        auto getPhase(int = 0)       { return this->solution(_phi);  }
        auto getPhase(int = 0) const { return this->solution(_phi);  }

        using Super::problem;
        using Element = typename Grid::template Codim<0>::Entity;
        using OverlappingMerge = Dune::GridGlue::OverlappingMerge<dim,dim,dim>;
        using Extractor = Dune::GridGlue::Codim0Extractor<typename Grid::LeafGridView>;
        using Glue = Dune::GridGlue::GridGlue<Extractor, Extractor>;

    public:
        int n_ = 0;
        Grid const *grid_, *gridOld_;
        Glue glue_;
        Problem problem_;                                           // the underlying ProblemStat
        std::string name_;
        ElementVector<Grid,int> elementLevels;
        ElementVector<Grid,double> elementSizes;
        double vMin_; //minimum element size
    };

} // end namespace AMDiS

#include "InterpolateGrids.impl.hpp"
