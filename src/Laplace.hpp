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

namespace AMDiS
{
    template <class Grid>
    using LaplaceTraits = LagrangeBasis<Grid,1,1>;

    template <class Grid>
    class Laplace
            : public BaseProblem<LaplaceTraits<Grid>>
    {
        using Super = BaseProblem<LaplaceTraits<Grid>>;
        using AdaptiveGrid = AdaptiveGrid_t<Grid>;

    public:
        using GlobalBasis = typename Super::GlobalBasis;

    public:
        Laplace(std::string const& name, std::shared_ptr<AdaptiveGrid> grid, std::vector<int> partitions_, std::vector<int> boundaryIDs_);

        void initTimestep(AdaptInfo& adaptInfo) override;

        void initData(AdaptInfo& adaptInfo) override;
        void fillOperators(AdaptInfo& adaptInfo) override;
        void fillBoundaryConditions(AdaptInfo& adaptInfo) override;
        void iteration(AdaptInfo& adaptInfo) {
            Dune::Timer t;
            this->initTimestep(adaptInfo);
            this->beginIteration(adaptInfo);
            this->oneIteration(adaptInfo);
            this->endIteration(adaptInfo);
            this->closeTimestep(adaptInfo);
            AMDiS::info(2, "Laplace iteration needed {} seconds", t.elapsed());
        }

        using Super::problem;

    protected:
        int n_ = 0;
        std::vector<int> partitions, boundaryIDs;
    };

} // end namespace AMDiS

#include "Laplace.impl.hpp"
