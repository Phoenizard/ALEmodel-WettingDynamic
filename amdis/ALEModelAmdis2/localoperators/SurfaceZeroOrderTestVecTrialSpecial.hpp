#pragma once

#include <type_traits>

#include "SurfaceZeroOrderTestTrialvecSpecial.hpp"

namespace AMDiS
{
    /**
     * \addtogroup operators
     * @{
     **/

    namespace tag
    {
        template <class IS = int>
        struct surface_testvec_trial_special {
            int comp;
            std::size_t inside = 0; // integrate over the shell from inside (0) or outside (1)
            IS const &indexSet = {}; //the indexSet of the leafGridView of the grid!
            std::vector<int> const&  partitions = {}; //partition numbers of each element
            std::vector<int> const&  facets = {}; //facet number of the surface facet (-1 if element i not touching the surface)
        };
    }


    /// zero-order operator \f$ \langle\Psi, \mathbf{b}\,\phi\rangle \f$
    template <class LC, class IS>
    struct GridFunctionOperatorRegistry<tag::surface_testvec_trial_special<IS>, LC>
    {
        static constexpr int degree = 0;
        static tag::surface_test_trialvec_special<IS> transposedTag(tag::surface_testvec_trial_special<IS> t)
        {
            return {t.comp, t.inside, t.indexSet, t.partitions, t.facets};
        }
        using type = SurfaceZeroOrderTestTrialvecSpecial<IS>;
    };

    /** @} **/

} // end namespace AMDiS
