// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_SUBDOMAINBASIS_HH
#define DUNE_SUBDOMAINBASIS_HH

#include <vector>

#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/indices.hh>

#include <dune/functions/common/utility.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/concepts.hh>

#include "subdomaingridview.hh"
#include "subdomainindexset.hh"


namespace Dune {
namespace Subdomains {

template <class PB, class NodeIndexSet>
class SubdomainNodeIndexSet;

template <class GV, class PreBasis>
class SubdomainPreBasis
    : public PreBasis
{
  using PGridView = typename PreBasis::GridView;
  using PIndexSet = typename PGridView::IndexSet;

public:
  static const std::size_t maxMultiIndexSize = PreBasis::maxMultiIndexSize;
  static const std::size_t minMultiIndexSize = PreBasis::minMultiIndexSize;
  static const std::size_t multiIndexBufferSize = PreBasis::multiIndexBufferSize;
  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Template mapping root tree path to type of created tree node index set
  //using IndexSet = SubdomainNodeIndexSet<SubdomainPreBasis, typename PreBasis::IndexSet>;

public:
  //! Constructor for given pre-basis objects
  template <class PB>
  SubdomainPreBasis (PB&& preBasis, std::shared_ptr<PIndexSet> indexSet)
    : PreBasis(std::forward<PB>(preBasis))
    , indexSet_(std::move(indexSet))
  {}

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView () const
  {
    return PreBasis::gridView().impl();
  }

  //! Initialize the global indices
  void initializeIndices ()
  {
      indexSet_->update(gridView());
      PreBasis::initializeIndices();
  }

  //! Initialize the global indices
  void initializeIndices (std::vector<int> newSubdomains)
  {
      indexSet_->update(gridView(), newSubdomains);
      PreBasis::initializeIndices();
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update (const GridView& gv)
  {
    indexSet_->update(gv);
    PreBasis::update(PGridView{gv, indexSet_});
  }

  template <class Entity>
  void setSubdomain (Entity const& entity, int partition)
  {
    static_assert(Entity::codimension == 0);
    indexSet_->setSubdomain(entity, partition);
  }

private:
  std::shared_ptr<PIndexSet> indexSet_;
};



template <class PB, class NodeIndexSet>
class SubdomainNodeIndexSet
    : public NodeIndexSet
{
public:
  using PreBasis = PB;

public:
  SubdomainNodeIndexSet (NodeIndexSet const& nodeIndexSet)
    : NodeIndexSet{nodeIndexSet}
  {}

  SubdomainNodeIndexSet (NodeIndexSet&& nodeIndexSet)
    : NodeIndexSet{std::move(nodeIndexSet)}
  {}
};


} // end namespace Subdomains


namespace Functions {
namespace BasisFactory {
namespace Imp {

template <class PreBasisFactory>
class SubdomainPreBasisFactory
{
public:
  //static const std::size_t requiredMultiIndexSize = std::decay_t<PreBasisFactory>::minMultiIndexSize;

  template <class PBF,
    disableCopyMove<SubdomainPreBasisFactory,PBF> = 0>
  SubdomainPreBasisFactory(PBF&& preBasisFactory)
    : preBasisFactory_(std::forward<PBF>(preBasisFactory))
  {}

  template <class PBF>
  SubdomainPreBasisFactory(PBF&& preBasisFactory, std::vector<int> const& partitions)
    : preBasisFactory_(std::forward<PBF>(preBasisFactory))
    , subdomains_(std::move(partitions))
  {}

  template <class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
      using PIndexSet = Subdomains::SubdomainIndexSet<typename GridView::IndexSet>;
      auto pIndexSet = std::make_shared<PIndexSet>(gridView.indexSet(), subdomains_);

      using PGridView = Subdomains::SubdomainGridView<GridView>;
      auto preBasis = preBasisFactory_(PGridView{gridView, pIndexSet});

      using PreBasis = Subdomains::SubdomainPreBasis<GridView, decltype(preBasis)>;
      return PreBasis{std::move(preBasis), pIndexSet};
  }

  template <class GridView>
  auto operator()(const GridView& gridView) const
  {
      using PIndexSet = Subdomains::SubdomainIndexSet<typename GridView::IndexSet>;
      auto pIndexSet = std::make_shared<PIndexSet>(gridView.indexSet(), subdomains_);

      using PGridView = Subdomains::SubdomainGridView<GridView>;
      auto preBasis = preBasisFactory_(PGridView{gridView, pIndexSet});

      using PreBasis = Subdomains::SubdomainPreBasis<GridView, decltype(preBasis)>;
      return PreBasis{std::move(preBasis), pIndexSet};
  }

private:
  PreBasisFactory preBasisFactory_;
  std::vector<int> const& subdomains_;
};

} // end namespace BasisFactory::Imp


/**
 * \brief Create a pre-basis factory that can build a SubdomainPreBasis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam numPartitions   Number of partitions
 * \param preBasisFactory  Pre-basis factory to wrap
 */
template <class PreBasisFactory>
auto subdomains(PreBasisFactory const& preBasisFactory)
{
  using Factory = Imp::SubdomainPreBasisFactory<PreBasisFactory>;
  return Factory{preBasisFactory};
}

template <class PreBasisFactory>
auto subdomains(PreBasisFactory const& preBasisFactory, std::vector<int> const& partitions)
{
  using Factory = Imp::SubdomainPreBasisFactory<PreBasisFactory>;
  return Factory{preBasisFactory, std::move(partitions)};
}

} // end namespace BasisFactory
} // end namespace Functions
} // end namespace Dune


#endif // DUNE_SUBDOMAINBASIS_HH