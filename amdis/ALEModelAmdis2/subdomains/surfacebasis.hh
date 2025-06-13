// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_SUBDOMAINS_SURFACEBASIS_HH
#define DUNE_SUBDOMAINS_SURFACEBASIS_HH

#include <vector>

#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/indices.hh>

#include <dune/functions/common/utility.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/concepts.hh>

#include "surfacegridview.hh"
#include "surfaceindexset.hh"


namespace Dune {
namespace Subdomains {

template <class PB, class NodeIndexSet>
class SurfaceNodeIndexSet;

template <class GV, class PreBasis>
class SurfacePreBasis
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

  //! Type used for indices and size information
  using size_type = std::size_t;

public:
  //! Constructor for given pre-basis objects
  template <class PB>
  SurfacePreBasis (PB&& preBasis, std::shared_ptr<PIndexSet> indexSet)
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

  //! Update the stored grid view, to be called if the grid has changed
  void update (const GridView& gv)
  {
    indexSet_->update(gv);
    PreBasis::update(PGridView{gv, indexSet_});
  }

  template <class Entity>
  void setSurface (Entity const& entity, int facet)
  {
    static_assert(Entity::codimension == 0);
    indexSet_->setSurface(entity, facet);
  }

private:
  std::shared_ptr<PIndexSet> indexSet_;
};

template <class PB, class NodeIndexSet>
class SurfaceNodeIndexSet
: public NodeIndexSet
{
public:
    using PreBasis = PB;

public:
    SurfaceNodeIndexSet (NodeIndexSet const& nodeIndexSet)
    : NodeIndexSet{nodeIndexSet}
    {}

    SurfaceNodeIndexSet (NodeIndexSet&& nodeIndexSet)
    : NodeIndexSet{std::move(nodeIndexSet)}
    {}
};
} // end namespace Subdomains


namespace Functions {
namespace BasisFactory {
namespace Imp {

template <class PreBasisFactory>
class SurfacePreBasisFactory
{
public:

  template <class PBF,
    disableCopyMove<SurfacePreBasisFactory,PBF> = 0>
  SurfacePreBasisFactory(PBF&& preBasisFactory)
    : preBasisFactory_(std::forward<PBF>(preBasisFactory))
  {}

  template <class PBF>
  SurfacePreBasisFactory(PBF&& preBasisFactory, std::vector<int>const& facets)
    : preBasisFactory_(std::forward<PBF>(preBasisFactory))
    , facets_(std::move(facets))
  {}

  template <class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
      using PIndexSet = Subdomains::SurfaceIndexSet<typename GridView::IndexSet,typename GridView::Grid::LevelGridView::IndexSet>;
      auto pIndexSet = std::make_shared<PIndexSet>(gridView.indexSet(),gridView.grid().levelGridView(0).indexSet(), facets_);

      using PGridView = Subdomains::SurfaceGridView<GridView>;
      auto preBasis = preBasisFactory_(PGridView{gridView, pIndexSet});

      using PreBasis = Subdomains::SurfacePreBasis<GridView, decltype(preBasis)>;
      return PreBasis{std::move(preBasis), pIndexSet};
  }

  template <class GridView>
  auto operator()(const GridView& gridView) const
  {
      using PIndexSet = Subdomains::SurfaceIndexSet<typename GridView::IndexSet,typename GridView::Grid::LevelGridView::IndexSet>;
      auto pIndexSet = std::make_shared<PIndexSet>(gridView.indexSet(),gridView.grid().levelGridView(0).indexSet(), facets_);

      using PGridView = Subdomains::SurfaceGridView<GridView>;
      auto preBasis = preBasisFactory_(PGridView{gridView, pIndexSet});

      using PreBasis = Subdomains::SurfacePreBasis<GridView, decltype(preBasis)>;
      return PreBasis{std::move(preBasis), pIndexSet};
  }

private:
  PreBasisFactory preBasisFactory_;
  std::vector<int> const& facets_;
};

} // end namespace BasisFactory::Imp


/**
 * \brief Create a pre-basis factory that can build a SurfacePreBasis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam numPartitions   Number of partitions
 * \param preBasisFactory  Pre-basis factory to wrap
 */
template <class PreBasisFactory>
auto surface(PreBasisFactory const& preBasisFactory)
{
  using Factory = Imp::SurfacePreBasisFactory<PreBasisFactory>;
  return Factory{preBasisFactory};
}

template <class PreBasisFactory>
auto surface(PreBasisFactory const& preBasisFactory, std::vector<int> const& facets)
{
  using Factory = Imp::SurfacePreBasisFactory<PreBasisFactory>;
  return Factory{preBasisFactory, std::move(facets)};
}

} // end namespace BasisFactory
} // end namespace Functions
} // end namespace Dune


#endif // DUNE_SUBDOMAINS_SURFACEBASIS_HH