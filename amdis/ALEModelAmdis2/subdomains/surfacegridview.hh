#ifndef DUNE_SUBDOMAINS_SURFACEGRIDVIEW_HH
#define DUNE_SUBDOMAINS_SURFACEGRIDVIEW_HH

#include <memory>
#include <type_traits>
#include <utility>

#include "surfaceindexset.hh"

namespace Dune::Subdomains {

template <class GV>
class SurfaceGridView
{
  using GridView = GV;

public:
  /** \brief type of the grid */
  using Grid = typename GridView::Grid;

  /** \brief type of the index set */
  using IndexSet = SurfaceIndexSet<typename GridView::IndexSet,typename GridView::Grid::LevelGridView::IndexSet>;

  // /** \brief type of the intersection */
   using Intersection = typename GridView::Intersection;

  // /** \brief type of the intersection iterator */
  using IntersectionIterator = typename GridView::IntersectionIterator;

  /** \brief type of the collective communication */
  using Communication = typename GridView::Communication;

  /** \brief A struct that collects all associated types from the Traits class. */
  template <int cd>
  using Codim = typename GridView::template Codim<cd>;

  /** \brief type used for coordinates in grid */
  using ctype = typename Grid::ctype;

  /** \brief EntitySeed type of grid elements */
  using Element = typename Grid::template Codim<0>::Entity;

  /** \brief Traits class */
  using Traits = typename GridView::Traits;

  enum {
    dimension = Grid::dimension,          //< The dimension of the grid
    dimensionworld = Grid::dimensionworld,  //< The dimension of the world the grid lives in
    conforming = Traits::conforming         //< Export if this grid view is conforming
  };

public:
  SurfaceGridView(GridView const& gridView, std::shared_ptr<IndexSet> indexSet)
    : gridView_(gridView)
    , indexSet_(std::move(indexSet))
  {}

public:
  /** \brief obtain a const reference to the underlying hierarchic grid */
  Grid const& grid () const
  {
    return impl().grid();
  }

  /** \brief obtain the index set
   *
   * The lifetime of the returned index set is bound to the lifetime of the
   * grid view. Keep a copy of the grid view to prevent the index set from
   * becoming a dangling reference.
   */
  IndexSet const& indexSet () const
  {
    return *indexSet_;
  }

  /** \brief obtain number of entities in a given codimension */
  int size (int codim) const
  {
    return indexSet().size(codim);
  }

  /** \brief obtain number of entities with a given geometry type */
  int size (GeometryType const& type) const
  {
    return indexSet().size(type);
  }

  /** \brief Return true if the given entity is contained in this grid view
   *
   * \note If the input element e is not an element of the grid, then
   *       the result of contains() is undefined.
   */
  template <class EntityType>
  bool contains (EntityType const& e) const
  {
    return impl().indexSet().contains(e);
  }

  /** \brief obtain begin iterator for this view */
  template <int cd>
  typename Codim<cd>::Iterator begin () const
  {
    return impl().template begin<cd>();
  }

  /** \brief obtain end iterator for this view */
  template <int cd>
  typename Codim<cd>::Iterator end () const
  {
    return impl().template end<cd>();
  }

  /** \brief obtain begin iterator for this view */
  template <int cd, PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator begin () const
  {
    return impl().template begin<cd,pitype>();
  }

  /** \brief obtain end iterator for this view */
  template <int cd, PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator end () const
  {
    return impl().template end<cd,pitype>();
  }

  /** \brief obtain begin intersection iterator with respect to this view */
  IntersectionIterator ibegin (typename Codim<0>::Entity const& entity) const
  {
    return impl().ibegin(entity);
  }

  /** \brief obtain end intersection iterator with respect to this view */
  IntersectionIterator iend (typename Codim<0>::Entity const& entity) const
  {
    return impl().iend(entity);
  }

  /** \brief obtain collective communication object */
  Communication const& comm () const
  {
    return impl().comm();
  }

  /** \brief Return size of the overlap region for a given codim on the grid view.  */
  int overlapSize (int codim) const
  {
    return impl().overlapSize(codim);
  }

  /** \brief Return size of the ghost region for a given codim on the grid view.  */
  int ghostSize (int codim) const
  {
    return impl().ghostSize(codim);
  }

  /** \brief Communicate data on this view */
  template <class DataHandleImp, class DataType>
  auto communicate (CommDataHandleIF<DataHandleImp, DataType>& data,
                    InterfaceType iftype,
                    CommunicationDirection dir) const
  {
    using CommFuture = decltype(impl().communicate(data,iftype,dir));
    return communicate(data, iftype, dir, std::is_same<CommFuture, void>{});
  }

  /**
   * \brief access to the underlying implementation
   *
   * \warning Implementation details may change without prior notification.
   **/
  GridView& impl () { return gridView_; }

  /**
   * \brief access to the underlying implementation
   *
   * \warning Implementation details may change without prior notification.
   **/
  GridView const& impl () const { return gridView_; }

private:
  GridView gridView_;
  std::shared_ptr<IndexSet> indexSet_;
};

}

#endif // DUNE_SUBDOMAINS_SURFACEGRIDVIEW_HH
