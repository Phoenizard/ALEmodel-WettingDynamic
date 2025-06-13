#ifndef DUNE_SUBDOMAINS_SURFACEINDEXSET_HH
#define DUNE_SUBDOMAINS_SURFACEINDEXSET_HH

#include <array>
#include <cassert>
#include <map>
#include <vector>

#include <dune/geometry/referenceelement.hh>
#include <dune/geometry/type.hh>

namespace Dune::Subdomains {

template <class IS, class IS0>
class SurfaceIndexSet
{
    using IndexSet = IS;
    using IndexSet0 = IS0;

public:
  /** \brief Export the type of the entity used as parameter in the index(...) method */
  template <int cc>
  using Codim = typename IndexSet::template Codim<cc>;

  /** \brief The type used for the indices */
  using IndexType = typename IndexSet::IndexType;

  /** \brief iterator range for geometry types in domain */
  using Types = typename IndexSet::Types;

  enum {
    dimension = IndexSet::dimension //< dimension of the grid
  };

  static inline constexpr IndexType invalid = IndexType(-1);

public:
    SurfaceIndexSet (IndexSet const& indexSet, IndexSet0 const& indexSet0, std::vector<int> const& facets)
    : indexSet_(indexSet)
    , indexSet0_(indexSet0)
    //, facets_(std::move(facets))
    , facetsCoarse_(facets)
    , numNeighborfacets_(indexSet_.size(0),2)
  {
    update();
  }

  SurfaceIndexSet (IndexSet const& indexSet)
    : indexSet_(indexSet)
    , indexSet0_(indexSet)
    //, facets_(indexSet.size(0), 0)
    , facetsCoarse_(indexSet0_.size(0))
    , numNeighborfacets_(indexSet_.size(),2)
  {
    update();
  }

  template <class GridView>
  void update (GridView const& gridView)
  {
    if (!facetsCorrected)
        correctFacets(gridView);

    update();

    // get info if shell is closed from initfile (default: no)
    auto closedShell = AMDiS::Parameters::get<int>("closed shell").value_or(1);
    numNeighborfacets_.resize(impl().size(0),2); //default = 2 neighboring facets
    if (!closedShell) numNeighborfacets_ = numNeighborFacets(gridView); //to be called if and only if shell is not closed!

    std::array<IndexType,dimension> indices{};
    for (auto const& e : elements(gridView)) {
        auto father = e;
        while (father.hasFather())
            father = father.father();
        std::size_t idx = impl().index(e);
        std::size_t idx0 = impl0().index(father);

      int facet = facetsCoarse_[idx0]; //facets = 0 or = 1 if vertex on element touches shell, -1 else
      auto& subIndexMap = indexMap_[-1];
      for (int d = 0; d < dimension; ++d)
          subIndexMap[d].resize(impl().size(dimension-d), invalid);

      auto refElem = referenceElement(e);
      for (int d = 0; d < dimension; ++d) {
          if (facet == 1) {
              indexMap_[facet][d].resize(impl().size(dimension-d), invalid);
              for (auto const& is : intersections(gridView,e)) {
                  if (!is.neighbor())
                      continue;

                  auto fatherOut = is.outside();
                  while (fatherOut.hasFather())
                      fatherOut = fatherOut.father();

                  if (facetsCoarse_[impl0().index(fatherOut)] == 0) { //now I am on the desired intersection
                      int seIdx = is.indexInInside();
                      const int codim = dimension-d;
                      auto& index = indices[d];
                      auto subRefElem = referenceElement(e.template subEntity<1>(seIdx));
                      for (int i = 0; i < refElem.size(seIdx,1,codim); ++i) {
                          auto& newIndex = indexMap_[facet][d][impl().subIndex(e,refElem.subEntity(seIdx,1,i,codim),codim)];
                          if (newIndex == invalid &&
                                (numNeighborfacets_[indexSet_.index(e)] == 2 || codim < 2)) { // don't add extra indices on the end of surface
                              newIndex = index++;
                              sizes_[subRefElem.type(i,codim-1)]++;
                          }
                      }
                  }
              }
          }
          const int codim = dimension-d;
          auto& index = indices[d];
          for (int i = 0; i < refElem.size(codim); ++i) {
              auto& newIndex = subIndexMap[d][impl().subIndex(e,i,codim)];
              if (newIndex == invalid) {
                  newIndex = index++;
                  sizes_[refElem.type(i,codim)]++;
              }
          }
      }
    }
  }

  void update ()
  {
    //facets_.resize(impl().size(0), -1);
    for (int d = 0; d < dimension; ++d) {
        for (auto& entry : indexMap_)
            entry.second[d].clear();

        for (auto const& t : impl().types(dimension-d))
            sizes_[t] = 0;
    }
  }

  // correct facetsCoarse_ s.t. elements, where only one vertex touches the interface are also labeled with 1
  template<class GV>
  void correctFacets(GV const& gridView) { // todo: make correct for open grids if more than three elements on one side touch a shell point
        std::vector<int> indexes;
        for (auto const& e : elements(gridView.grid().levelGridView(0))) {
            int p = facetsCoarse_[impl0().index(e)];
            for (auto const& is : intersections(gridView.grid().levelGridView(0),e)) {
                if (!is.neighbor())
                    continue;

                int q = facetsCoarse_[impl0().index(is.outside())];
                if (p == -1 && q == 1)
                    indexes.push_back(impl0().index(e));
            }
        }
        // make facets = 1 also for elements, where only a vertex touches the shell
        for (int i = 0; i < indexes.size(); i++) {
            facetsCoarse_[indexes[i]] = 1;
        }
        facetsCorrected = true;
    }


  // find number of neighboring facets of each interface facet
  // loop over all elements. for each element, find if there is an intersection of the shell
  // if so, loop again over all elements to find neighboring intersections
  // Remark: done with facetsCoarse_: a vector filled with 1 (left) and 0 (right) for elements
  // of the unrefined base grid touching the membrane, -1 else
  template<class GridView>
  std::vector<int> numNeighborFacets(GridView const& gridView) {
        std::vector<int> neighbors(impl().size(0),0);
        for (const auto& e : elements(gridView)) {
            auto father = e;
            while (father.hasFather())
                father = father.father();
            std::size_t idx = impl().index(e);
            std::size_t idx0 = impl0().index(father);
            if (facetsCoarse_[idx0] != 1) //only loop through elements touching the shell (from one side)
                continue;

            int isBoundaryEl = 0;
            for (const auto& is : intersections(gridView,e)) {
                if (is.boundary())  //if the element, where the facet belongs to, touches the domain boundary, set neighbors to 2
                    isBoundaryEl = 1;

                if (!is.neighbor())
                    continue;
                auto fatherOut = is.outside();
                while (fatherOut.hasFather())
                    fatherOut = fatherOut.father();
                if (facetsCoarse_[impl0().index(fatherOut)] == 0) { //now I am on the desired Intersection and try to find neighbor Intersections
                    for (const auto& e2 : elements(gridView)) {
                        auto father2 = e2;
                        while (father2.hasFather())
                            father2 = father2.father();
                        std::size_t idx2 = impl().index(e2);
                        std::size_t idx2_0 = impl0().index(father2);

                        if (facetsCoarse_[idx2_0] != 1 || idx2 == idx) // check only elements with coarse index 1 and different from e
                            continue;

                        for (const auto& is2 : intersections(gridView,e2)) {
                            if (!is2.neighbor())
                                continue;
                            auto fatherOut2 = is2.outside();
                            while (fatherOut2.hasFather())
                                fatherOut2 = fatherOut2.father();

                            if (facetsCoarse_[impl0().index(fatherOut2)] == 0) { // candidate for the neighbor
                                for (int i = 0; i < is.geometry().corners(); i++) {
                                    for (int j = 0; j < is2.geometry().corners(); j++) {
                                        if ((is.geometry().corner(i)-is2.geometry().corner(j)).two_norm() < 1e-5)
                                            neighbors[idx]++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (isBoundaryEl && neighbors[idx] > 0)
                neighbors[idx] = 2;
        }
        return neighbors;
    }

  //todo: make usable for refinement
  template <class Entity>
  void setSurface (Entity const& entity, int facet)
  {
    auto father = entity;
    while(father.hasFather())
        father = father.father();
    facetsCoarse_[impl0().index(father)] = facet;
  }

  //todo: make usable for refinement
  template <class Entity>
  int surface (Entity const& entity) const
  {
    auto father = entity;
    while(father.hasFather())
        father = father.father();
    return facetsCoarse_[impl().index(entity)];
  }

public:
  /** \brief Map entity to index. The result of calling this method with an entity that is not
   * in the index set is undefined.
   *
   * \param e  Reference to codim cc entity, where cc is the template parameter of the function.
   * \return   An index in the range 0 ... Max number of entities in set - 1.
   */
  template <int cc>
  IndexType index (typename Codim<cc>::Entity const& e) const
  {
    static_assert(cc == 0);
    return impl().template index<cc>(e);
  }

  /** \brief Map entity to index. Easier to use than the above because codimension template
   * parameter need not be supplied explicitly.
   * The result of calling this method with an entity that is not in the index set is undefined.
   *
   * \param e  Reference to codim cc entity. Since entity knows its codimension, automatic
   *           extraction is possible.
   * \return   An index in the range 0 ... Max number of entities in set - 1.
   */
  template <class Entity>
  IndexType index (Entity const& e) const
  {
    enum { cc = Entity::codimension };
    static_assert(cc == 0);
    return impl().template index<cc>(e);
  }

  /** \brief Map a subentity to an index.
   *
   *  The result of calling this method with an entity that is not in the
   *  index set is undefined.
   *
   *  \tparam  cc  codimension of the entity
   *
   *  \param[in]  e      reference to codimension cc entity
   *  \param[in]  i      number subentity of e within the codimension
   *  \param[in]  codim  codimension of the subentity we're interested in
   *
   *  \note The parameter <tt>codim</tt> denotes the codimension with respect
   *        to the grid, i.e., it must satisfy cc <= codim <= dimension.
   *
   *  \return An index in the range 0 ... Max number of entities in set - 1.
   */
  template <int cc>
  IndexType subIndex (typename Codim<cc>::Entity const& e, int i, unsigned int codim) const
  {
      static_assert(cc == 0);
      if (codim == 0)
          return impl().template index<cc>(e);

      auto father = e;
      while (father.hasFather())
          father = father.father();
      IndexType const index_e = impl().template index<cc>(e);
      IndexType const index_e0 = impl0().template index<cc>(father);
      IndexType const index_s = impl().template subIndex<cc>(e,i,codim);

      int p = facetsCoarse_[index_e0];
      std::size_t d = dimension - codim;

      if (p == 1 && indexMap_.at(p)[d][index_s] != invalid) {
          return indexMap_.at(p)[d][index_s];
      }
      assert(indexMap_.at(-1)[d][index_s] != invalid);
      return indexMap_.at(-1)[d][index_s];
  }

  /** \brief Map a subentity to an index.
   *
   *  The result of calling this method with an entity that is not in the
   *  index set is undefined.
   *
   *  \note This method exists for convenience only.
   *        It extracts the codimension from the type of the entity, which can
   *        be guessed by the compiler.
   *
   *  \tparam  Entity  type of entity (must be GridImp::Codim< cc >::Entity
   *                   for some cc)
   *
   *  \param[in]  e      reference to entity
   *  \param[in]  i      number subentity of e within the codimension
   *  \param[in]  codim  codimension of the subentity we're interested in
   *
   *  \note The parameter <tt>codim</tt> denotes the codimension with respect
   *        to the grid, i.e., it must satisfy cc <= codim <= dimension.
   *
   *  \return An index in the range 0 ... Max number of entities in set - 1.
   */
  template <class Entity>
  IndexType subIndex (Entity const& e, int i, unsigned int codim) const
  {
      enum { cc = Entity::codimension };
      return subIndex<cc>(e,i,codim);
  }

  /**
   * \brief obtain all geometry types of entities in domain
   *
   * This method returns an iterator range (something that behaves like
   * Dune::IteratorRange) visiting all geometry types of codimension codim
   * in the domain of the index map exactly once.
   * The iterator must implement the concept of a forward iterator (in the
   * sense of the STL).
   * The elements in the iterator range are required to be of type
   * Dune::GeometryType.
   *
   * \param[in]  codim  a valid codimension
   *
   * \return iterator range over Const reference to a vector of geometry types.
   */
  Types types (int codim) const
  {
    return impl().types(codim);
  }

  /** \brief Return total number of entities of given geometry type in entity set \f$E\f$.
   *
   * \param[in] type A valid geometry type.
   * \return         number of entities.
   */
  IndexType size (GeometryType const& type) const
  {
    return sizes_.count(type) ? sizes_.at(type) : 0;
  }

  /** \brief Return total number of entities of given codim in the entity set \f$E\f$. This
   * is simply a sum over all geometry types.
   *
   * \param[in] codim A valid codimension
   * \return    number of entities.
   */
  IndexType size (int codim) const
  {
    IndexType s = 0;
    for (GeometryType const& t : impl().types(codim))
      s += this->size(t);
    return s;
  }

  /** \brief Return true if the given entity is contained in \f$E\f$.
   *
   * \note If the input element e is not an element of the grid, then
   *       the result of contains() is undefined.
   */
  template <class Entity>
  bool contains (Entity const& e) const
  {
    return impl().contains(e) && index(e) != invalid;
  }

  IndexSet const& impl () const { return indexSet_; }
  IndexSet0 const& impl0 () const { return indexSet0_; }


private:
  IndexSet const& indexSet_;
  IndexSet0 const& indexSet0_;
//  std::vector<int> const& facets_;
  std::vector<int> facetsCoarse_;
  std::vector<int> numNeighborfacets_;
  std::map<int, std::array<std::vector<IndexType>,dimension>> indexMap_;

  std::map<GeometryType, IndexType> sizes_;
  bool facetsCorrected = false;
};

} // end namespace Dune::Subdomains

#endif // DUNE_SUBDOMAINS_SURFACEINDEXSET_HH
