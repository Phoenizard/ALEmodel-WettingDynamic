#ifndef DUNE_SUBDOMAININDEXSET_HH
#define DUNE_SUBDOMAININDEXSET_HH

#include <array>
#include <cassert>
#include <map>
#include <vector>

#include <dune/geometry/referenceelement.hh>
#include <dune/geometry/type.hh>

namespace Dune::Subdomains {

template <class IS>
class SubdomainIndexSet
{
  using IndexSet = IS;

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
  SubdomainIndexSet (IndexSet const& indexSet, std::vector<int> const& subdomains)
    : indexSet_(indexSet)
    , subdomains_(std::move(subdomains))
  {
    update();
  }

  SubdomainIndexSet (IndexSet const& indexSet)
    : indexSet_(indexSet)
    , subdomains_(indexSet.size(0), 0)
  {
    update();
  }

  template <class GridView>
  void update (GridView const& gridView)
  {
    update();
    std::array<IndexType,dimension> indices{};
    for (auto const& e : elements(gridView)) {
      auto& subIndexMap = indexMap_[subdomains_[impl().index(e)]];
      for (int d = 0; d < dimension; ++d)
        subIndexMap[d].resize(impl().size(dimension-d), invalid);

      for (int d = 0; d < dimension; ++d) {
        const int codim = dimension-d;
        auto& index = indices[d];

        auto refElem = referenceElement(e);
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
      //subdomains_.resize(impl().size(0), 0);
      for (int d = 0; d < dimension; ++d) {
          for (auto& entry : indexMap_)
              entry.second[d].clear();

          for (auto const& t : impl().types(dimension-d))
              sizes_[t] = 0;
      }
  }

  template <class Entity>
  void setSubdomain (Entity const& entity, int subdomain)
  {
    subdomains_[impl().index(entity)] = subdomain;
  }

  template <class Entity>
  int subdomain (Entity const& entity) const
  {
    return subdomains_[impl().index(entity)];
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

    IndexType const index_e = impl().template index<cc>(e);
    IndexType const index_s = impl().template subIndex<cc>(e,i,codim);

    int p = subdomains_[index_e];
    std::size_t d = dimension - codim;

    assert(indexMap_.at(p)[d][index_s] != invalid);
    return indexMap_.at(p)[d][index_s];
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
    return impl().contains(e);
  }

  IndexSet const& impl () const { return indexSet_; }

private:
  IndexSet const& indexSet_;
  std::vector<int> const& subdomains_;
  std::map<int, std::array<std::vector<IndexType>,dimension>> indexMap_;
  std::map<GeometryType, IndexType> sizes_;
};

} // end namespace Dune::Subdomains

#endif // DUNE_SUBDOMAININDEXSET_HH
