// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    codim0extractor.hh
 *  Version:     1.0
 *  Created on:  Jun 23, 2009
 *  Author:      Oliver Sander, Christian Engwer
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: base class for grid extractors extracting surface grids
 *
 */
/**
 * @file
 * @brief Mesh grid extractor base class
 */

#ifndef DUNE_GRIDGLUE_EXTRACTORS_CODIM0EXTRACTOR_HH
#define DUNE_GRIDGLUE_EXTRACTORS_CODIM0EXTRACTOR_HH

#include <deque>
#include <functional>

#include <dune/common/deprecated.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include "extractor.hh"

namespace Dune {

  namespace GridGlue {

/**
 * extractor for codim-0 entities (elements)
 */
template<typename GV>
class Codim0Extractor : public Extractor<GV,0>
{

public:

  /*  E X P O R T E D  T Y P E S   A N D   C O N S T A N T S  */
  using Extractor<GV,0>::codim;
  typedef typename Extractor<GV,0>::ctype ctype;
  using Extractor<GV,0>::dim;
  using Extractor<GV,0>::dimworld;
  typedef typename Extractor<GV,0>::IndexType IndexType;

  typedef typename GV::Traits::template Codim<dim>::Entity Vertex;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef std::function<bool(const Element&, unsigned int subentity)> Predicate;

  // import typedefs from base class
  typedef typename Extractor<GV,0>::SubEntityInfo SubEntityInfo;
  typedef typename Extractor<GV,0>::ElementInfo ElementInfo;
  typedef typename Extractor<GV,0>::VertexInfo VertexInfo;
  typedef typename Extractor<GV,0>::CoordinateInfo CoordinateInfo;
  typedef typename Extractor<GV,0>::VertexInfoMap VertexInfoMap;

  /**
   * @brief Constructor
   * @param gv the grid view object to work with
   * @param predicate a predicate to mark entities for extraction (unary functor returning bool)
   */
  Codim0Extractor(const GV& gv, const Predicate& predicate)
    :  Extractor<GV,0>(gv), positiveNormalDirection_(false)
  {
    std::cout << "This is Codim0Extractor on a <"
              << GV::dimension << "," << GV::dimensionworld << "> grid!" << std::endl;
    update(predicate);
  }

  bool & positiveNormalDirection() { return positiveNormalDirection_; }
  const bool & positiveNormalDirection() const { return positiveNormalDirection_; }

protected:
  bool positiveNormalDirection_;
private:
  void update(const Predicate& predicate);
};


template<typename GV>
void Codim0Extractor<GV>::update(const Predicate& predicate)
{
  // In this first pass iterate over all entities of codim 0.
  // Get its corner vertices, find resp. store them together with their associated index,
  // and remember the indices of the corners.

  // free everything there is in this object
  this->clear();

  // several counter for consecutive indexing are needed
  size_t element_index = 0;
  size_t vertex_index = 0;

  // a temporary container where newly acquired face
  // information can be stored at first
  std::deque<SubEntityInfo> temp_faces;

  //create a coordinate vector
  std::vector<FieldVector<double,dim>> coords(this->gv_.size(dim));
  // iterate over all codim 0 elements on the grid
  for (const auto& elmt : elements(this->gv_, Partitions::interior))
  {
    const auto geometry = elmt.geometry();
    IndexType eindex = this->cellMapper_.index(elmt);

    // only do sth. if this element is "interesting"
    // implicit cast is done automatically
    if (predicate(elmt,0))
    {
      // add an entry to the element info map, the index will be set properly later
      this->elmtInfo_.emplace(eindex, ElementInfo(element_index, elmt, 1));

      unsigned int numCorners = elmt.subEntities(dim);
      unsigned int vertex_indices[numCorners];       // index in global vector
      unsigned int vertex_numbers[numCorners];       // index in parent entity

      // try for each of the faces vertices whether it is already inserted or not
      for (unsigned int i = 0; i < numCorners; ++i)
      {
        vertex_numbers[i] = i;

        // get the vertex pointer and the index from the index set
        const Vertex vertex = elmt.template subEntity<dim>(vertex_numbers[i]);
        IndexType vindex = this->gv_.indexSet().template index<dim>(vertex);

        // if the vertex is not yet inserted in the vertex info map
        // it is a new one -> it will be inserted now!
        typename VertexInfoMap::iterator vimit = this->vtxInfo_.find(vindex);
        if (vimit == this->vtxInfo_.end())
        {
          // insert into the map
          this->vtxInfo_.emplace(vindex, VertexInfo(vertex_index, vertex));
          // remember this vertex' index
          vertex_indices[i] = vertex_index;
          // increase the current index
          vertex_index++;
          coords[vindex] = elmt.geometry().corner(vertex_numbers[i]);
        }
        else
        {
          // only remember the vertex' index
          vertex_indices[i] = vimit->second.idx;
        }
      }

      // flip cell if necessary
      {
        switch (int(dim))
        {
        case 0 :
          break;
        case 1 :
        {
          // The following test only works if the zero-th coordinate is the
          // one that defines the orientation.  A sufficient condition for
          // this is dimworld == 1
          /* assert(dimworld==1); */
          bool elementNormalDirection =
            (geometry.corner(1)[0] < geometry.corner(0)[0]);
          if ( positiveNormalDirection_ != elementNormalDirection )
          {
            std::swap(vertex_indices[0], vertex_indices[1]);
            std::swap(vertex_numbers[0], vertex_numbers[1]);
          }
          break;
        }
        case 2 :
        {
          Dune::FieldVector<ctype, dimworld>
          v0 = geometry.corner(1),
            v1 = geometry.corner(2);
          v0 -= geometry.corner(0);
          v1 -= geometry.corner(0);
          ctype normal_sign = v0[0]*v1[1] - v0[1]*v1[0];
          bool elementNormalDirection = (normal_sign < 0);
          if ( positiveNormalDirection_ != elementNormalDirection )
          {
            if (elmt.type().isCube())
            {
              for (int i = 0; i < (1<<dim); i+=2)
              {
                // swap i and i+1
                std::swap(vertex_indices[i], vertex_indices[i+1]);
                std::swap(vertex_numbers[i], vertex_numbers[i+1]);
              }
            } else if (elmt.type().isSimplex()) {
              std::swap(vertex_indices[0], vertex_indices[1]);
              std::swap(vertex_numbers[0], vertex_numbers[1]);
            } else {
              DUNE_THROW(Dune::Exception, "Unexpected Geometrytype");
            }
          }
          break;
        }
        }
      }

      // add a new face to the temporary collection
      temp_faces.emplace_back(eindex, 0, elmt.type());
      element_index++;
      for (unsigned int i=0; i<numCorners; i++) {
        temp_faces.back().corners[i].idx = vertex_indices[i];
        // remember the vertices' numbers in parent element's vertices
        temp_faces.back().corners[i].num = vertex_numbers[i];
      }

    }
  }   // end loop over elements

  // allocate the array for the face specific information...
  this->subEntities_.resize(element_index);
  // ...and fill in the data from the temporary containers
  copy(temp_faces.begin(), temp_faces.end(), this->subEntities_.begin());

  // now first write the array with the coordinates...
  this->coords_.resize(this->vtxInfo_.size());
  for (const auto& vinfo : this->vtxInfo_)
  {
    // get a pointer to the associated info object
    CoordinateInfo* current = &this->coords_[vinfo.second.idx];
    // store this coordinates index // NEEDED?
    current->index = vinfo.second.idx;
    // store the vertex' index for the index2vertex mapping
    current->vtxindex = vinfo.first;
    // store the vertex' coordinates under the associated index
    // in the coordinates array
    const auto vtx = this->grid().entity(vinfo.second.p);
    current->coord = coords[vinfo.first];
  }

}

}  // namespace GridGlue

}  // namespace Dune

#endif // DUNE_GRIDGLUE_EXTRACTORS_CODIM0EXTRACTOR_HH
