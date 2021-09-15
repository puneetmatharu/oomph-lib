// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//    Version 1.0; svn revision $LastChangedRevision$
// LIC//
// LIC// $LastChangedDate$
// LIC//
// LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================

#ifndef OOMPH_TETGEN_MESH_HEADER
#define OOMPH_TETGEN_MESH_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#ifdef OOMPH_HAS_MPI
// mpi headers
#include "mpi.h"
#endif

#include "generic/tetgen_scaffold_mesh.h"
#include "generic/tet_mesh.h"

namespace oomph
{
  //=========start of TetgenMesh class======================================
  /// \short  Unstructured tet mesh based on output from Tetgen:
  /// http://wias-berlin.de/software/tetgen/
  //========================================================================
  template<class ELEMENT>
  class TetgenMesh : public virtual TetMeshBase
  {
  public:
    /// \short Empty constructor
    TetgenMesh()
    {
      // Mesh can only be built with 3D Telements.
      MeshChecker::assert_geometric_element<TElementGeometricBase, ELEMENT>(3);
    }

    /// \short Constructor with the input files
    TetgenMesh(const std::string& node_file_name,
               const std::string& element_file_name,
               const std::string& face_file_name,
               TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper,
               const bool& use_attributes = false)
      : Tetgenio_exists(false), Tetgenio_pt(0)
    {
      // Mesh can only be built with 3D Telements.
      MeshChecker::assert_geometric_element<TElementGeometricBase, ELEMENT>(3);

      // Store the attributes
      Use_attributes = use_attributes;

      // Store timestepper used to build elements
      Time_stepper_pt = time_stepper_pt;

      // Build scaffold
      Tmp_mesh_pt = new TetgenScaffoldMesh(
        node_file_name, element_file_name, face_file_name);

      // Convert mesh from scaffold to actual mesh
      build_from_scaffold(time_stepper_pt, use_attributes);

      // Kill the scaffold
      delete Tmp_mesh_pt;
      Tmp_mesh_pt = 0;

      // Setup boundary coordinates
      unsigned nb = nboundary();
      for (unsigned b = 0; b < nb; b++)
      {
        bool switch_normal = false;
        setup_boundary_coordinates<ELEMENT>(b, switch_normal);
      }
    }

    /// \short Constructor with tetgenio data structure
    TetgenMesh(tetgenio& tetgen_data,
               TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper,
               const bool& use_attributes = false)

    {
      // Mesh can only be built with 3D Telements.
      MeshChecker::assert_geometric_element<TElementGeometricBase, ELEMENT>(3);

      // Store the attributes
      Use_attributes = use_attributes;

      // Store timestepper used to build elements
      Time_stepper_pt = time_stepper_pt;

      // We do not have a tetgenio representation
      Tetgenio_exists = false;
      Tetgenio_pt = 0;

      // Build scaffold
      Tmp_mesh_pt = new TetgenScaffoldMesh(tetgen_data);

      // Convert mesh from scaffold to actual mesh
      build_from_scaffold(time_stepper_pt, use_attributes);

      // Kill the scaffold
      delete Tmp_mesh_pt;
      Tmp_mesh_pt = 0;

      // Setup boundary coordinates
      unsigned nb = nboundary();
      for (unsigned b = 0; b < nb; b++)
      {
        bool switch_normal = false;
        setup_boundary_coordinates<ELEMENT>(b, switch_normal);
      }
    }

    /// \short Constructor with the input files. Setting the boolean
    /// flag to true splits "corner" elements, i.e. elements that
    /// that have at least three faces on a domain boundary. The
    /// relevant elements are split without introducing hanging
    /// nodes so the sons have a "worse" shape than their fathers.
    /// However, this step avoids otherwise-hard-to-diagnose
    /// problems in fluids problems where the application of
    /// boundary conditions at such "corner" elements can
    /// overconstrain the solution.
    TetgenMesh(const std::string& node_file_name,
               const std::string& element_file_name,
               const std::string& face_file_name,
               const bool& split_corner_elements,
               TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper,
               const bool& use_attributes = false)

    {
      // Mesh can only be built with 3D Telements.
      MeshChecker::assert_geometric_element<TElementGeometricBase, ELEMENT>(3);

      // Store the attributes
      Use_attributes = use_attributes;

      // Store timestepper used to build elements
      Time_stepper_pt = time_stepper_pt;

      // We do not have a tetgenio representation
      this->Tetgenio_exists = false;
      this->Tetgenio_pt = 0;

      // Build scaffold
      Tmp_mesh_pt = new TetgenScaffoldMesh(
        node_file_name, element_file_name, face_file_name);

      // Convert mesh from scaffold to actual mesh
      build_from_scaffold(time_stepper_pt, use_attributes);

      // Kill the scaffold
      delete Tmp_mesh_pt;
      Tmp_mesh_pt = 0;

      // Split corner elements
      if (split_corner_elements)
      {
        split_elements_in_corners<ELEMENT>();
      }

      // Setup boundary coordinates
      unsigned nb = nboundary();
      for (unsigned b = 0; b < nb; b++)
      {
        bool switch_normal = false;
        setup_boundary_coordinates<ELEMENT>(b, switch_normal);
      }
    }

    /// \short Constructor with tetgen data structure Setting the boolean
    /// flag to true splits "corner" elements, i.e. elements that
    /// that have at least three faces on a domain boundary. The
    /// relevant elements are split without introducing hanging
    /// nodes so the sons have a "worse" shape than their fathers.
    /// However, this step avoids otherwise-hard-to-diagnose
    /// problems in fluids problems where the application of
    /// boundary conditions at such "corner" elements can
    /// overconstrain the solution.
    TetgenMesh(tetgenio& tetgen_data,
               const bool& split_corner_elements,
               TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper,
               const bool& use_attributes = false)

    {
      // Mesh can only be built with 3D Telements.
      MeshChecker::assert_geometric_element<TElementGeometricBase, ELEMENT>(3);

      // Store the attributes
      Use_attributes = use_attributes;

      // Store timestepper used to build elements
      Time_stepper_pt = time_stepper_pt;

      // We do not have a tetgenio representation
      this->Tetgenio_exists = false;
      ;
      this->Tetgenio_pt = 0;

      // Build scaffold
      Tmp_mesh_pt = new TetgenScaffoldMesh(tetgen_data);

      // Convert mesh from scaffold to actual mesh
      build_from_scaffold(time_stepper_pt, use_attributes);

      // Kill the scaffold
      delete Tmp_mesh_pt;
      Tmp_mesh_pt = 0;

      // Split corner elements
      if (split_corner_elements)
      {
        split_elements_in_corners<ELEMENT>();
      }

      // Setup boundary coordinates
      unsigned nb = nboundary();
      for (unsigned b = 0; b < nb; b++)
      {
        bool switch_normal = false;
        setup_boundary_coordinates<ELEMENT>(b, switch_normal);
      }
    }

    /// \short Build mesh, based on a TetgenMeshFactedClosedSurface that
    /// specifies the outer boundary of the domain and any number of internal
    /// boundaries, specified by TetMeshFacetedSurfaces.
    /// Also specify target size for uniform element size.
    TetgenMesh(TetMeshFacetedClosedSurface* const& outer_boundary_pt,
               Vector<TetMeshFacetedSurface*>& internal_surface_pt,
               const double& element_volume,
               TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper,
               const bool& use_attributes = false,
               const bool& split_corner_elements = false)
    {
      // Mesh can only be built with 3D Telements.
      MeshChecker::assert_geometric_element<TElementGeometricBase, ELEMENT>(3);

      // Store the attributes
      Use_attributes = use_attributes;

      // Store timestepper used to build elements
      Time_stepper_pt = time_stepper_pt;

      // Copy across
      Outer_boundary_pt = outer_boundary_pt;

      // Setup reverse lookup scheme
      {
        unsigned n_facet = Outer_boundary_pt->nfacet();
        for (unsigned f = 0; f < n_facet; f++)
        {
          unsigned b = Outer_boundary_pt->one_based_facet_boundary_id(f);
          if (b != 0)
          {
            Tet_mesh_faceted_surface_pt[b - 1] = Outer_boundary_pt;
            Tet_mesh_facet_pt[b - 1] = Outer_boundary_pt->facet_pt(f);
          }
          else
          {
            std::ostringstream error_message;
            error_message << "Boundary IDs have to be one-based. Yours is " << b
                          << "\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }

      // Store the internal boundary
      Internal_surface_pt = internal_surface_pt;

      // Setup reverse lookup scheme
      {
        unsigned n = Internal_surface_pt.size();
        for (unsigned i = 0; i < n; i++)
        {
          unsigned n_facet = Internal_surface_pt[i]->nfacet();
          for (unsigned f = 0; f < n_facet; f++)
          {
            unsigned b = Internal_surface_pt[i]->one_based_facet_boundary_id(f);
            if (b != 0)
            {
              Tet_mesh_faceted_surface_pt[b - 1] = Internal_surface_pt[i];
              Tet_mesh_facet_pt[b - 1] = Internal_surface_pt[i]->facet_pt(f);
            }
            else
            {
              std::ostringstream error_message;
              error_message << "Boundary IDs have to be one-based. Yours is "
                            << b << "\n";
              throw OomphLibError(error_message.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
          }
        }
      }

      // Tetgen data structure for the input and output
      tetgenio in;
      this->build_tetgenio(outer_boundary_pt, internal_surface_pt, in);

      // Now tetrahedralise
      std::stringstream input_string;
      input_string << "pAa" << element_volume; // Andrew's parameters

      // input_string << "Vpq1.414Aa" << element_volume;
      // << "Q"; // Q for quiet!
      // << "V"; // V for verbose incl. quality output!

      // If any of the boundaries should not be split add the "Y" flag
      bool can_boundaries_be_split =
        outer_boundary_pt->boundaries_can_be_split_in_tetgen();
      {
        unsigned n_internal = internal_surface_pt.size();
        for (unsigned i = 0; i < n_internal; i++)
        {
          can_boundaries_be_split &=
            internal_surface_pt[i]->boundaries_can_be_split_in_tetgen();
        }
      }

      // If we can't split the boundaries add the flag
      if (can_boundaries_be_split == false)
      {
        input_string << "Y";
      }

      // Now convert to a C-style string
      char tetswitches[100];
      sprintf(tetswitches, "%s", input_string.str().c_str());

      // Make a new tetgen representation
      this->Tetgenio_exists = true;
      Tetgenio_pt = new tetgenio;
      tetrahedralize(tetswitches, &in, this->Tetgenio_pt);

      // Build scaffold
      Tmp_mesh_pt = new TetgenScaffoldMesh(*this->Tetgenio_pt);

      // If any of the objects are different regions then we need to use
      // the atributes
      bool regions_exist = false;
      {
        unsigned n_internal = internal_surface_pt.size();
        for (unsigned i = 0; i < n_internal; i++)
        {
          TetMeshFacetedClosedSurface* srf_pt =
            dynamic_cast<TetMeshFacetedClosedSurface*>(internal_surface_pt[i]);
          if (srf_pt != 0)
          {
            unsigned n_int_pts = srf_pt->ninternal_point_for_tetgen();
            for (unsigned j = 0; j < n_int_pts; j++)
            {
              regions_exist |=
                srf_pt->internal_point_identifies_region_for_tetgen(j);
            }
          }
        }
      }

      // If there are regions, use the attributes
      if (regions_exist)
      {
        Use_attributes = true;
      }

      // Convert mesh from scaffold to actual mesh
      build_from_scaffold(time_stepper_pt, Use_attributes);

      // Kill the scaffold
      delete Tmp_mesh_pt;
      Tmp_mesh_pt = 0;

      // Split corner elements
      if (split_corner_elements)
      {
        split_elements_in_corners<ELEMENT>();
      }

      // Setup boundary coordinates
      unsigned nb = nboundary();
      for (unsigned b = 0; b < nb; b++)
      {
        bool switch_normal = false;
        setup_boundary_coordinates<ELEMENT>(b, switch_normal);
      }

      // Now snap onto geometric objects associated with triangular facets
      // (if any!)
      snap_nodes_onto_geometric_objects();
    }

    ///\short Build tetgenio object from the TetMeshFacetedSurfaces
    void build_tetgenio(TetMeshFacetedSurface* const& outer_boundary_pt,
                        Vector<TetMeshFacetedSurface*>& internal_surface_pt,
                        tetgenio& tetgen_io)
    {
      // Pointer to Tetgen facet
      tetgenio::facet* f;
      // Pointer to Tetgen polygon
      tetgenio::polygon* p;

      // Start all indices from 1 (it's a choice and we've made it
      tetgen_io.firstnumber = 1;
      /// ALH: This may not be needed
      tetgen_io.useindex = true;

      // Find the number of internal surfaces
      const unsigned n_internal = internal_surface_pt.size();

      // Find the number of points on the outer boundary
      const unsigned n_outer_vertex = outer_boundary_pt->nvertex();
      tetgen_io.numberofpoints = n_outer_vertex;

      // Find the number of points on the inner boundaries and add to the totals
      Vector<unsigned> n_internal_vertex(n_internal);
      Vector<unsigned> internal_vertex_offset(n_internal);
      for (unsigned h = 0; h < n_internal; ++h)
      {
        n_internal_vertex[h] = internal_surface_pt[h]->nvertex();
        internal_vertex_offset[h] = tetgen_io.numberofpoints;
        tetgen_io.numberofpoints += n_internal_vertex[h];
      }

      // Read the data into the point list
      tetgen_io.pointlist = new double[tetgen_io.numberofpoints * 3];
      tetgen_io.pointmarkerlist = new int[tetgen_io.numberofpoints];
      unsigned counter = 0;
      for (unsigned n = 0; n < n_outer_vertex; ++n)
      {
        for (unsigned i = 0; i < 3; ++i)
        {
          tetgen_io.pointlist[counter] =
            outer_boundary_pt->vertex_coordinate(n, i);
          ++counter;
        }
      }
      for (unsigned h = 0; h < n_internal; ++h)
      {
        const unsigned n_inner = n_internal_vertex[h];
        for (unsigned n = 0; n < n_inner; ++n)
        {
          for (unsigned i = 0; i < 3; ++i)
          {
            tetgen_io.pointlist[counter] =
              internal_surface_pt[h]->vertex_coordinate(n, i);
            ++counter;
          }
        }
      }

      // Set up the pointmarkers
      counter = 0;
      for (unsigned n = 0; n < n_outer_vertex; ++n)
      {
        tetgen_io.pointmarkerlist[counter] =
          outer_boundary_pt->one_based_vertex_boundary_id(n);
        ++counter;
      }
      for (unsigned h = 0; h < n_internal; ++h)
      {
        const unsigned n_inner = n_internal_vertex[h];
        for (unsigned n = 0; n < n_inner; ++n)
        {
          tetgen_io.pointmarkerlist[counter] =
            internal_surface_pt[h]->one_based_vertex_boundary_id(n);
          ++counter;
        }
      }

      // Now the facets
      unsigned n_outer_facet = outer_boundary_pt->nfacet();
      tetgen_io.numberoffacets = n_outer_facet;
      Vector<unsigned> n_inner_facet(n_internal);
      for (unsigned h = 0; h < n_internal; ++h)
      {
        n_inner_facet[h] = internal_surface_pt[h]->nfacet();
        tetgen_io.numberoffacets += n_inner_facet[h];
      }

      tetgen_io.facetlist = new tetgenio::facet[tetgen_io.numberoffacets];
      tetgen_io.facetmarkerlist = new int[tetgen_io.numberoffacets];

      counter = 0;
      for (unsigned n = 0; n < n_outer_facet; ++n)
      {
        // Set pointer to facet
        f = &tetgen_io.facetlist[counter];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];

        Vector<unsigned> facet = outer_boundary_pt->vertex_index_in_tetgen(n);

        p->numberofvertices = facet.size();
        p->vertexlist = new int[p->numberofvertices];
        for (int i = 0; i < p->numberofvertices; ++i)
        {
          // The offset here is because we have insisted on 1-based indexing
          p->vertexlist[i] = facet[i] + 1;
        }

        // Set up the boundary markers
        tetgen_io.facetmarkerlist[counter] =
          outer_boundary_pt->one_based_facet_boundary_id(n);
        // Increase the counter
        ++counter;
      }

      // Initialise the number of holes
      tetgen_io.numberofholes = 0;
      // and the number of regions
      tetgen_io.numberofregions = 0;

      // Loop over the internal stuff
      for (unsigned h = 0; h < n_internal; ++h)
      {
        for (unsigned n = 0; n < n_inner_facet[h]; ++n)
        {
          // Set pointer to facet
          f = &tetgen_io.facetlist[counter];
          f->numberofpolygons = 1;
          f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
          f->numberofholes = 0;
          f->holelist = NULL;
          p = &f->polygonlist[0];

          Vector<unsigned> facet =
            internal_surface_pt[h]->vertex_index_in_tetgen(n);

          p->numberofvertices = facet.size();
          p->vertexlist = new int[p->numberofvertices];
          for (int i = 0; i < p->numberofvertices; ++i)
          {
            // Add the 1-based and vertex offsets to get these number correct
            p->vertexlist[i] = facet[i] + internal_vertex_offset[h] + 1;
          }
          // Set up the boundary markers
          tetgen_io.facetmarkerlist[counter] =
            internal_surface_pt[h]->one_based_facet_boundary_id(n);
          ++counter;
        }

        // If it's a hole/region add it
        TetMeshFacetedClosedSurface* srf_pt =
          dynamic_cast<TetMeshFacetedClosedSurface*>(internal_surface_pt[h]);
        if (srf_pt != 0)
        {
          unsigned n_int_pts = srf_pt->ninternal_point_for_tetgen();
          for (unsigned j = 0; j < n_int_pts; j++)
          {
            if (srf_pt->internal_point_identifies_hole_for_tetgen(j))
            {
              ++tetgen_io.numberofholes;
            }
            // Otherwise it may be region
            else
            {
              if (srf_pt->internal_point_identifies_region_for_tetgen(j))
              {
                ++tetgen_io.numberofregions;
              }
            }
          }
        }
      }

      // Set storage for the holes
      tetgen_io.holelist = new double[3 * tetgen_io.numberofholes];

      // Loop over all the internal boundaries again
      counter = 0;
      for (unsigned h = 0; h < n_internal; ++h)
      {
        TetMeshFacetedClosedSurface* srf_pt =
          dynamic_cast<TetMeshFacetedClosedSurface*>(internal_surface_pt[h]);
        if (srf_pt != 0)
        {
          unsigned n_int_pts = srf_pt->ninternal_point_for_tetgen();
          for (unsigned j = 0; j < n_int_pts; j++)
          {
            if (srf_pt->internal_point_identifies_hole_for_tetgen(j))
            {
              for (unsigned i = 0; i < 3; ++i)
              {
                tetgen_io.holelist[counter] =
                  srf_pt->internal_point_for_tetgen(j, i);
                ++counter;
              }
            }
          }
        }
      }

      // Set storage for the regions
      tetgen_io.regionlist = new double[5 * tetgen_io.numberofregions];
      // Loop over all the internal boundaries again
      counter = 0;
      for (unsigned h = 0; h < n_internal; ++h)
      {
        TetMeshFacetedClosedSurface* srf_pt =
          dynamic_cast<TetMeshFacetedClosedSurface*>(internal_surface_pt[h]);
        if (srf_pt != 0)
        {
          unsigned n_int_pts = srf_pt->ninternal_point_for_tetgen();
          for (unsigned j = 0; j < n_int_pts; j++)
          {
            if (srf_pt->internal_point_identifies_region_for_tetgen(j))
            {
              for (unsigned i = 0; i < 3; ++i)
              {
                tetgen_io.regionlist[counter] =
                  srf_pt->internal_point_for_tetgen(j, i);
                ++counter;
              }
              tetgen_io.regionlist[counter] =
                static_cast<double>(srf_pt->region_id_for_tetgen(j));
              ++counter;
              // Area target
              tetgen_io.regionlist[counter] = 0.0;
              ++counter;
            }
          }
        }
      }
    }

    /// Empty destructor
    ~TetgenMesh()
    {
      if (Tetgenio_exists)
      {
        delete Tetgenio_pt;
      }
    }

    /// \short Overload set_mesh_level_time_stepper so that the stored
    /// time stepper now corresponds to the new timestepper
    void set_mesh_level_time_stepper(TimeStepper* const& time_stepper_pt,
                                     const bool& preserve_existing_data)
    {
      this->Time_stepper_pt = time_stepper_pt;
    }

    /// Boolen defining whether tetgenio object has been built or not
    bool tetgenio_exists() const
    {
      return Tetgenio_exists;
    }

    /// Access to the triangulateio representation of the mesh
    tetgenio*& tetgenio_pt()
    {
      return Tetgenio_pt;
    }

    /// Set the tetgen pointer by a deep copy
    void set_deep_copy_tetgenio_pt(tetgenio* const& tetgenio_pt)
    {
      // Delete the existing data
      if (Tetgenio_exists)
      {
        delete Tetgenio_pt;
      }
      this->Tetgenio_pt = new tetgenio;
      // Create a deep copy of tetgenio_pt and store the result in
      // Tetgenio_pt
      this->deep_copy_of_tetgenio(tetgenio_pt, this->Tetgenio_pt);
    }

    /// Transfer tetgenio data from the input to the output
    /// The output is assumed to have been constructed and "empty"
    void deep_copy_of_tetgenio(tetgenio* const& input_pt, tetgenio*& output_pt);

  protected:
    /// Build mesh from scaffold
    void build_from_scaffold(TimeStepper* time_stepper_pt,
                             const bool& use_attributes);

    /// Temporary scaffold mesh
    TetgenScaffoldMesh* Tmp_mesh_pt;

    /// \short Boolean to indicate whether a tetgenio representation of the
    /// mesh exists
    bool Tetgenio_exists;

    /// \short Tetgen representation of mesh
    tetgenio* Tetgenio_pt;

    /// \short Boolean flag to indicate whether to use attributes or not
    /// (required for multidomain meshes)
    bool Use_attributes;

  }; // end class

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

  //==============start_mesh=================================================
  /// Tetgen-based mesh upgraded to become a solid mesh. Automatically
  /// enumerates all boundaries.
  //=========================================================================
  template<class ELEMENT>
  class SolidTetgenMesh : public virtual TetgenMesh<ELEMENT>,
                          public virtual SolidMesh
  {
  public:
    /// \short Constructor. Boundary coordinates are setup
    /// automatically.
    SolidTetgenMesh(const std::string& node_file_name,
                    const std::string& element_file_name,
                    const std::string& face_file_name,
                    const bool& split_corner_elements,
                    TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper,
                    const bool& use_attributes = false)
      : TetgenMesh<ELEMENT>(node_file_name,
                            element_file_name,
                            face_file_name,
                            split_corner_elements,
                            time_stepper_pt,
                            use_attributes)
    {
      // Assign the Lagrangian coordinates
      set_lagrangian_nodal_coordinates();
    }

    /// \short Constructor. Boundary coordinates are re-setup
    /// automatically, with the orientation of the outer unit
    /// normal determined by switch_normal.
    SolidTetgenMesh(const std::string& node_file_name,
                    const std::string& element_file_name,
                    const std::string& face_file_name,
                    const bool& split_corner_elements,
                    const bool& switch_normal,
                    TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper,
                    const bool& use_attributes = false)
      : TetgenMesh<ELEMENT>(node_file_name,
                            element_file_name,
                            face_file_name,
                            split_corner_elements,
                            time_stepper_pt,
                            use_attributes)
    {
      // Assign the Lagrangian coordinates
      set_lagrangian_nodal_coordinates();

      // Re-setup boundary coordinates for all boundaries with specified
      // orientation of nnormal
      unsigned nb = this->nboundary();
      for (unsigned b = 0; b < nb; b++)
      {
        this->template setup_boundary_coordinates<ELEMENT>(b, switch_normal);
      }
    }

    /// Empty Destructor
    virtual ~SolidTetgenMesh() {}
  };

} // namespace oomph

#include "tetgen_mesh.tpp"
#endif