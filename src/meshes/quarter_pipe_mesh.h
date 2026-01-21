// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2026 Matthias Heil and Andrew Hazel
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
// Include guards
#ifndef OOMPH_QUARTER_PIPE_MESH_TEMPLATE_HEADER
#define OOMPH_QUARTER_PIPE_MESH_TEMPLATE_HEADER

#include "generic/element_factory.h"

// Generic oomph-lib includes
#include "generic/mesh.h"
#include "generic/brick_mesh.h"
#include "generic/refineable_brick_mesh.h"
#include "simple_cubic_mesh.h"
#include "quarter_pipe_domain.h"

#include "generic/macro_element.h"
#include "generic/domain.h"

namespace oomph
{
  //================================================================
  /// Non refineable quarter pipe mesh class
  /// Deform a simple cubic mesh into a quarter pipe
  /// r: radial direction
  /// theta: azimuthal direction
  /// z: axis direction
  //================================================================
  class QuarterPipeTFMesh : public virtual SimpleCubicTFMesh
  {
  public:
    /// Constructor: Pass number of elements in various directions,
    /// the inner and outer radius and the length of the tube
    QuarterPipeTFMesh(
      const unsigned& ntheta,
      const unsigned& nr,
      const unsigned& nz,
      const double& rmin,
      const double& rmax,
      const double& length,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// Empty Destructor
    virtual ~QuarterPipeTFMesh()
    {
      delete Domain_pt;
    }

    /// Access function to domain
    QuarterPipeDomain* domain_pt()
    {
      return Domain_pt;
    }

    /// Access function to underlying domain
    QuarterPipeDomain* domain_pt() const
    {
      return Domain_pt;
    }

  protected:
    /// Number of elements azimuthal direction
    unsigned Ntheta;

    /// Number of elements radial direction
    unsigned Nr;

    /// Number of elements axial direction
    unsigned Nz;

    /// Inner radius
    double Rmin;

    /// Outer radius
    double Rmax;

    /// Length
    double Length;

    /// Pointer to domain
    QuarterPipeDomain* Domain_pt;

  }; // endofclass


  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  //================================================================
  /// Refineable quarter pipe mesh class
  //================================================================
  class RefineableQuarterPipeTFMesh : public virtual QuarterPipeTFMesh,
                                      public RefineableBrickTFMesh
  {
  public:
    /// Constructor: Pass number of elements in various directions,
    /// the inner and outer radius and the length of the tube
    RefineableQuarterPipeTFMesh(
      const unsigned& ntheta,
      const unsigned& nr,
      const unsigned& nz,
      const double& rmin,
      const double& rmax,
      const double& length,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : SimpleCubicTFMesh(ntheta, nr, nz, 1.0, 1.0, 1.0, time_stepper_pt),
        QuarterPipeTFMesh(ntheta, nr, nz, rmin, rmax, length, time_stepper_pt)
    {
      // Setup Octree forest: Turn elements into individual octrees
      // and plant in forest
      Vector<TreeRoot*> trees_pt;
      for (unsigned iel = 0; iel < (nr * ntheta * nz); iel++)
      {
        OcTreeRoot* octree_root_pt = new OcTreeRoot(
          dynamic_cast<RefineableElement*>(this->element_pt(iel)));
        trees_pt.push_back(octree_root_pt);
      }

      this->Forest_pt = new OcTreeForest(trees_pt);
    }

    /// Destructor -- delete forest
    virtual ~RefineableQuarterPipeTFMesh()
    {
      delete this->Forest_pt;
    }

  }; // endofclass


  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  //================================================================
  /// Non refineable elastic quarter pipe mesh class
  /// setup lagrangian coordinates for solid mechanics problems
  //================================================================
  class ElasticQuarterPipeTFMesh : public virtual QuarterPipeTFMesh,
                                   public virtual SolidTFMesh
  {
  public:
    /// Constructor:  Pass number of elements in various directions,
    /// the inner and outer radius and the length of the tube.
    /// Builds mesh and copies Eulerian coords to Lagrangian
    /// ones so that the initial configuration is the stress-free one.
    ElasticQuarterPipeTFMesh(
      const unsigned& ntheta,
      const unsigned& nr,
      const unsigned& nz,
      const double& rmin,
      const double& rmax,
      const double& length,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : SimpleCubicTFMesh(ntheta, nr, nz, 1.0, 1.0, 1.0, time_stepper_pt),
        QuarterPipeTFMesh(ntheta, nr, nz, rmin, rmax, length, time_stepper_pt)
    {
      /// Make the current configuration the undeformed one by
      /// setting the nodal Lagrangian coordinates to their current
      /// Eulerian ones
      set_lagrangian_nodal_coordinates();
    }
  };


  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  //================================================================
  /// Refineable elastic quarter pipe mesh class
  //================================================================
  class ElasticRefineableQuarterPipeTFMesh
    : public virtual ElasticQuarterPipeTFMesh,
      public RefineableBrickTFMesh
  {
  public:
    /// Constructor:  Pass number of elements in various directions,
    /// the inner and outer radius and the length of the tube.
    /// Builds mesh and copies Eulerian coords to Lagrangian
    /// ones so that the initial configuration is the stress-free one.
    ElasticRefineableQuarterPipeTFMesh(
      const unsigned& ntheta,
      const unsigned& nr,
      const unsigned& nz,
      const double& rmin,
      const double& rmax,
      const double& length,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : SimpleCubicTFMesh(ntheta, nr, nz, 1.0, 1.0, 1.0, time_stepper_pt),
        QuarterPipeTFMesh(ntheta, nr, nz, rmin, rmax, length, time_stepper_pt),
        ElasticQuarterPipeTFMesh(
          ntheta, nr, nz, rmin, rmax, length, time_stepper_pt)
    {
      // Setup Octree forest: Turn elements into individual octrees
      // and plant in forest
      Vector<TreeRoot*> trees_pt;
      for (unsigned iel = 0; iel < (nr * ntheta * nz); iel++)
      {
        OcTreeRoot* octree_root_pt = new OcTreeRoot(
          dynamic_cast<RefineableElement*>(this->element_pt(iel)));
        trees_pt.push_back(octree_root_pt);
      }
      this->Forest_pt = new OcTreeForest(trees_pt);

      // Loop over all elements and set the undeformed macro element pointer
      unsigned n_element = this->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        // Get pointer to full element type
        auto* el_pt =
          dynamic_cast<RefineableSolidElement*>(this->element_pt(e));

        // Set pointer to macro element so the curvlinear boundaries
        // of the undeformed mesh/domain get picked up during adaptive
        // mesh refinement
        el_pt->set_undeformed_macro_elem_pt(
          this->Domain_pt->macro_element_pt(e));

        // Use MacroElement representation for
        // Lagrangian coordinates of newly created
        // nodes
        el_pt
          ->enable_use_of_undeformed_macro_element_for_new_lagrangian_coords();
      }
    }
  };

} // namespace oomph

#endif
