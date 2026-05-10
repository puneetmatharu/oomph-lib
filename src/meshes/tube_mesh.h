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
#ifndef OOMPH_TUBE_MESH_HEADER
#define OOMPH_TUBE_MESH_HEADER

#include "generic/element_factory.h"

// Headers
#include "generic/refineable_brick_mesh.h"
#include "generic/macro_element.h"
#include "generic/domain.h"
#include "generic/algebraic_elements.h"
#include "generic/brick_mesh.h"
#include "generic/macro_element_node_update_element.h"

// Include the headers file for domain
#include "tube_domain.h"

namespace oomph
{
  //====================================================================
  /// 3D tube mesh class.
  /// The domain is specified by the GeomObject that identifies
  /// the entire volume. Non-refineable base version!
  ///
  /// The mesh boundaries are numbered as follows:
  /// - Boundary 0: "Inflow" cross section; located along the
  ///               line parametrised by \f$ \xi_0 =  \xi_0^{lo} \f$.
  /// - Boundary 1: The outer wall, represetned by \f$\xi_2 = 1\f$.
  /// - Boundary 2: The out flow, represented by \f$\xi_0 = \xi_0^{hi}\f$.
  ///
  //====================================================================
  class TubeTFMesh : public virtual BrickMeshBase
  {
  public:
    /// Constructor: Pass pointer to geometric object that
    /// specifies the volume, start and end coordinates for the centreline
    /// on the geometric object. Values of theta at which dividing lines
    /// are to be placed, fractions of the radius for the central box
    /// at the dividing lines, the number of layers
    /// and the timestepper.
    /// Timestepper defaults to Steady dummy timestepper.
    TubeTFMesh(GeomObject* wall_pt,
               const Vector<double>& centreline_limits,
               const Vector<double>& theta_positions,
               const Vector<double>& radius_box,
               const unsigned& nlayer,
               TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// Constructor with runtime element factory.
    TubeTFMesh(const ElementFactoryPt& element_factory_pt,
               GeomObject* wall_pt,
               const Vector<double>& centreline_limits,
               const Vector<double>& theta_positions,
               const Vector<double>& radius_box,
               const unsigned& nlayer,
               TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
    {
      set_element_factory_pt(element_factory_pt);
      // Delegate to existing ctor body via placement-new style: simulate by
      // calling the same logic via a helper. We can't delegate to another
      // ctor here because of virtual base; just inline the work.
      build_tube_mesh(wall_pt,
                      centreline_limits,
                      theta_positions,
                      radius_box,
                      nlayer,
                      time_stepper_pt);
    }

    /// Destructor: empty
    virtual ~TubeTFMesh()
    {
      delete Domain_pt;
    }

  private:
    /// Helper that contains the body of the original constructor.
    void build_tube_mesh(GeomObject* wall_pt,
                         const Vector<double>& centreline_limits,
                         const Vector<double>& theta_positions,
                         const Vector<double>& radius_box,
                         const unsigned& nlayer,
                         TimeStepper* time_stepper_pt);

  public:
    /// Access function to GeomObject representing wall
    GeomObject*& volume_pt()
    {
      return Volume_pt;
    }

    /// Access function to domain
    TubeDomain* domain_pt()
    {
      return Domain_pt;
    }

    /// Access function to underlying domain
    TubeDomain* domain_pt() const
    {
      return Domain_pt;
    }

  protected:
    /// Pointer to domain
    TubeDomain* Domain_pt;

    /// Pointer to the geometric object that represents the curved wall
    GeomObject* Volume_pt;
  };


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  //=============================================================
  /// Adaptative version of the TubeTFMesh base mesh.
  /// The domain is specified by the GeomObject that identifies
  /// the entire volume
  ///
  /// The mesh boundaries are numbered as follows:
  /// - Boundary 0: "Inflow" cross section; located along the
  ///               line parametrised by \f$ \xi_0 =  \xi_0^{lo} \f$.
  /// - Boundary 1: The outer wall, represetned by \f$\xi_2 = 1\f$.
  /// - Boundary 2: The out flow, represented by \f$\xi_0 = \xi_0^{hi}\f$.
  ///
  //=============================================================
  class RefineableTubeTFMesh : public TubeTFMesh, public RefineableBrickTFMesh

  {
  public:
    /// Constructor for adaptive deformable quarter tube mesh class.
    /// Pass pointer to geometric object that
    /// specifies the volume, start and end coordinates for the centreline
    /// on the geometric object. Values of theta at which dividing lines
    /// are to be placed, fractions of the radius for the central box
    /// at the dividing lines, the number of layers
    /// and the timestepper.
    /// Timestepper defaults to Steady dummy timestepper.
    RefineableTubeTFMesh(
      GeomObject* wall_pt,
      const Vector<double>& centreline_limits,
      const Vector<double>& theta_positions,
      const Vector<double>& radius_box,
      const unsigned& nlayer,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : TubeTFMesh(wall_pt,
                   centreline_limits,
                   theta_positions,
                   radius_box,
                   nlayer,
                   time_stepper_pt)
    {
      // Loop over all elements and set macro element pointer
      for (unsigned ielem = 0; ielem < TubeTFMesh::nelement(); ielem++)
      {
        dynamic_cast<RefineableQElement<3>*>(this->element_pt(ielem))
          ->set_macro_elem_pt(this->Domain_pt->macro_element_pt(ielem));
      }

      // Setup Octree forest: Turn elements into individual octrees
      // and plant in forest
      Vector<TreeRoot*> trees_pt;
      for (unsigned iel = 0; iel < this->nelement(); iel++)
      {
        FiniteElement* el_pt = this->finite_element_pt(iel);
        OcTreeRoot* octree_root_pt =
          new OcTreeRoot(dynamic_cast<RefineableElement*>(el_pt));
        trees_pt.push_back(octree_root_pt);
      }
      this->Forest_pt = new OcTreeForest(trees_pt);

#ifdef PARANOID
      // Run self test
      unsigned success_flag =
        dynamic_cast<OcTreeForest*>(this->Forest_pt)->self_test();
      if (success_flag == 0)
      {
        oomph_info << "Successfully built octree forest " << std::endl;
      }
      else
      {
        throw OomphLibError("Trouble in building octree forest ",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }

    /// Constructor with runtime element factory.
    RefineableTubeTFMesh(
      const ElementFactoryPt& element_factory_pt,
      GeomObject* wall_pt,
      const Vector<double>& centreline_limits,
      const Vector<double>& theta_positions,
      const Vector<double>& radius_box,
      const unsigned& nlayer,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : TubeTFMesh(element_factory_pt,
                   wall_pt,
                   centreline_limits,
                   theta_positions,
                   radius_box,
                   nlayer,
                   time_stepper_pt)
    {
      for (unsigned ielem = 0; ielem < TubeTFMesh::nelement(); ielem++)
      {
        dynamic_cast<RefineableQElement<3>*>(this->element_pt(ielem))
          ->set_macro_elem_pt(this->Domain_pt->macro_element_pt(ielem));
      }
      Vector<TreeRoot*> trees_pt;
      for (unsigned iel = 0; iel < this->nelement(); iel++)
      {
        FiniteElement* el_pt = this->finite_element_pt(iel);
        OcTreeRoot* octree_root_pt =
          new OcTreeRoot(dynamic_cast<RefineableElement*>(el_pt));
        trees_pt.push_back(octree_root_pt);
      }
      this->Forest_pt = new OcTreeForest(trees_pt);
    }

    /// Destructor: empty
    virtual ~RefineableTubeTFMesh() {}
  };


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //=========================================================================
  /// Backwards-compatible templated TubeMesh wrappers.
  //=========================================================================
  template<class ELEMENT>
  class TubeMesh : public TubeTFMesh
  {
  public:
    TubeMesh(GeomObject* wall_pt,
             const Vector<double>& centreline_limits,
             const Vector<double>& theta_positions,
             const Vector<double>& radius_box,
             const unsigned& nlayer,
             TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : TubeTFMesh(make_element_factory_pt<ELEMENT>(),
                   wall_pt,
                   centreline_limits,
                   theta_positions,
                   radius_box,
                   nlayer,
                   time_stepper_pt)
    {
    }
  };


  template<class ELEMENT>
  class RefineableTubeMesh : public RefineableTubeTFMesh
  {
  public:
    RefineableTubeMesh(
      GeomObject* wall_pt,
      const Vector<double>& centreline_limits,
      const Vector<double>& theta_positions,
      const Vector<double>& radius_box,
      const unsigned& nlayer,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : RefineableTubeTFMesh(make_element_factory_pt<ELEMENT>(),
                             wall_pt,
                             centreline_limits,
                             theta_positions,
                             radius_box,
                             nlayer,
                             time_stepper_pt)
    {
      // RefineableTubeTFMesh's ctor builds the octree forest manually
      // without propagating the factory to each element. Do it here.
      unsigned n_element = this->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        if (auto* ref_el_pt =
              dynamic_cast<RefineableElement*>(this->element_pt(e)))
        {
          ref_el_pt->set_element_factory_pt(this->element_factory_pt);
        }
      }
    }
  };

} // namespace oomph
#endif
