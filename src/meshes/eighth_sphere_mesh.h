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
#ifndef OOMPH_EIGHTH_SPHERE_MESH_HEADER
#define OOMPH_EIGHTH_SPHERE_MESH_HEADER

#include "generic/element_factory.h"

// Headers
#include "generic/refineable_brick_mesh.h"
#include "generic/macro_element.h"
#include "generic/domain.h"
#include "generic/algebraic_elements.h"
#include "generic/brick_mesh.h"

// Include the headers file for domain
#include "eighth_sphere_domain.h"

namespace oomph
{
  //======================================================================
  /// Eight of a sphere brick mesh, based on the EightSphereDomain
  /// Non-refineable version with four brick elements.
  /// The eighth-sphere is located in the positive octant,
  /// centred at the origin. The mesh boundaries are numbered
  /// as follows:
  /// - Boundary 0: Plane x=0
  /// - Boundary 1: Plane y=0
  /// - Boundary 2: Plane z=0
  /// - Boundary 3: The surface of the sphere.
  //======================================================================
  class EighthSphereTFMesh : public virtual BrickMeshBase
  {
  public:
    /// Constructor: Pass radius and timestepper; defaults to
    /// static default timestepper
    EighthSphereTFMesh(
      const ElementFactoryPt& element_factory_pt,
      const double& radius,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// Destructor
    ~EighthSphereTFMesh()
    {
      delete Domain_pt;
      Domain_pt = 0;
    }

  protected:
    /// Pointer to the domain
    Domain* Domain_pt;

    /// Radius of the sphere
    double Radius;
  };

  //======================================================================
  /// Refineable version of the eight of a sphere brick mesh.
  /// The eighth-sphere is located in the positive octant,
  /// centred at the origin. The mesh boundaries are numbered
  /// as follows:
  /// - Boundary 0: Plane x=0
  /// - Boundary 1: Plane y=0
  /// - Boundary 2: Plane z=0
  /// - Boundary 3: The surface of the sphere.
  //======================================================================
  class RefineableEighthSphereTFMesh : public EighthSphereTFMesh,
                                       public virtual RefineableBrickTFMesh
  {
  public:
    /// Constructor: Pass radius and timestepper; defaults to
    /// static default timestepper
    RefineableEighthSphereTFMesh(
      const ElementFactoryPt& element_factory_pt,
      const double& radius,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : EighthSphereTFMesh(element_factory_pt, radius, time_stepper_pt)
    {
      // Loop over all elements and set macro element pointer
      unsigned nel = this->nelement();
      for (unsigned ielem = 0; ielem < nel; ielem++)
      {
        this->finite_element_pt(ielem)->set_macro_elem_pt(
          this->Domain_pt->macro_element_pt(ielem));
      }

      // Associate the elements with octrees and plant in forest
      this->setup_octree_forest();

#ifdef PARANOID
      // Run self test on octree forest
      dynamic_cast<OcTreeForest*>(this->Forest_pt)->self_test();
#endif
    }
  };


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //=========================================================================
  /// Backwards-compatible templated EighthSphereMesh wrappers.
  //=========================================================================
  template<class ELEMENT>
  class EighthSphereMesh : public EighthSphereTFMesh
  {
  public:
    EighthSphereMesh(const double& radius,
                     TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : EighthSphereTFMesh(
          make_element_factory_pt<ELEMENT>(), radius, time_stepper_pt)
    {
    }
  };


  template<class ELEMENT>
  class RefineableEighthSphereMesh : public RefineableEighthSphereTFMesh
  {
  public:
    RefineableEighthSphereMesh(
      const double& radius,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : EighthSphereTFMesh(
          make_element_factory_pt<ELEMENT>(), radius, time_stepper_pt),
        RefineableEighthSphereTFMesh(
          make_element_factory_pt<ELEMENT>(), radius, time_stepper_pt)
    {
    }
  };

} // namespace oomph

#endif
