// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_ONE_D_MESH_HEADER
#define OOMPH_ONE_D_MESH_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// oomph-lib headers
#include "../generic/line_mesh.h"
#include "../generic/refineable_line_mesh.h"

namespace oomph
{
  //====================================================================
  /// 1D mesh consisting of N one-dimensional elements from the
  /// QElement family.
  /// \f[ x \in [Xmin,Xmax] \f]
  /// The mesh has two boundaries:
  /// - Boundary 0 is at \f$x=Xmin\f$.
  /// - Boundary 1 is at \f$x=Xmax\f$.
  /// .
  /// There is one node on each of these boundaries.
  //====================================================================
  template<class ELEMENT>
  class OneDMesh : public virtual LineMeshBase
  {
  public:
    /// Constructor: Pass number of elements, n_element, length of
    /// domain, length, and pointer to timestepper (defaults to a Steady
    /// timestepper so we don't need to specify one in problems without
    /// time-dependence).
    OneDMesh(const unsigned& n_element,
             const double& length,
             TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : Xmin(0.0), Xmax(length), N(n_element)
    {
      check_1d();

      build_mesh(time_stepper_pt);
    }

    /// Constructor: Pass number of elements, n_element, minimum
    /// coordinate, xmin, maximum coordinate, xmax, and a pointer to a
    /// timestepper.
    OneDMesh(const unsigned& n_element,
             const double& xmin,
             const double& xmax,
             TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : Xmin(xmin), Xmax(xmax), N(n_element)
    {
      check_1d();

      build_mesh(time_stepper_pt);
    }

  protected:
    /// Mesh can only be built with 1D elements (but can be either T or Q so
    /// can't use the normal assert_geometric_element function.
    void check_1d() const
    {
#ifdef PARANOID
      FiniteElement* el_pt = new ELEMENT;
      if (el_pt->dim() != 1)
      {
        std::string err = "OneDMesh is only for 1D elements";
        throw OomphLibError(
          err, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      delete el_pt;
      el_pt = 0;
#endif
    }

    /// Minimum coordinate
    double Xmin;

    /// Maximum coordinate
    double Xmax;

    /// Length of the domain
    double Length;

    /// Number of elements
    unsigned N;

    /// Generic mesh constuction routine, called by all constructors
    void build_mesh(TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);
  };


  //====================================================================
  /// Refineable version of the OneDMesh
  //====================================================================
  template<class ELEMENT>
  class RefineableOneDMesh : public virtual OneDMesh<ELEMENT>,
                             public RefineableLineMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass number of elements, n_element, length of
    /// domain, length, and pointer to timestepper (defaults to Steady)
    RefineableOneDMesh(
      const unsigned& n_element,
      const double& length,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : OneDMesh<ELEMENT>(n_element, length, time_stepper_pt)
    {
      // Nodal positions etc. were created in constructor for OneDMesh<...>
      // so only need to set up binary tree forest
      this->setup_binary_tree_forest();
    }

    /// Constructor that allows the specification of minimum and
    /// maximum values of x coordinates. Also pass pointer to timestepper
    /// (defaults to Steady).
    RefineableOneDMesh(
      const unsigned& n_element,
      const double& xmin,
      const double& xmax,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : OneDMesh<ELEMENT>(n_element, xmin, xmax, time_stepper_pt)
    {
      // Nodal positions etc. were created in constructor for OneDMesh<...>
      // so only need to set up binary tree forest
      this->setup_binary_tree_forest();
    }
  };


} // namespace oomph

#endif
