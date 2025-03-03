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

#ifndef OOMPH_REFINEABLE_YOUNGLAPLACE_ELEMENTS_HEADER
#define OOMPH_REFINEABLE_YOUNGLAPLACE_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


// oomph-lib headers
#include "../generic/refineable_quad_element.h"
#include "../generic/error_estimator.h"
#include "young_laplace_elements.h"


namespace oomph
{
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////


  //======================================================================
  /// Refineable version of YoungLaplace equations
  ///
  ///
  //======================================================================
  class RefineableYoungLaplaceEquations
    : public virtual YoungLaplaceEquations,
      public virtual RefineableElement,
      public virtual ElementWithZ2ErrorEstimator
  {
  public:
    /// Constructor: Pass refinement level to refineable element
    /// (default 0 = root)
    RefineableYoungLaplaceEquations(const int& refine_level = 0)
      : YoungLaplaceEquations(),
        RefineableElement(),
        ElementWithZ2ErrorEstimator()
    {
    }


    /// Broken copy constructor
    RefineableYoungLaplaceEquations(
      const RefineableYoungLaplaceEquations& dummy) = delete;

    /// Broken assignment operator
    void operator=(const RefineableYoungLaplaceEquations&) = delete;

    /// Compute element residual vector taking hanging nodes into account
    void fill_in_contribution_to_residuals(Vector<double>& residuals);

    /// Number of 'flux' terms for Z2 error estimation
    unsigned num_Z2_flux_terms()
    {
      return 2;
    }

    /// Get 'flux' for Z2 error recovery:  Standard flux
    /// from YoungLaplace equations
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
      this->get_flux(s, flux);
    }

    /// Get the function value u in Vector.
    /// Note: Given the generality of the interface (this function
    /// is usually called from black-box documentation or interpolation
    /// routines), the values Vector sets its own size in here.
    void get_interpolated_values(const Vector<double>& s,
                                 Vector<double>& values)
    {
      // Set size of Vector: u
      values.resize(1);

      // Find number of nodes
      unsigned n_node = nnode();

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Initialise value of u
      values[0] = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        values[0] += this->u(l) * psi[l];
      }
    }


    /// Get the function value u in Vector.
    /// Note: Given the generality of the interface (this function
    /// is usually called from black-box documentation or interpolation
    /// routines), the values Vector sets its own size in here.
    void get_interpolated_values(const unsigned& t,
                                 const Vector<double>& s,
                                 Vector<double>& values)
    {
      if (t != 0)
      {
        throw OomphLibError("These equations are steady => No time dependence",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        get_interpolated_values(s, values);
      }
    }


    /// Further build: Copy function pointers to spine and spine base
    /// functions from father. Kappa is passed across with set_kappa(...)
    /// to ensure that it's added to the element's external Data.
    void further_build()
    {
      // Copy kappa with set_kappa() to ensure that it's added to the
      // element's external Data.
      this->set_kappa(dynamic_cast<RefineableYoungLaplaceEquations*>(
                        this->father_element_pt())
                        ->kappa_pt());

      // Copy spine functions
      this->Spine_fct_pt = dynamic_cast<RefineableYoungLaplaceEquations*>(
                             this->father_element_pt())
                             ->spine_fct_pt();

      this->Spine_base_fct_pt = dynamic_cast<RefineableYoungLaplaceEquations*>(
                                  this->father_element_pt())
                                  ->spine_base_fct_pt();
    }
  };


  //======================================================================
  /// Refineable version of 2D QYoungLaplaceElement elements
  ///
  ///
  //======================================================================
  template<unsigned NNODE_1D>
  class RefineableQYoungLaplaceElement
    : public QYoungLaplaceElement<NNODE_1D>,
      public virtual RefineableYoungLaplaceEquations,
      public virtual RefineableQElement<2>
  {
  public:
    /// Constructor: Pass refinement level to refineable quad element
    /// (default 0 = root)
    RefineableQYoungLaplaceElement()
      : RefineableElement(),
        RefineableYoungLaplaceEquations(),
        RefineableQElement<2>(),
        QYoungLaplaceElement<NNODE_1D>()
    {
    }

    /// Broken copy constructor
    RefineableQYoungLaplaceElement(
      const RefineableQYoungLaplaceElement<NNODE_1D>& dummy) = delete;

    /// Broken assignment operator
    void operator=(const RefineableQYoungLaplaceElement<NNODE_1D>&) = delete;

    /// Number of continuously interpolated values: 1
    unsigned ncont_interpolated_values() const
    {
      return 1;
    }

    /// Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return QYoungLaplaceElement<NNODE_1D>::nvertex_node();
    }

    /// Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return QYoungLaplaceElement<NNODE_1D>::vertex_node_pt(j);
    }

    /// Rebuild from sons: empty
    void rebuild_from_sons(Mesh*& mesh_pt) {}

    /// Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return (NNODE_1D - 1);
    }

    ///  Perform additional hanging node procedures for variables
    /// that are not interpolated by all nodes. Empty.
    void further_setup_hanging_nodes() {}
  };

  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Face geometry for the RefineableQuadYoungLaplaceElement elements: The
  /// spatial dimension of the face elements is one lower than that of the bulk
  /// element but they have the same number of points along their 1D edges.
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<RefineableQYoungLaplaceElement<NNODE_1D>>
    : public virtual QElement<1, NNODE_1D>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : QElement<1, NNODE_1D>() {}
  };


} // namespace oomph

#endif
