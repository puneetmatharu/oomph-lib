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
#ifndef OOMPH_COLLAPSIBLE_CHANNEL_DOMAIN_HEADER
#define OOMPH_COLLAPSIBLE_CHANNEL_DOMAIN_HEADER

#include "generic/element_factory.h"

// Generic oomph-lib includes
#include "generic/quadtree.h"
#include "generic/domain.h"
#include "generic/geom_objects.h"

namespace oomph
{
  //==================================================================
  /// Collapsible channel domain
  //==================================================================
  class CollapsibleChannelDomain : public Domain
  {
  public:
    /// Constructor: Pass the number of (macro-)elements,
    /// the domain lengths in the x- and y-direction
    /// and the pointer to the geometric object that specifies
    /// the shape of the "collapsible" segment.
    CollapsibleChannelDomain(const unsigned& nup,
                             const unsigned& ncollapsible,
                             const unsigned& ndown,
                             const unsigned& ny,
                             const double& lup,
                             const double& lcollapsible,
                             const double& ldown,
                             const double& ly,
                             GeomObject* wall_pt)
      : BL_squash_fct_pt(&default_BL_squash_fct),
        Axial_spacing_fct_pt(&default_axial_spacing_fct),
        Rotate_domain(false)
    {
      Nup = nup;
      Ncollapsible = ncollapsible;
      Ndown = ndown;
      Ny = ny;
      Lup = lup;
      Lcollapsible = lcollapsible;
      Ldown = ldown;
      Ly = ly;
      Wall_pt = wall_pt;

      // Total number of macro elements
      unsigned nmacro = (Nup + Ncollapsible + Ndown) * Ny;

      // Build the macro elements
      Macro_element_pt.resize(nmacro);
      for (unsigned i = 0; i < nmacro; i++)
      {
        Macro_element_pt[i] = new QMacroElement<2>(this, i);
      }
    }

    /// Destructor: emtpy; cleanup done in base class
    ~CollapsibleChannelDomain() {}

    /// Number of vertical columns of macro elements the upstream section
    unsigned nup()
    {
      return Nup;
    }

    /// Number of vertical clumns of macro elements in the "collapsible" segment
    unsigned ncollapsible()
    {
      return Ncollapsible;
    }

    /// Number of vertical columns of macro elements in the downstream section
    unsigned ndown()
    {
      return Ndown;
    }

    /// Number of macro-elements across the channel
    unsigned ny()
    {
      return Ny;
    }

    /// Length of upstream section
    double l_up()
    {
      return Lup;
    }

    /// Length of collapsible segment
    double l_collapsible()
    {
      return Lcollapsible;
    }

    /// Length of downstream section
    double l_down()
    {
      return Ldown;
    }

    /// Width of channel
    double l_y()
    {
      return Ly;
    }

    /// Access to pointer to the geometric object that parametrises
    /// the collapsible wall
    GeomObject*& wall_pt()
    {
      return Wall_pt;
    }


    /// Access to pointer to the geometric object that parametrises
    /// the collapsible wall (const version)
    GeomObject* wall_pt() const
    {
      return Wall_pt;
    }

    /// Typedef for function pointer for function that squashes
    /// the macro elements near the wall to help resolution of any
    /// wall boundary layers.
    typedef double (*BLSquashFctPt)(const double& s);


    /// Default for function that squashes
    /// the macro elements near the walls. Identity.
    static double default_BL_squash_fct(const double& s)
    {
      return s;
    }

    /// Function pointer for function that squashes
    /// the macro elements near wall. Default mapping (identity)
    /// leaves the y-coordinate of the nodal points unchanged.
    BLSquashFctPt& bl_squash_fct_pt()
    {
      return BL_squash_fct_pt;
    }

    /// Function that squashes the macro elements near the wall.
    /// Input argument should vary between 0 and 1; function should return
    /// stretched/squashed coordinate in the same range. Default implementation
    /// is the identity; can be overloaded by specifying a different
    /// function pointer with bl_squash_fct_pt().
    double s_squash(const double& s)
    {
      return BL_squash_fct_pt(s);
    }

    /// Typedef for function pointer for function that implements
    /// axial spacing of macro elements
    typedef double (*AxialSpacingFctPt)(const double& xi);

    /// Function pointer for function that  implements
    /// axial spacing of macro elements
    AxialSpacingFctPt& axial_spacing_fct_pt()
    {
      return Axial_spacing_fct_pt;
    }

    /// Function that implements
    /// axial spacing of macro elements
    double axial_spacing_fct(const double& xi)
    {
      return Axial_spacing_fct_pt(xi);
    }


    /// Vector representation of the  imacro-th macro element
    /// boundary idirect (N/S/W/E) at time level t
    /// (t=0: present; t>0: previous): \f$ {\bf r}({\bf zeta}) \f$
    /// Note that the local coordinate \b zeta is a 1D
    /// Vector rather than a scalar -- this is unavoidable because
    /// this function implements the pure virtual function in the
    /// Domain base class.
    void macro_element_boundary(const unsigned& t,
                                const unsigned& imacro,
                                const unsigned& idirect,
                                const Vector<double>& zeta,
                                Vector<double>& r);

    void enable_rotate_domain()
    {
      Rotate_domain = true;
    }

    /// Undo rotation of the domain (for axisymmetric problems)
    void disable_rotate_domain()
    {
      Rotate_domain = false;
    }

  private:
    /// Northern boundary of the macro element imacro in the
    /// upstream (part=0) or downstream (part=1) sections
    void r_N_straight(const Vector<double>& zeta,
                      Vector<double>& r,
                      const unsigned& imacro,
                      const unsigned& part);

    /// Western boundary of the macro element imacro in the
    /// upstream (part=0) or downstream (part=1) sections
    void r_W_straight(const Vector<double>& zeta,
                      Vector<double>& r,
                      const unsigned& imacro,
                      const unsigned& part);

    /// Southern boundary of the macro element imacro in the
    ///  upstream (part=0) or downstream (part=1) sections
    void r_S_straight(const Vector<double>& zeta,
                      Vector<double>& r,
                      const unsigned& imacro,
                      const unsigned& part);

    /// Eastern boundary of the macro element imacro in the
    /// upstream (part=0) or downstream (part=1) sections
    void r_E_straight(const Vector<double>& zeta,
                      Vector<double>& r,
                      const unsigned& imacro,
                      const unsigned& part);

    /// Northern boundary of the macro element imacro in the collapsible
    /// section
    void r_N_collapsible(const unsigned& t,
                         const Vector<double>& zeta,
                         Vector<double>& r,
                         const unsigned& imacro);

    /// Western boundary of the macro element imacro in the collapsible
    /// section
    void r_W_collapsible(const unsigned& t,
                         const Vector<double>& zeta,
                         Vector<double>& r,
                         const unsigned& imacro);

    /// Southern boundary of the macro element imacro in the collapsible
    /// section
    void r_S_collapsible(const unsigned& t,
                         const Vector<double>& zeta,
                         Vector<double>& r,
                         const unsigned& imacro);

    /// Eastern boundary of the macro element imacro in the collapsible
    /// section
    void r_E_collapsible(const unsigned& t,
                         const Vector<double>& zeta,
                         Vector<double>& r,
                         const unsigned& imacro);


    /// Function pointer for function that squashes
    /// the macro elements near the walls
    BLSquashFctPt BL_squash_fct_pt;

    /// Function pointer for function that implements
    /// axial spacing of macro elements
    AxialSpacingFctPt Axial_spacing_fct_pt;

    /// Default for function that  implements
    /// axial spacing of macro elements
    static double default_axial_spacing_fct(const double& xi)
    {
      return xi;
    }

    /// Number of vertical element columns in upstream section
    unsigned Nup;

    /// Number of vertical element columns in "collapsible" section
    unsigned Ncollapsible;

    /// Number of vertical element columns in downstream section
    unsigned Ndown;

    /// Number of macro elements across channel
    unsigned Ny;

    /// x-length in the upstream part of the channel
    double Lup;

    /// x-length in the "collapsible" part of the channel
    double Lcollapsible;

    /// x-length in the downstream part of the channel
    double Ldown;

    /// Width
    double Ly;

    /// Pointer to the geometric object that parametrises the collapsible wall
    GeomObject* Wall_pt;

    /// Rotate domain (for axisymmetric problems, say)
    bool Rotate_domain;
  };

} // namespace oomph

#endif
