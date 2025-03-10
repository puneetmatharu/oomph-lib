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
#ifndef OOMPH_QUARTER_CIRCLE_SECTOR_DOMAIN_HEADER
#define OOMPH_QUARTER_CIRCLE_SECTOR_DOMAIN_HEADER

#include "generic/element_factory.h"

// Generic oomph-lib includes
#include "generic/quadtree.h"
#include "generic/domain.h"
#include "generic/geom_objects.h"

namespace oomph
{
  //=================================================================
  /// Circular sector as domain. Domain is bounded by
  /// curved boundary which is represented by a GeomObject. Domain is
  /// parametrised by three macro elements.
  //=================================================================
  class QuarterCircleSectorDomain : public Domain
  {
  public:
    /// Constructor: Pass boundary object and start and end coordinates
    /// and fraction along boundary object where outer ring is divided.
    QuarterCircleSectorDomain(GeomObject* boundary_geom_object_pt,
                              const double& xi_lo,
                              const double& fract_mid,
                              const double& xi_hi)
      : Xi_lo(xi_lo),
        Fract_mid(fract_mid),
        Xi_hi(xi_hi),
        Wall_pt(boundary_geom_object_pt),
        BL_squash_fct_pt(&default_BL_squash_fct)
    {
      // There are three macro elements
      unsigned nmacro = 3;

      // Resize
      Macro_element_pt.resize(nmacro);

      // Create macro elements
      for (unsigned i = 0; i < nmacro; i++)
      {
        Macro_element_pt[i] = new QMacroElement<2>(this, i);
      }
    }

    /// Broken copy constructor
    QuarterCircleSectorDomain(const QuarterCircleSectorDomain&) = delete;

    /// Broken assignment operator
    void operator=(const QuarterCircleSectorDomain&) = delete;

    /// Destructor: empty; cleanup done in base class
    ~QuarterCircleSectorDomain() {}

    /// Typedef for function pointer for function that squashes
    /// the outer two macro elements towards
    /// the wall by mapping the input value of the "radial" macro element
    /// coordinate to the return value
    typedef double (*BLSquashFctPt)(const double& s);

    /// Function pointer for function that squashes
    /// the outer two macro elements towards
    /// the wall by mapping the input value of the "radial" macro element
    /// coordinate to the return value
    BLSquashFctPt& bl_squash_fct_pt()
    {
      return BL_squash_fct_pt;
    }

    /// Function that squashes the outer two macro elements towards
    /// the wall by mapping the input value of the "radial" macro element
    /// coordinate to the return value
    double s_squashed(const double& s)
    {
      return BL_squash_fct_pt(s);
    }

    /// Vector representation of the  i_macro-th macro element
    /// boundary i_direct (N/S/W/E) at time level t
    /// (t=0: present; t>0: previous):
    /// f(s). Note that the local coordinate \b s is a 1D
    /// Vector rather than a scalar -- this is unavoidable because
    /// this function implements the pure virtual function in the
    /// Domain base class.
    void macro_element_boundary(const unsigned& t,
                                const unsigned& i_macro,
                                const unsigned& i_direct,
                                const Vector<double>& s,
                                Vector<double>& f);

  private:
    /// Lower limit for the (1D) coordinates along the wall
    double Xi_lo;

    /// Fraction along wall where outer ring is to be divided
    double Fract_mid;

    /// Upper limit for the (1D) coordinates along the wall
    double Xi_hi;

    /// Pointer to geometric object that represents the curved wall
    GeomObject* Wall_pt;

    /// Function pointer for function that squashes
    /// the outer two macro elements towards
    /// the wall by mapping the input value of the "radial" macro element
    /// coordinate to the return value
    BLSquashFctPt BL_squash_fct_pt;

    /// Default for function that squashes
    /// the outer two macro elements towards
    /// the wall by mapping the input value of the "radial" macro element
    /// coordinate to the return value: Identity.
    static double default_BL_squash_fct(const double& s)
    {
      return s;
    }

    /// Boundary of top left macro element zeta \f$ \in [-1,1] \f$
    void r_top_left_N(const unsigned& t,
                      const Vector<double>& zeta,
                      Vector<double>& f);

    /// Boundary of top left macro element zeta \f$ \in [-1,1] \f$
    void r_top_left_W(const unsigned& t,
                      const Vector<double>& zeta,
                      Vector<double>& f);

    /// Boundary of top left macro element zeta \f$ \in [-1,1] \f$
    void r_top_left_S(const unsigned& t,
                      const Vector<double>& zeta,
                      Vector<double>& f);

    /// Boundary of top left macro element zeta \f$ \in [-1,1] \f$
    void r_top_left_E(const unsigned& t,
                      const Vector<double>& zeta,
                      Vector<double>& f);

    /// Boundary of bottom right macro element zeta \f$ \in [-1,1] \f$
    void r_bot_right_N(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f);

    /// Boundary of bottom right macro element zeta \f$ \in [-1,1] \f$
    void r_bot_right_W(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f);

    /// Boundary of bottom right macro element zeta \f$ \in [-1,1] \f$
    void r_bot_right_S(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f);

    /// Boundary of bottom right macro element zeta \f$ \in [-1,1] \f$
    void r_bot_right_E(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f);

    /// Boundary of central box macro element zeta \f$ \in [-1,1] \f$
    void r_centr_N(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of central box macro element zeta \f$ \in [-1,1] \f$
    void r_centr_E(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of central box macro element zeta \f$ \in [-1,1] \f$
    void r_centr_S(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of central box macro element zeta \f$ \in [-1,1] \f$
    void r_centr_W(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);
  };

} // namespace oomph

#endif
