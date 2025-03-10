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
#ifndef OOMPH_QUARTER_TUBE_DOMAIN_HEADER
#define OOMPH_QUARTER_TUBE_DOMAIN_HEADER

#include "generic/element_factory.h"

// Generic oomph-lib includes
#include "generic/quadtree.h"
#include "generic/domain.h"
#include "generic/geom_objects.h"

namespace oomph
{
  //=================================================================
  /// Quarter tube as domain. Domain is bounded by
  /// curved boundary which is represented by a GeomObject. Domain is
  /// parametrised by three macro elements in each of the nlayer slices
  //=================================================================
  class QuarterTubeDomain : public Domain
  {
  public:
    /// Constructor: Pass boundary object and start and end coordinates
    /// and fraction along boundary object where outer ring is divided.
    /// We form nlayer axial slices.
    QuarterTubeDomain(GeomObject* boundary_geom_object_pt,
                      const Vector<double>& xi_lo,
                      const double& fract_mid,
                      const Vector<double>& xi_hi,
                      const unsigned& nlayer)
      : Xi_lo(xi_lo),
        Fract_mid(fract_mid),
        Xi_hi(xi_hi),
        Nlayer(nlayer),
        Wall_pt(boundary_geom_object_pt),
        BL_squash_fct_pt(&default_BL_squash_fct),
        Axial_spacing_fct_pt(&default_axial_spacing_fct)
    {
      // There are three macro elements
      unsigned nmacro = 3 * nlayer;

      // Resize
      Macro_element_pt.resize(nmacro);

      // Create macro elements
      for (unsigned i = 0; i < nmacro; i++)
      {
        Macro_element_pt[i] = new QMacroElement<3>(this, i);
      }
    }

    /// Broken copy constructor
    QuarterTubeDomain(const QuarterTubeDomain&) = delete;

    /// Broken assignment operator
    void operator=(const QuarterTubeDomain&) = delete;

    /// Destructor: empty; cleanup done in base class
    ~QuarterTubeDomain() {}

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
    /// coordinate to the return value.
    double s_squashed(const double& s)
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


    /// Vector representation of the  i_macro-th macro element
    /// boundary i_direct (L/R/D/U/B/F) at time level t
    /// (t=0: present; t>0: previous):
    /// f(s).
    void macro_element_boundary(const unsigned& t,
                                const unsigned& i_macro,
                                const unsigned& i_direct,
                                const Vector<double>& s,
                                Vector<double>& f);

  private:
    /// Lower limit for the coordinates along the wall
    Vector<double> Xi_lo;

    /// Fraction along wall where outer ring is to be divided
    double Fract_mid;

    /// Upper limit for the coordinates along the wall
    Vector<double> Xi_hi;

    /// Number of layers
    unsigned Nlayer;

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


    /// Function pointer for function that implements
    /// axial spacing of macro elements
    AxialSpacingFctPt Axial_spacing_fct_pt;


    /// Default for function that  implements
    /// axial spacing of macro elements
    static double default_axial_spacing_fct(const double& xi)
    {
      return xi;
    }


    /// Boundary of central box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_L(const unsigned& t,
                   const Vector<double>& zeta,
                   const unsigned& i_layer,
                   Vector<double>& f);

    /// Boundary of central box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_R(const unsigned& t,
                   const Vector<double>& zeta,
                   const unsigned& i_layer,
                   Vector<double>& f);

    /// Boundary of central box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_D(const unsigned& t,
                   const Vector<double>& zeta,
                   const unsigned& i_layer,
                   Vector<double>& f);

    /// Boundary of central box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_U(const unsigned& t,
                   const Vector<double>& zeta,
                   const unsigned& i_layer,
                   Vector<double>& f);

    /// Boundary of central box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_B(const unsigned& t,
                   const Vector<double>& zeta,
                   const unsigned& i_layer,
                   Vector<double>& f);

    /// Boundary of central box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_F(const unsigned& t,
                   const Vector<double>& zeta,
                   const unsigned& i_layer,
                   Vector<double>& f);


    /// Boundary of bottom right box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_bot_right_L(const unsigned& t,
                       const Vector<double>& zeta,
                       const unsigned& i_layer,
                       Vector<double>& f);

    /// Boundary of bottom right box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_bot_right_R(const unsigned& t,
                       const Vector<double>& zeta,
                       const unsigned& i_layer,
                       Vector<double>& f);

    /// Boundary of bottom right box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_bot_right_D(const unsigned& t,
                       const Vector<double>& zeta,
                       const unsigned& i_layer,
                       Vector<double>& f);

    /// Boundary of bottom right box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_bot_right_U(const unsigned& t,
                       const Vector<double>& zeta,
                       const unsigned& i_layer,
                       Vector<double>& f);

    /// Boundary of bottom right box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_bot_right_B(const unsigned& t,
                       const Vector<double>& zeta,
                       const unsigned& i_layer,
                       Vector<double>& f);

    /// Boundary of bottom right box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_bot_right_F(const unsigned& t,
                       const Vector<double>& zeta,
                       const unsigned& i_layer,
                       Vector<double>& f);


    /// Boundary of top left box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_top_left_L(const unsigned& t,
                      const Vector<double>& zeta,
                      const unsigned& i_layer,
                      Vector<double>& f);

    /// Boundary of top left box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_top_left_R(const unsigned& t,
                      const Vector<double>& zeta,
                      const unsigned& i_layer,
                      Vector<double>& f);

    /// Boundary of top left box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_top_left_D(const unsigned& t,
                      const Vector<double>& zeta,
                      const unsigned& i_layer,
                      Vector<double>& f);

    /// Boundary of top left box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_top_left_U(const unsigned& t,
                      const Vector<double>& zeta,
                      const unsigned& i_layer,
                      Vector<double>& f);

    /// Boundary of top left box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_top_left_B(const unsigned& t,
                      const Vector<double>& zeta,
                      const unsigned& i_layer,
                      Vector<double>& f);

    /// Boundary of top left box macro element in layer i_layer
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_top_left_F(const unsigned& t,
                      const Vector<double>& zeta,
                      const unsigned& i_layer,
                      Vector<double>& f);
  };

} // namespace oomph

#endif
