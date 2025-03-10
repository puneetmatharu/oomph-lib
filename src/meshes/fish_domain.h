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
#ifndef OOMPH_FISH_DOMAIN_HEADER
#define OOMPH_FISH_DOMAIN_HEADER

#include "generic/element_factory.h"

// Generic oomph-lib includes
#include "generic/quadtree.h"
#include "generic/domain.h"
#include "generic/geom_objects.h"

namespace oomph
{
  //===========start_of_fish_domain=======================================
  /// Fish shaped domain, represented by four
  /// MacroElements. Shape is parametrised by GeomObject
  /// that represents the fish's back.
  //=======================================================================
  class FishDomain : public Domain
  {
  public:
    /// Constructor: Pass pointer to GeomObject that represents the
    /// (upper) curved boundary of the fish's body, and the start and end values
    /// of the Lagrangian coordinates along the GeomObject.
    FishDomain(GeomObject* back_pt,
               const double& xi_nose,
               const double& xi_tail)
      : Xi_nose(xi_nose), Xi_tail(xi_tail), Back_pt(back_pt)
    {
      // Set values for private data members that are describe
      // geometric features of the fish: x-coordinate of the fin,
      // (half-)height of the fin, and x-position of the mouth.
      X_fin = 1.7;
      Y_fin = 0.9;
      X_mouth = 0.0;

      // There are four macro elements
      unsigned nmacro = 4;
      Macro_element_pt.resize(nmacro);

      // Build them
      for (unsigned i = 0; i < nmacro; i++)
      {
        Macro_element_pt[i] = new QMacroElement<2>(this, i);
      }
    } // end of constructor

    /// Broken copy constructor
    FishDomain(const FishDomain&) = delete;

    /// Broken assignment operator
    void operator=(const FishDomain&) = delete;

    /// Destructor for FishDomain: Empty; cleanup done in base class
    virtual ~FishDomain() {}

    /// x-position of fin tip
    double& x_fin()
    {
      return X_fin;
    }

    /// y-position of fin tip
    double& y_fin()
    {
      return Y_fin;
    }

    /// x-position of mouth
    double& x_mouth()
    {
      return X_mouth;
    }

    /// Start coordinate on wall (near nose)
    double& xi_nose()
    {
      return Xi_nose;
    }

    /// End coordinate on wall (near tail)
    double& xi_tail()
    {
      return Xi_tail;
    }

    /// Vector representation of the  i_macro-th macro element
    /// boundary i_direct (N/S/W/E) at the discrete time level t
    /// (t=0: present; t>0: previous): \f$ {\bf r}({\bf zeta}) \f$
    /// Note that the local coordinate \b zeta is a 1D
    /// Vector rather than a scalar -- this is unavoidable because
    /// this function implements the pure virtual function in the
    /// Domain base class.
    void macro_element_boundary(const unsigned& t,
                                const unsigned& i_macro,
                                const unsigned& i_direct,
                                const Vector<double>& zeta,
                                Vector<double>& r);

  private:
    /// "Nose" limit for the (1D) coordinates along the wall
    double Xi_nose;

    /// "Tail" limit for the (1D) coordinates along the wall
    double Xi_tail;

    /// X coordinate of fin tip
    double X_fin;

    /// Y coordinate of fin tip
    double Y_fin;

    /// X coordinate of corner of mouth
    double X_mouth;

    /// Pointer to the fish's back
    GeomObject* Back_pt;

    /// Boundary of upper body macro element zeta \f$ \in [-1,1] \f$
    void r_upper_body_N(const unsigned& t,
                        const Vector<double>& zeta,
                        Vector<double>& f);

    /// Boundary of upper body macro element zeta \f$ \in [-1,1] \f$
    void r_upper_body_W(const unsigned& t,
                        const Vector<double>& zeta,
                        Vector<double>& f);

    /// Boundary of upper body macro element zeta \f$ \in [-1,1] \f$
    void r_upper_body_S(const unsigned& t,
                        const Vector<double>& zeta,
                        Vector<double>& f);

    /// Boundary of upper body macro element zeta \f$ \in [-1,1] \f$
    void r_upper_body_E(const unsigned& t,
                        const Vector<double>& zeta,
                        Vector<double>& f);

    /// Boundary of upper fin macro element zeta \f$ \in [-1,1] \f$
    void r_upper_fin_N(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f);

    /// Boundary of upper fin macro element zeta \f$ \in [-1,1] \f$
    void r_upper_fin_W(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f);

    /// Boundary of upper fin macro element zeta \f$ \in [-1,1] \f$
    void r_upper_fin_S(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f);

    /// Boundary of upper fin macro element zeta \f$ \in [-1,1] \f$
    void r_upper_fin_E(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f);


    /// Boundary of lower body macro element zeta \f$ \in [-1,1] \f$
    void r_lower_body_N(const unsigned& t,
                        const Vector<double>& zeta,
                        Vector<double>& f)
    {
      // North of lower body is element is south of upper one.
      // Direction of the coordinate stays the same.
      r_upper_body_S(t, zeta, f);
      // Reflect vertical position
      f[1] = -f[1];
    }


    /// Boundary of lower body macro element zeta \f$ \in [-1,1] \f$
    void r_lower_body_W(const unsigned& t,
                        const Vector<double>& zeta,
                        Vector<double>& f)
    {
      // West of lower body is element is west of upper one.
      // Direction of the coordinate is inverted
      Vector<double> zeta_new(1);
      zeta_new[0] = -zeta[0];
      r_upper_body_W(t, zeta_new, f);
      // Vertical coordinate is reflected
      f[1] = -f[1];
    }

    /// Southern boundary of lower body macro element zeta \f$\in [-1,1] \f$
    void r_lower_body_S(const unsigned& t,
                        const Vector<double>& zeta,
                        Vector<double>& f)
    {
      // South of lower body is element is north of upper one.
      // Direction of the coordinate stays the same.
      r_upper_body_N(t, zeta, f);
      // Reflect vertical position
      f[1] = -f[1];
    }

    /// Boundary of lower body macro element zeta \f$ \in [-1,1] \f$
    void r_lower_body_E(const unsigned& t,
                        const Vector<double>& zeta,
                        Vector<double>& f)
    {
      // East of lower body is element is east of upper one.
      // Direction of the coordinate is inverted.
      Vector<double> zeta_new(1);
      zeta_new[0] = -zeta[0];
      r_upper_body_E(t, zeta_new, f);
      // Vertical coordinate is reflected
      f[1] = -f[1];
    }


    /// Boundary of lower fin macro element zeta \f$ \in [-1,1] \f$
    void r_lower_fin_N(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f)
    {
      // North of lower fin is element is south of upper one.
      // Direction of the coordinate stays the same.
      r_upper_fin_S(t, zeta, f);
      // Reflect vertical position
      f[1] = -f[1];
    }


    /// Boundary of lower fin macro element zeta \f$ \in [-1,1] \f$
    void r_lower_fin_W(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f)
    {
      // West of lower fin is element is west of upper one.
      // Direction of the coordinate is inverted
      Vector<double> zeta_new(1);
      zeta_new[0] = -zeta[0];
      r_upper_fin_W(t, zeta_new, f);
      // Vertical coordinate is reflected
      f[1] = -f[1];
    }

    /// Boundary of lower fin macro element zeta \f$ \in [-1,1] \f$
    void r_lower_fin_S(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f)
    {
      // South of lower fin is element is north of upper one.
      // Direction of the coordinate stays the same.
      r_upper_fin_N(t, zeta, f);
      // Reflect vertical position
      f[1] = -f[1];
    }

    /// Boundary of lower fin macro element zeta \f$ \in [-1,1] \f$
    void r_lower_fin_E(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f)
    {
      // East of lower fin is element is east of upper one.
      // Direction of the coordinate is inverted.
      Vector<double> zeta_new(1);
      zeta_new[0] = -zeta[0];
      r_upper_fin_E(t, zeta_new, f);
      // Vertical coordinate is reflected
      f[1] = -f[1];
    }
  };

} // namespace oomph

#endif
