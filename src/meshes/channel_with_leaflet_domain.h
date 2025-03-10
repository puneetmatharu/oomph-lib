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
#ifndef OOMPH_CHANNEL_WITH_LEAFLET_DOMAIN_HEADER
#define OOMPH_CHANNEL_WITH_LEAFLET_DOMAIN_HEADER

#include "generic/element_factory.h"

// Generic includes
#include "generic/geom_objects.h"
#include "generic/domain.h"
#include "generic/macro_element.h"

namespace oomph
{
  //===========================================================
  /// Rectangular domain with a leaflet blocking the lower
  /// half
  //===========================================================
  class ChannelWithLeafletDomain : public Domain
  {
  public:
    /// Constructor: Pass pointer to GeomObject that represents the
    /// leaflet,
    /// the length of the domain to left and right of the leaflet, the
    /// height of the leaflet and the overall height of the channel,
    /// the number of element columns to the left and right of the leaflet,
    /// the number of rows of elements from the bottom of the channel to
    /// the end of the leaflet, the number of rows of elements above the
    /// end of the leaflet.
    ChannelWithLeafletDomain(GeomObject* leaflet_pt,
                             const double& lleft,
                             const double& lright,
                             const double& hleaflet,
                             const double& htot,
                             const unsigned& nleft,
                             const unsigned& nright,
                             const unsigned& ny1,
                             const unsigned& ny2)
    {
      // Copy assignments into private storage
      Lleft = lleft;
      Lright = lright;
      Hleaflet = hleaflet;
      Htot = htot;
      Nleft = nleft;
      Nright = nright;
      Ny1 = ny1;
      Ny2 = ny2;
      Leaflet_pt = leaflet_pt;

      // Store origin of leaflet for fast reference
      Vector<double> zeta(1);
      zeta[0] = 0.0;
      Vector<double> r(2);
      Leaflet_pt->position(zeta, r);
      X_0 = r[0];

      // Number of macro elements
      unsigned nmacro = (Ny1 + Ny2) * (Nleft + Nright);
      Macro_element_pt.resize(nmacro);

      // Build the 2D macro elements
      for (unsigned i = 0; i < nmacro; i++)
      {
        Macro_element_pt[i] = new QMacroElement<2>(this, i);
      }
    }

    /// Destructor: Empty; cleanup done in base class
    ~ChannelWithLeafletDomain() {}

    /// Total height of domain (width of channel)
    double htot()
    {
      return Htot;
    }

    /// Height of leaflet
    double hleaflet()
    {
      return Hleaflet;
    }

    /// Length of domain to the left of leaflet
    double lleft()
    {
      return Lleft;
    }

    /// Length of domain to the right of leaflet
    double lright()
    {
      return Lright;
    }

    /// Pointer to the wall
    GeomObject*& leaflet_pt()
    {
      return Leaflet_pt;
    };

    /// Parametrisation of macro element boundaries
    void macro_element_boundary(const unsigned& t,
                                const unsigned& imacro,
                                const unsigned& idirect,
                                const Vector<double>& zeta,
                                Vector<double>& r);

  protected:
    /// Helper function
    void macro_bound_I_N(const unsigned& t,
                         const Vector<double>& zeta,
                         Vector<double>& r,
                         const unsigned& i,
                         const unsigned& j);

    /// Helper function
    void macro_bound_I_S(const unsigned& t,
                         const Vector<double>& zeta,
                         Vector<double>& r,
                         const unsigned& i,
                         const unsigned& j);

    /// Helper function
    void macro_bound_I_W(const unsigned& t,
                         const Vector<double>& zeta,
                         Vector<double>& r,
                         const unsigned& i,
                         const unsigned& j);

    /// Helper function
    void macro_bound_I_E(const unsigned& t,
                         const Vector<double>& zeta,
                         Vector<double>& r,
                         const unsigned& i,
                         const unsigned& j);

    /// Helper function
    void macro_bound_II_N(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r,
                          const unsigned& i,
                          const unsigned& j);

    /// Helper function
    void macro_bound_II_S(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r,
                          const unsigned& i,
                          const unsigned& j);

    /// Helper function
    void macro_bound_II_W(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r,
                          const unsigned& i,
                          const unsigned& j);

    /// Helper function
    void macro_bound_II_E(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r,
                          const unsigned& i,
                          const unsigned& j);

    /// Helper function
    void macro_bound_III_N(const unsigned& t,
                           const Vector<double>& zeta,
                           Vector<double>& r,
                           const unsigned& i,
                           const unsigned& j);

    /// Helper function
    void macro_bound_III_S(const unsigned& t,
                           const Vector<double>& zeta,
                           Vector<double>& r,
                           const unsigned& i,
                           const unsigned& j);

    /// Helper function
    void macro_bound_III_W(const unsigned& t,
                           const Vector<double>& zeta,
                           Vector<double>& r,
                           const unsigned& i,
                           const unsigned& j);

    /// Helper function
    void macro_bound_III_E(const unsigned& t,
                           const Vector<double>& zeta,
                           Vector<double>& r,
                           const unsigned& i,
                           const unsigned& j);

    /// Helper function
    void macro_bound_IV_N(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r,
                          const unsigned& i,
                          const unsigned& j);

    /// Helper function
    void macro_bound_IV_S(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r,
                          const unsigned& i,
                          const unsigned& j);

    /// Helper function
    void macro_bound_IV_W(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r,
                          const unsigned& i,
                          const unsigned& j);

    /// Helper function
    void macro_bound_IV_E(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r,
                          const unsigned& i,
                          const unsigned& j);
    /// Helper function
    void slanted_bound_up(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r);

    /// Length of the domain to the right of the leaflet
    double Lright;

    /// Length of the domain to the left of the leaflet
    double Lleft;

    /// Lagrangian coordinate at end of leaflet
    double Hleaflet;

    /// Total width of the channel
    double Htot;

    /// Number of macro element columnns to the right of the leaflet
    unsigned Nright;

    /// Number of macro element columns to the left of the leaflet
    unsigned Nleft;

    /// Number of macro element rows up to the end of the leaflet
    unsigned Ny1;

    /// Number of macro element rows above the leaflet
    unsigned Ny2;

    /// Center of the domain : origin of the leaflet, extracted
    /// from GeomObject and stored for fast access.
    double X_0;

    /// Pointer to leaflet
    GeomObject* Leaflet_pt;
  };

} // namespace oomph

#endif
