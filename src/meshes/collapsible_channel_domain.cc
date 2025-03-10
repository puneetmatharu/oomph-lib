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
// Include guards

#include "collapsible_channel_domain.h"

namespace oomph
{

  //===================================================================
  /// Vector representation of the  imacro-th macro element
  /// boundary idirect (N/S/W/E) at time level t
  /// (t=0: present; t>0: previous): \f$ {\bf r}({\bf zeta}) \f$
  /// Note that the local coordinate \b zeta is a 1D
  /// Vector rather than a scalar -- this is unavoidable because
  /// this function implements the pure virtual function in the
  /// Domain base class.
  //=================================================================
  void CollapsibleChannelDomain::macro_element_boundary(
    const unsigned& t,
    const unsigned& imacro,
    const unsigned& idirect,
    const Vector<double>& zeta,
    Vector<double>& r)
  {
    using namespace QuadTreeNames;

#ifdef WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES
    // Warn about time argument being moved to the front
    OomphLibWarning(
      "Order of function arguments has changed between versions 0.8 and 0.85",
      "CollapsibleChannelDomain::macro_element_boundary(...)",
      OOMPH_EXCEPTION_LOCATION);
#endif

    // Total number of vertical columns of (macro-)elements
    unsigned n_x = Nup + Ncollapsible + Ndown;

    // Which direction?
    if (idirect == N)
    {
      // Upstream part of the channel
      if ((imacro % n_x) < Nup)
      {
        r_N_straight(zeta, r, imacro, 0);
      }
      // Downstream part of channel
      else if ((imacro % n_x) >= Nup + Ncollapsible)
      {
        r_N_straight(zeta, r, imacro, 1);
      }
      // Collapsible segment
      else if (((imacro % n_x) < Nup + Ncollapsible) && ((imacro % n_x) >= Nup))
      {
        r_N_collapsible(t, zeta, r, imacro);
      }
      else
      {
        std::ostringstream error_stream;
        error_stream << "Never get here! imacro, idirect: " << imacro << " "
                     << idirect << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
    else if (idirect == S)
    {
      // Upstream part
      if ((imacro % n_x) < Nup)
      {
        r_S_straight(zeta, r, imacro, 0);
      }
      // Downstream part
      else if ((imacro % n_x) >= Nup + Ncollapsible)
      {
        r_S_straight(zeta, r, imacro, 1);
      }
      // "Collapsible" bit
      else if (((imacro % n_x) < Nup + Ncollapsible) && ((imacro % n_x) >= Nup))
      {
        r_S_collapsible(t, zeta, r, imacro);
      }
      else
      {
        std::ostringstream error_stream;
        error_stream << "Never get here! imacro, idirect: " << imacro << " "
                     << idirect << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
    else if (idirect == E)
    {
      // Upstream bit
      if ((imacro % n_x) < Nup)
      {
        r_E_straight(zeta, r, imacro, 0);
      }
      // Downstream bit
      else if ((imacro % n_x) >= Nup + Ncollapsible)
      {
        r_E_straight(zeta, r, imacro, 1);
      }
      // "Collapsible" bit
      else if (((imacro % n_x) < Nup + Ncollapsible) && ((imacro % n_x) >= Nup))
      {
        r_E_collapsible(t, zeta, r, imacro);
      }
      else
      {
        std::ostringstream error_stream;
        error_stream << "Never get here! imacro, idirect: " << imacro << " "
                     << idirect << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }

    else if (idirect == W)
    {
      // Upstream bit
      if ((imacro % n_x) < Nup)
      {
        r_W_straight(zeta, r, imacro, 0);
      }
      // Downstream bit
      else if ((imacro % n_x) >= Nup + Ncollapsible)
      {
        r_W_straight(zeta, r, imacro, 1);
      }
      // "Collapsible" bit
      else if (((imacro % n_x) < Nup + Ncollapsible) && ((imacro % n_x) >= Nup))
      {
        r_W_collapsible(t, zeta, r, imacro);
      }
      else
      {
        std::ostringstream error_stream;
        error_stream << "Never get here! imacro, idirect: " << imacro << " "
                     << idirect << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }

    // Rotate?
    if (Rotate_domain)
    {
      double radius = r[1];
      double z = r[0];

      r[0] = radius;
      r[1] = -z;
    }
  }

  //===========================================================================
  /// Western edge of the  macro element in the upstream (part=0)
  /// or downstream (part=1) parts of the channel; \f$ \zeta \in [-1,1] \f$
  //===========================================================================
  void CollapsibleChannelDomain::r_W_straight(const Vector<double>& zeta,
                                              Vector<double>& r,
                                              const unsigned& imacro,
                                              const unsigned& part)
  {
    // Determines the "coordinates" of the macro-element
    unsigned x = unsigned(imacro % (Nup + Ncollapsible + Ndown));
    unsigned y = unsigned(double(imacro) / double(Nup + Ncollapsible + Ndown));

    // Where are we?
    switch (part)
    {
      case 0: // in the upstream part of the channel

        // Parametrize the boundary
        r[0] = axial_spacing_fct(double(x) * (Lup / double(Nup)));
        r[1] = (double(y) + (0.5 * (1.0 + zeta[0]))) * (Ly / double(Ny));

        // Map it via squash fct
        r[1] = Ly * s_squash(r[1] / Ly);

        break;

      case 1: // in the downstream part of the channel

        // Parametrizes the boundary
        r[0] = axial_spacing_fct(double(x - Nup - Ncollapsible) *
                                   (Ldown / double(Ndown)) +
                                 Lup + Lcollapsible);
        r[1] = (double(y) + (0.5 * (1.0 + zeta[0]))) * (Ly / double(Ny));

        // Map it via squash fct
        r[1] = Ly * s_squash(r[1] / Ly);

        break;

      default:

        std::ostringstream error_stream;
        error_stream << "Never get here! part=" << part << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }

  //===========================================================================
  /// Eastern  edge of the  macro element in the straight parts
  /// of the channel; \f$ \zeta \in [-1,1] \f$
  /// part=0 in the upstream part, part=1 in the downstream part
  //===========================================================================
  void CollapsibleChannelDomain::r_E_straight(const Vector<double>& zeta,
                                              Vector<double>& r,
                                              const unsigned& imacro,
                                              const unsigned& part)
  {
    // Determines the "coordinates" of the macro-element
    unsigned x = unsigned(imacro % (Nup + Ncollapsible + Ndown));
    unsigned y = unsigned(double(imacro) / double(Nup + Ncollapsible + Ndown));

    // Where are we?
    switch (part)
    {
      case 0: // in the upstream part of the channel

        // Parametrizes the boundary
        r[0] = axial_spacing_fct((double(x) + 1.0) * (Lup / double(Nup)));
        r[1] = (double(y) + (0.5 * (1.0 + zeta[0]))) * (Ly / double(Ny));

        // Map it via squash fct
        r[1] = Ly * s_squash(r[1] / Ly);

        break;

      case 1: // in the downstream part of the channel

        // Parametrizes the boundary
        r[0] = axial_spacing_fct((double(x - Nup - Ncollapsible) + 1.0) *
                                   (Ldown / double(Ndown)) +
                                 Lup + Lcollapsible);
        r[1] = (double(y) + (0.5 * (1.0 + zeta[0]))) * (Ly / double(Ny));

        // Map it via squash fct
        r[1] = Ly * s_squash(r[1] / Ly);

        break;

      default:

        std::ostringstream error_stream;
        error_stream << "Never get here! part=" << part << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }

  //==========================================================================
  /// Northern edge of the  macro element in the straight parts of
  /// the channel; \f$ \zeta \in [-1,1] \f$
  /// part=0 in the left part, part=1 in the right part
  //==========================================================================
  void CollapsibleChannelDomain::r_N_straight(const Vector<double>& zeta,
                                              Vector<double>& r,
                                              const unsigned& imacro,
                                              const unsigned& part)
  {
    // Determines the "coordinates" of the macro-element
    unsigned x = unsigned(imacro % (Nup + Ncollapsible + Ndown));
    unsigned y = unsigned(double(imacro) / double(Nup + Ncollapsible + Ndown));

    // Where are we?
    switch (part)
    {
      case 0: // in the upstream part of the channel

        // Parametrizes the boundary
        r[0] = axial_spacing_fct((double(x) + (0.5 * (1.0 + zeta[0]))) *
                                 (Lup / double(Nup)));
        r[1] = (double(y) + 1.0) * (Ly / double(Ny));

        // Map it via squash fct
        r[1] = Ly * s_squash(r[1] / Ly);

        break;

      case 1: // in the downstream part of the channel

        // Parametrizes the boundary
        r[0] = axial_spacing_fct(
          (double(x - Nup - Ncollapsible) + (0.5 * (1.0 + zeta[0]))) *
            (Ldown / double(Ndown)) +
          Lup + Lcollapsible);
        r[1] = (double(y) + 1.0) * (Ly / double(Ny));

        // Map it via squash fct
        r[1] = Ly * s_squash(r[1] / Ly);

        break;

      default:

        std::ostringstream error_stream;
        error_stream << "Never get here! part=" << part << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }

  //=========================================================================
  /// Southern  edge of the  macro element in the straight parts of
  /// the channel; \f$ \zeta \in [-1,1] \f$
  /// part=0 in the left part, part=1 in the right part
  //=========================================================================
  void CollapsibleChannelDomain::r_S_straight(const Vector<double>& zeta,
                                              Vector<double>& r,
                                              const unsigned& imacro,
                                              const unsigned& part)
  {
    // Determines the "coordinates" of the macro-element
    unsigned x = unsigned(imacro % (Nup + Ncollapsible + Ndown));
    unsigned y = unsigned(double(imacro) / double(Nup + Ncollapsible + Ndown));

    // Where are we?
    switch (part)
    {
      case 0: // in the upstream bit

        // Parametrizes the boundary
        r[0] = axial_spacing_fct((double(x) + (0.5 * (1 + zeta[0]))) *
                                 (Lup / double(Nup)));
        r[1] = double(y) * (Ly / double(Ny));

        // Map it via squash fct
        r[1] = Ly * s_squash(r[1] / Ly);

        break;

      case 1: // in the downstream bit

        // Parametrizes the boundary
        r[0] = axial_spacing_fct(
          (double(x - Nup - Ncollapsible) + (0.5 * (1 + zeta[0]))) *
            (Ldown / double(Ndown)) +
          Lup + Lcollapsible);
        r[1] = double(y) * (Ly / double(Ny));

        // Map it via squash fct
        r[1] = Ly * s_squash(r[1] / Ly);

        break;

      default:

        std::ostringstream error_stream;
        error_stream << "Never get here! part=" << part << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }

  //========================================================================
  /// Western edge of the  macro element in the collapsible part of the
  /// channel; \f$ \zeta \in [-1,1] \f$.
  //========================================================================
  void CollapsibleChannelDomain::r_W_collapsible(const unsigned& t,
                                                 const Vector<double>& zeta,
                                                 Vector<double>& r,
                                                 const unsigned& imacro)
  {
    // Determines the "coordinates" of the macro-element
    unsigned x = unsigned(imacro % (Nup + Ncollapsible + Ndown));
    unsigned y = unsigned(double(imacro) / double(Nup + Ncollapsible + Ndown));

    // Vector of Lagrangian coordinates
    Vector<double> xi(1);
    xi[0] = double(x - Nup) * (Lcollapsible / double(Ncollapsible));

    // Position vector on upper wall:
    Vector<double> r_wall(2);
    Wall_pt->position(t, xi, r_wall);

    // Point will be located on straight line from bottom to top wall
    double fract = (double(y) + (0.5 * (1.0 + zeta[0]))) / double(Ny);

    // Map it via squash fct
    fract = s_squash(fract);

    // x-cooordinate -- straight line from fixed position on the bottom
    // wall to moving position on the top wall
    r[0] = Lup + xi[0] + (r_wall[0] - (xi[0] + Lup)) * fract;

    // y-coordinate
    r[1] = r_wall[1] * fract;
  }

  //=========================================================================
  /// Eastern  edge of the  macro element in the collapsible part of the
  /// channel; \f$ \zeta \in [-1,1] \f$
  //=========================================================================
  void CollapsibleChannelDomain::r_E_collapsible(const unsigned& t,
                                                 const Vector<double>& zeta,
                                                 Vector<double>& r,
                                                 const unsigned& imacro)
  {
    // Determines the "coordinates" of the macro-element
    unsigned x = unsigned(imacro % (Nup + Ncollapsible + Ndown));
    unsigned y = unsigned(double(imacro) / double(Nup + Ncollapsible + Ndown));

    // Vector of Lagrangian coordinates
    Vector<double> xi(1);
    xi[0] = (double(x - Nup) + 1.0) * (Lcollapsible / double(Ncollapsible));

    // Position vector on upper wall:
    Vector<double> r_wall(2);
    Wall_pt->position(t, xi, r_wall);

    // Point will be located on straight line from bottom to top wall
    double fract = (double(y) + (0.5 * (1.0 + zeta[0]))) / double(Ny);

    // Map it via squash fct
    fract = s_squash(fract);

    // x-cooordinate -- straight line from fixed position on the bottom
    // wall to moving position on the top wall
    r[0] = Lup + xi[0] + (r_wall[0] - (xi[0] + Lup)) * fract;

    // y-coordinate
    r[1] = r_wall[1] * fract;
  }

  //==========================================================================
  /// Northern edge of the  macro element in the collapsible part of the
  /// channel; \f$ \zeta \in [-1,1] \f$
  //==========================================================================
  void CollapsibleChannelDomain::r_N_collapsible(const unsigned& t,
                                                 const Vector<double>& zeta,
                                                 Vector<double>& r,
                                                 const unsigned& imacro)
  {
    // Determines the "coordinates" of the macro-element
    unsigned x = unsigned(imacro % (Nup + Ncollapsible + Ndown));
    unsigned y = unsigned(double(imacro) / double(Nup + Ncollapsible + Ndown));

    // Vector of Lagrangian coordinates
    Vector<double> xi(1);
    xi[0] = (double(x - Nup) + (0.5 * (1.0 + zeta[0]))) *
            (Lcollapsible / double(Ncollapsible));

    // Position vector on upper wall:
    Vector<double> r_wall(2);
    Wall_pt->position(t, xi, r_wall);

    // Point will be located on straight line from bottom to top wall
    double fract = (double(y) + 1.0) / double(Ny);

    // Map it via squash fct
    fract = s_squash(fract);

    // x-cooordinate -- straight line from fixed position on the bottom
    // wall to moving position on the top wall
    r[0] = Lup + xi[0] + (r_wall[0] - (xi[0] + Lup)) * fract;

    // y-coordinate
    r[1] = r_wall[1] * fract;
  }

  //========================================================================
  /// Southern  edge of the  macro element in the collapsible part of the
  /// channel; \f$ \zeta \in [-1,1] \f$
  //========================================================================
  void CollapsibleChannelDomain::r_S_collapsible(const unsigned& t,
                                                 const Vector<double>& zeta,
                                                 Vector<double>& r,
                                                 const unsigned& imacro)
  {
    // Determines the "coordinates" of the macro-element
    unsigned x = unsigned(imacro % (Nup + Ncollapsible + Ndown));
    unsigned y = unsigned(double(imacro) / double(Nup + Ncollapsible + Ndown));

    // Vector of Lagrangian coordinates
    Vector<double> xi(1);
    xi[0] = (double(x - Nup) + (0.5 * (1.0 + zeta[0]))) *
            (Lcollapsible / double(Ncollapsible));

    // Position vector on upper wall:
    Vector<double> r_wall(2);
    Wall_pt->position(t, xi, r_wall);

    // Point will be located on straight line from bottom to top wall
    double fract = double(y) / double(Ny);

    // Map it via squash fct
    fract = s_squash(fract);

    // x-cooordinate -- straight line from fixed position on the bottom
    // wall to moving position on the top wall
    r[0] = Lup + xi[0] + (r_wall[0] - (xi[0] + Lup)) * fract;

    // y-coordinate
    r[1] = r_wall[1] * fract;
  }

} // namespace oomph
