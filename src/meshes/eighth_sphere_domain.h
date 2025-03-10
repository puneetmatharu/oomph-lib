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
#ifndef OOMPH_EIGHTH_SPHERE_DOMAIN_HEADER
#define OOMPH_EIGHTH_SPHERE_DOMAIN_HEADER

#include "generic/element_factory.h"

// Generic oomph-lib includes
#include "generic/quadtree.h"
#include "generic/domain.h"
#include "generic/geom_objects.h"

namespace oomph
{
  //=================================================================
  /// Eighth sphere as domain. Domain is
  /// parametrised by four macro elements
  //=================================================================
  class EighthSphereDomain : public Domain
  {
  public:
    /// Constructor: Pass the radius of the sphere.
    EighthSphereDomain(const double& radius) : Radius(radius)
    {
      // There are four macro elements
      unsigned nmacro = 4;

      // Resize
      Macro_element_pt.resize(nmacro);

      // Create the 3D Q macro elements
      for (unsigned i = 0; i < nmacro; i++)
      {
        Macro_element_pt[i] = new QMacroElement<3>(this, i);
      }
    }

    /// Broken copy constructor
    EighthSphereDomain(const EighthSphereDomain&) = delete;

    /// Broken assignment operator
    void operator=(const EighthSphereDomain&) = delete;

    /// Destructor: Empty; cleanup done in base class
    ~EighthSphereDomain() {}


    /// Vector representation of the  imacro-th macro element
    /// boundary idirect (L/R/D/U/B/F) at time level t
    /// (t=0: present; t>0: previous):
    /// f(s).
    void macro_element_boundary(const unsigned& t,
                                const unsigned& imacro,
                                const unsigned& idirect,
                                const Vector<double>& s,
                                Vector<double>& f)
    {
      using namespace OcTreeNames;

#ifdef WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES
      // Warn about time argument being moved to the front
      OomphLibWarning(
        "Order of function arguments has changed between versions 0.8 and 0.85",
        "EighthSphereDomain::macro_element_boundary(...)",
        OOMPH_EXCEPTION_LOCATION);
#endif

      // Which macro element?
      // --------------------
      switch (imacro)
      {
          // Macro element 0: Central box
        case 0:

          if (idirect == L)
          {
            r_centr_L(t, s, f);
          }
          else if (idirect == R)
          {
            r_centr_R(t, s, f);
          }
          else if (idirect == D)
          {
            r_centr_D(t, s, f);
          }
          else if (idirect == U)
          {
            r_centr_U(t, s, f);
          }
          else if (idirect == B)
          {
            r_centr_B(t, s, f);
          }
          else if (idirect == F)
          {
            r_centr_F(t, s, f);
          }
          else
          {
            std::ostringstream error_message;
            error_message << "idirect is " << OcTree::Direct_string[idirect]
                          << "not one of L, R, U, D, B, F" << std::endl;

            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }

          break;

          // Macro element 1:right
        case 1:

          // Which direction?
          if (idirect == L)
          {
            r_right_L(t, s, f);
          }
          else if (idirect == R)
          {
            r_right_R(t, s, f);
          }
          else if (idirect == D)
          {
            r_right_D(t, s, f);
          }
          else if (idirect == U)
          {
            r_right_U(t, s, f);
          }
          else if (idirect == B)
          {
            r_right_B(t, s, f);
          }
          else if (idirect == F)
          {
            r_right_F(t, s, f);
          }
          else
          {
            std::ostringstream error_message;
            error_message << "idirect is " << OcTree::Direct_string[idirect]
                          << "not one of L, R, U, D, B, F" << std::endl;

            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }

          break;

          // Macro element 2: Up
        case 2:

          // Which direction?
          if (idirect == L)
          {
            r_up_L(t, s, f);
          }
          else if (idirect == R)
          {
            r_up_R(t, s, f);
          }
          else if (idirect == D)
          {
            r_up_D(t, s, f);
          }
          else if (idirect == U)
          {
            r_up_U(t, s, f);
          }
          else if (idirect == B)
          {
            r_up_B(t, s, f);
          }
          else if (idirect == F)
          {
            r_up_F(t, s, f);
          }
          else
          {
            std::ostringstream error_message;
            error_message << "idirect is " << OcTree::Direct_string[idirect]
                          << "not one of L, R, U, D, B, F" << std::endl;

            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }

          break;

          // Macro element 3: Front
        case 3:
          // Which direction?
          if (idirect == L)
          {
            r_front_L(t, s, f);
          }
          else if (idirect == R)
          {
            r_front_R(t, s, f);
          }
          else if (idirect == D)
          {
            r_front_D(t, s, f);
          }
          else if (idirect == U)
          {
            r_front_U(t, s, f);
          }
          else if (idirect == B)
          {
            r_front_B(t, s, f);
          }
          else if (idirect == F)
          {
            r_front_F(t, s, f);
          }
          else
          {
            std::ostringstream error_message;
            error_message << "idirect is " << OcTree::Direct_string[idirect]
                          << "not one of L, R, U, D, B, F" << std::endl;

            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }

          break;

        default:

          std::ostringstream error_message;
          error_message << "imacro is " << OcTree::Direct_string[idirect]
                        << ", but should be in the range 0-3" << std::endl;

          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }

  private:
    // Radius of the sphere
    double Radius;

    /// Boundary of central box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_L(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of central box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_R(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of central box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_D(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of central box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_U(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of central box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_B(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of central box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_centr_F(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of right box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_right_L(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of right box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_right_R(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of  right box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_right_D(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of right box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_right_U(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of  right box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_right_B(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of  right box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_right_F(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_up_L(const unsigned& t,
                const Vector<double>& zeta,
                Vector<double>& f);

    /// Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_up_R(const unsigned& t,
                const Vector<double>& zeta,
                Vector<double>& f);

    /// Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_up_D(const unsigned& t,
                const Vector<double>& zeta,
                Vector<double>& f);

    /// Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_up_U(const unsigned& t,
                const Vector<double>& zeta,
                Vector<double>& f);

    /// Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_up_B(const unsigned& t,
                const Vector<double>& zeta,
                Vector<double>& f);

    /// Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_up_F(const unsigned& t,
                const Vector<double>& zeta,
                Vector<double>& f);

    /// Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_front_L(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_front_R(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_front_D(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_front_U(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_front_B(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of top left box macro element
    /// zeta \f$ \in [-1,1]^2 \f$
    void r_front_F(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);
  };

} // namespace oomph

#endif
