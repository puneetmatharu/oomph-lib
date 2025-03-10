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


#include "full_circle_domain.h"

namespace oomph
{

  //=================================================================
  /// Vector representation of the  imacro-th macro element
  /// boundary idirect (N/S/W/E) at time level t
  /// (t=0: present; t>0: previous): f(s)
  //=================================================================
  void FullCircleDomain::macro_element_boundary(const unsigned& t,
                                                const unsigned& imacro,
                                                const unsigned& idirect,
                                                const Vector<double>& s,
                                                Vector<double>& f)
  {
    using namespace QuadTreeNames;

    // Get the coordinates of the corners of the box
    Vector<Vector<double>> Box(4);
    // Get the corresponding coordinates on the boundary
    Vector<Vector<double>> Wall(4);

    // Storage for position in the area
    Vector<double> zeta(2);

    // Loop over the angles
    for (unsigned j = 0; j < 4; j++)
    {
      // Set up the storage
      Box[j].resize(2);
      Wall[j].resize(2);

      // Set the other values of zeta
      zeta[0] = Radius_box[j];
      zeta[1] = Theta_positions[j];
      // Now get the values
      Area_pt->position(t, zeta, Box[j]);

      // Now get the position on the boundary
      zeta[0] = 1.0;
      Area_pt->position(t, zeta, Wall[j]);
    }

    // Define pi
    const double pi = MathematicalConstants::Pi;

    // Which macro element?
    // --------------------
    switch (imacro)
    {
        // Macro element 0: Central box
      case 0:

        // Choose a direction
        switch (idirect)
        {
          case N:
            // Linearly interpolate between the two corners of the box
            lin_interpolate(Box[3], Box[2], s[0], f);
            break;

          case S:
            // Linearly interpolate between the two corners of the box
            lin_interpolate(Box[0], Box[1], s[0], f);

          case W:
            // Linearly interpolate between the two corners of the box
            lin_interpolate(Box[0], Box[3], s[0], f);
            break;

          case E:
            // Linearly interpolate between the two corners of the box
            lin_interpolate(Box[1], Box[2], s[0], f);
            break;

          default:
            std::ostringstream error_stream;
            error_stream << "idirect is " << idirect << " not one of N, S, E, W"
                         << std::endl;

            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
            break;
        }

        break;

        // Macro element 1: Bottom
      case 1:

        // Choose a direction
        switch (idirect)
        {
          case N:
            // Linearly interpolate between the two corners of the box
            lin_interpolate(Box[0], Box[1], s[0], f);
            break;

          case S:
            // Get the position on the wall by linearly interpolating in theta
            zeta[0] = 1.0;
            zeta[1] =
              Theta_positions[0] +
              (Theta_positions[1] - Theta_positions[0]) * 0.5 * (s[0] + 1.0);

            Area_pt->position(t, zeta, f);
            break;

          case W:
            // Now linearly interpolate between the wall and the box
            lin_interpolate(Wall[0], Box[0], s[0], f);
            break;

          case E:
            // Linearly interpolate between the wall and the box
            lin_interpolate(Wall[1], Box[1], s[0], f);
            break;

          default:

            std::ostringstream error_stream;
            error_stream << "idirect is " << idirect << " not one of N, S, E, W"
                         << std::endl;

            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
            break;
        }

        break;

      // Macro element 2:Right
      case 2:

        // Which direction?
        switch (idirect)
        {
          case N:
            // Linearly interpolate between the box and the wall
            lin_interpolate(Box[2], Wall[2], s[0], f);
            break;

          case S:
            // Linearly interpolate between the box and the wall
            lin_interpolate(Box[1], Wall[1], s[0], f);
            break;

          case W:
            // Linearly interpolate between the two corners of the box
            lin_interpolate(Box[1], Box[2], s[0], f);
            break;

          case E:
            // Get the position on the wall by linearly interpolating in theta
            zeta[0] = 1.0;
            zeta[1] =
              Theta_positions[1] +
              (Theta_positions[2] - Theta_positions[1]) * 0.5 * (s[0] + 1.0);

            Area_pt->position(t, zeta, f);
            break;

          default:
            std::ostringstream error_stream;
            error_stream << "idirect is " << idirect << " not one of N, S, W, E"
                         << std::endl;

            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

        // Macro element 3: Top
      case 3:

        // Which direction?
        switch (idirect)
        {
          case N:
            // Get the position on the wall by linearly interpolating in theta
            zeta[0] = 1.0;
            zeta[1] =
              Theta_positions[3] +
              (Theta_positions[2] - Theta_positions[3]) * 0.5 * (s[0] + 1.0);

            Area_pt->position(t, zeta, f);
            break;

          case S:
            // Linearly interpolate between the two corners of the box
            lin_interpolate(Box[3], Box[1], s[0], f);
            break;

          case W:
            // Linearly interpolate between the box and the wall
            lin_interpolate(Box[3], Wall[3], s[0], f);
            break;

          case E:
            // Linearly interpolate between the box and the wall
            lin_interpolate(Box[2], Wall[2], s[0], f);
            break;

          default:
            std::ostringstream error_stream;
            error_stream << "idirect is " << idirect << " not one of N, S, E, W"
                         << std::endl;

            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }

        break;

        // Macro element 4: Left
      case 4:

        // Which direction?
        switch (idirect)
        {
          case N:
            // Linearly interpolate between the wall and the box
            lin_interpolate(Wall[3], Box[3], s[0], f);
            break;

          case S:
            // Linearly interpolate between the wall and the box
            lin_interpolate(Wall[0], Box[0], s[0], f);
            break;

          case W:
            // Entirely on the wall, Need to be careful
            // because this is the point at which theta passes through the
            //"branch cut"
            zeta[0] = 1.0;
            zeta[1] = Theta_positions[0] + 2.0 * pi +
                      (Theta_positions[3] - (Theta_positions[0] + 2.0 * pi)) *
                        0.5 * (s[0] + 1.0);

            Area_pt->position(t, zeta, f);
            break;

          case E:
            // Linearly interpolate between the two corners of the box
            lin_interpolate(Box[0], Box[3], s[0], f);
            break;

          default:
            std::ostringstream error_stream;
            error_stream << "idirect is " << idirect << " not one of N, S, W, E"
                         << std::endl;

            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }
        break;

      default:
        // Error
        std::ostringstream error_stream;
        error_stream << "Wrong imacro " << imacro << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }

} // namespace oomph
