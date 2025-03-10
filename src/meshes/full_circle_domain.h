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

#ifndef OOMPH_FULL_CIRCLE_DOMAIN_HEADER
#define OOMPH_FULL_CIRCLE_DOMAIN_HEADER

#include "generic/element_factory.h"

// Generic oomph-lib includes
#include "generic/quadtree.h"
#include "generic/domain.h"
#include "generic/geom_objects.h"

namespace oomph
{
  //=================================================================
  /// Topologically circular domain, e.g. a tube cross section.
  /// The entire domain must be defined by a GeomObject with the
  /// following convention: zeta[0] is the radial coordinate and
  /// zeta[1] is the theta coordinate around the cross-sectin.
  /// The outer boundary must lie at zeta[0] = 1.
  ///
  /// The domain is
  /// parametrised by five macro elements (a central box surrounded by
  /// four curved elements). The labelling of the macro elements is shown
  /// below.
  ///
  ///                       ----------------------------
  ///                       |\                        /|
  ///                       | \       Macro          / |
  ///                       |  3      Element 3     2  |
  ///                       |   \                  /   |
  ///                       |    \----------------/    |
  ///                       |     |               |    |
  ///                       | 4   |    Macro      |    |
  ///                       |     |    Element 0  | 2  |
  ///                       |     |               |    |
  ///                       |     -----------------    |
  ///                       |    /                 \   |
  ///                       |   0     Macro         1  |
  ///                       |  /      Element 1      \ |
  ///                       | /                       \|
  ///                       |/-------------------------|
  ///
  ///
  //=================================================================
  class FullCircleDomain : public Domain
  {
  public:
    /// Constructor: Pass geometric  object; the theta locations
    /// marking the division between
    /// the elements of the outer ring, labelled from the lower left to the
    /// upper left in order, theta should be in the range \f$-\pi\f$ to
    /// \f$\pi\f$; and  the corresponding fractions of the
    /// radius at which the central box is to be placed.
    FullCircleDomain(GeomObject* area_geom_object_pt,
                     const Vector<double>& theta_positions,
                     const Vector<double>& radius_box)
      : Theta_positions(theta_positions),
        Radius_box(radius_box),
        Area_pt(area_geom_object_pt)
    {
      // There are five macro elements
      const unsigned n_macro = 5;
      Macro_element_pt.resize(n_macro);

      // Create the macro elements
      for (unsigned i = 0; i < n_macro; i++)
      {
        Macro_element_pt[i] = new QMacroElement<2>(this, i);
      }
    }

    /// Broken copy constructor
    FullCircleDomain(const FullCircleDomain&) = delete;

    /// Broken assignment operator
    void operator=(const FullCircleDomain&) = delete;

    /// Destructor: Empty; cleanup done in base class
    ~FullCircleDomain() {}


    /// Vector representation of the  i_macro-th macro element
    /// boundary i_direct (N/S/W/E) at time level t
    /// (t=0: present; t>0: previous):
    /// f(s).
    void macro_element_boundary(const unsigned& t,
                                const unsigned& i_macro,
                                const unsigned& i_direct,
                                const Vector<double>& s,
                                Vector<double>& f);

  private:
    /// Storage for the dividing lines on the boundary
    /// starting from the lower left and proceeding anticlockwise to
    /// the upper left
    Vector<double> Theta_positions;

    /// Storage for the fraction of the radius at which the central box
    /// should be located corresponding to the chosen values of theta.
    Vector<double> Radius_box;

    /// Pointer to geometric object that represents the domain
    GeomObject* Area_pt;

    ///  A very little linear interpolation helper.
    /// Interpolate from the low
    /// point to the high point using the coordinate s, which is
    /// assumed to run from -1 to 1.
    void lin_interpolate(const Vector<double>& low,
                         const Vector<double>& high,
                         const double& s,
                         Vector<double>& f)
    {
      // Loop over all coordinates (we are working in 2D)
      for (unsigned i = 0; i < 2; i++)
      {
        f[i] = low[i] + (high[i] - low[i]) * 0.5 * (s + 1.0);
      }
    }
  };

} // namespace oomph

#endif
