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
#ifndef OOMPH_RECTANGLE_WITH_HOLE_DOMAIN_HEADER
#define OOMPH_RECTANGLE_WITH_HOLE_DOMAIN_HEADER

#include "generic/element_factory.h"

// Generic includes
#include "generic/quadtree.h"
#include "generic/geom_objects.h"
#include "generic/macro_element.h"
#include "generic/domain.h"

namespace oomph
{
  //===========================================================
  /// Rectangular domain with circular whole
  //===========================================================
  class RectangleWithHoleDomain : public Domain
  {
  public:
    /// Constructor. Pass pointer to geometric object that
    /// represents the cylinder, the length of the (square) domain.
    /// The GeomObject must be parametrised such that
    /// \f$\zeta \in [0,2\pi]\f$ sweeps around the circumference
    /// in anticlockwise direction.
    RectangleWithHoleDomain(GeomObject* cylinder_pt, const double& length);

    /// Destructor: Empty; cleanup done in base class
    ~RectangleWithHoleDomain() {}

    /// Helper function to interpolate linearly between the
    /// "right" and "left" points; \f$ s \in [-1,1] \f$
    void linear_interpolate(Vector<double> left,
                            Vector<double> right,
                            const double& s,
                            Vector<double>& f);


    /// Parametrisation of macro element boundaries: f(s) is the position
    /// vector to macro-element m's boundary in the specified direction
    /// [N/S/E/W] at the specfied discrete time level (time=0: present; time>0:
    /// previous)
    void macro_element_boundary(const unsigned& time,
                                const unsigned& m,
                                const unsigned& direction,
                                const Vector<double>& s,
                                Vector<double>& f);

  private:
    /// Lower left corner of rectangle
    Vector<double> Lower_left;

    /// Lower right corner of rectangle
    Vector<double> Lower_right;

    /// Where the "radial" line from circle meets lower boundary on left
    Vector<double> Lower_mid_left;

    /// Where the "radial" line from circle meets lower boundary on right
    Vector<double> Lower_mid_right;

    /// Upper left corner of rectangle
    Vector<double> Upper_left;

    /// Upper right corner of rectangle
    Vector<double> Upper_right;

    /// Where the "radial" line from circle meets upper boundary on left
    Vector<double> Upper_mid_left;

    /// Where the "radial" line from circle meets upper boundary on right
    Vector<double> Upper_mid_right;

    /// Pointer to geometric object that represents the central cylinder
    GeomObject* Cylinder_pt;
  };

} // namespace oomph
#endif
