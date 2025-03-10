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
#ifndef OOMPH_ANNULAR_DOMAIN_HEADER
#define OOMPH_ANNULAR_DOMAIN_HEADER

#include "generic/element_factory.h"

// Generic oomph-lib includes
#include "generic/quadtree.h"
#include "generic/domain.h"
#include "generic/geom_objects.h"

namespace oomph
{
  //=================================================================
  /// Annular domain
  //=================================================================
  class AnnularDomain : public Domain
  {
  public:
    /// Constructor: Specify azimuthal fraction (1.0 is 360 degrees)
    /// number of macro elements in azimuthal and radial direction,
    /// inner radius and thickness. Rotate mesh by angle phi.
    AnnularDomain(const double& azimuthal_fraction,
                  const unsigned& ntheta,
                  const unsigned& nr,
                  const double& a,
                  const double& h,
                  const double& phi);

    /// Broken copy constructor
    AnnularDomain(const AnnularDomain&) = delete;

    /// Broken assignment operator
    void operator=(const AnnularDomain&) = delete;

    /// Destructor: Empty; cleanup done in base class
    ~AnnularDomain() {}

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
    /// Azimuthal fraction
    double Azimuthal_fraction;

    /// Inner radius
    double Inner_radius;

    /// Thickness
    double Thickness;

    /// Number of macro elements in azimuthal direction
    unsigned Ntheta;

    /// Number of macro elements in radial direction
    unsigned Nr;

    /// Rotation angle
    double Phi;
  };

} // namespace oomph

#endif
