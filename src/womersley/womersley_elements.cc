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

// Non-inline functions for Womersley elements
#include "womersley_elements.h"


namespace oomph
{
  /// Default Womersley number
  template<unsigned DIM>
  double WomersleyEquations<DIM>::Default_ReSt_value = 0.0;


  //========================================================================
  /// Instantiation of static bool to suppress warning; initialise to false.
  //========================================================================
  bool TemplateFreeWomersleyMeshBase::Suppress_warning_about_unpinned_nst_dofs =
    false;


  //========================================================================
  /// Zero!
  //========================================================================
  double TemplateFreeWomersleyImpedanceTubeBase::Zero = 0.0;


  //======================================================================
  // Set the data for the number of Variables at each node
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  const unsigned QWomersleyElement<DIM, NNODE_1D>::Initial_Nvalue = 1;


  //======================================================================
  /// Compute element residual Vector and/or element Jacobian matrix
  ///
  /// flag=1: compute both
  /// flag=0: compute only residual Vector
  ///
  /// Pure version without hanging nodes
  //======================================================================
  template<unsigned DIM>
  void WomersleyEquations<DIM>::fill_in_generic_residual_contribution_womersley(
    Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
  {
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find the index at which the variable is stored
    unsigned u_nodal_index = u_index_womersley();

    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node, DIM), dtestdx(n_node, DIM);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Integers to hold the local equation and unknowns
    int local_eqn = 0, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++) s[i] = integral_pt()->knot(ipt, i);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_womersley(
        ipt, psi, dpsidx, test, dtestdx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Allocate memory for local quantities and initialise to zero
      double interpolated_u = 0.0;
      double dudt = 0.0;
      Vector<double> interpolated_x(DIM, 0.0);
      Vector<double> interpolated_dudx(DIM, 0.0);

      // Calculate function value and derivatives:
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Calculate the value at the nodes
        double u_value = raw_nodal_value(l, u_nodal_index);
        interpolated_u += u_value * psi(l);
        dudt += du_dt_womersley(l) * psi(l);
        // Loop over directions
        for (unsigned j = 0; j < DIM; j++)
        {
          interpolated_x[j] += raw_nodal_position(l, j) * psi(l);
          interpolated_dudx[j] += u_value * dpsidx(l, j);
        }
      }

      // Get pressure gradient
      //---------------------
      double dpdz;

      // If no pressure gradient has been set, use zero
      if (Pressure_gradient_data_pt == 0)
      {
        dpdz = 0.0;
      }
      else
      {
        // Get gressure gradient
        dpdz = Pressure_gradient_data_pt->value(0);
      }


      // Assemble residuals and Jacobian
      //--------------------------------

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        local_eqn = nodal_local_eqn(l, u_nodal_index);
        /*IF it's not a boundary condition*/
        if (local_eqn >= 0)
        {
          // Add dpdz term and time derivative
          residuals[local_eqn] += (re_st() * dudt + dpdz) * test(l) * W;

          // Laplace operator
          for (unsigned k = 0; k < DIM; k++)
          {
            residuals[local_eqn] += interpolated_dudx[k] * dtestdx(l, k) * W;
          }


          // Calculate the jacobian
          //-----------------------
          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              local_unknown = nodal_local_eqn(l2, u_nodal_index);
              // If at a non-zero degree of freedom add in the entry
              if (local_unknown >= 0)
              {
                // Mass matrix
                jacobian(local_eqn, local_unknown) +=
                  re_st() * test(l) * psi(l2) *
                  node_pt(l2)->time_stepper_pt()->weight(1, 0) * W;

                // Laplace operator & mesh velocity bit
                for (unsigned i = 0; i < DIM; i++)
                {
                  double tmp = dtestdx(l, i);
                  jacobian(local_eqn, local_unknown) += dpsidx(l2, i) * tmp * W;
                }
              }
            }

            // Derivative w.r.t. pressure gradient data (if it's
            // an unknown)
            if ((Pressure_gradient_data_pt != 0) &&
                (!Pressure_gradient_data_pt->is_pinned(0)))
            {
              local_unknown = external_local_eqn(0, 0);
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) += test(l) * W;
              }

              // Derivatives of the final eqn (volume flux constraint
              // w.r.t. to this unknown)
              unsigned final_local_eqn = external_local_eqn(0, 0);
              unsigned local_unknown = local_eqn; //  [from above just renaming
              //                                //   for clarity(!?)]
              jacobian(final_local_eqn, local_unknown) += psi(l) * W;
            }
          }
        }
      }


    } // End of loop over integration points
  }


  //======================================================================
  /// Compute volume flux through element
  //======================================================================
  template<unsigned DIM>
  double WomersleyEquations<DIM>::get_volume_flux()
  {
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find the index at which the variable is stored
    unsigned u_nodal_index = u_index_womersley();

    // Set up memory for the shape fcs
    Shape psi(n_node);
    DShape dpsidx(n_node, DIM);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Initialise flux
    double flux = 0.0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++) s[i] = integral_pt()->knot(ipt, i);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_eulerian_at_knot(ipt, psi, dpsidx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Allocate memory for local quantities and initialise to zero
      double interpolated_u = 0.0;

      // Calculate function value

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Calculate the value at the nodes (takes hanging node status
        // into account
        interpolated_u += nodal_value(l, u_nodal_index) * psi(l);
      }

      // Add to flux
      flux += interpolated_u * W;

    } // End of loop over integration points

    return flux;
  }


  //======================================================================
  /// Self-test:  Return 0 for OK
  //======================================================================
  template<unsigned DIM>
  unsigned WomersleyEquations<DIM>::self_test()
  {
    bool passed = true;

    // Check lower-level stuff
    if (FiniteElement::self_test() != 0)
    {
      passed = false;
    }

    // Return verdict
    if (passed)
    {
      return 0;
    }
    else
    {
      return 1;
    }
  }


  //======================================================================
  /// Output function:  x,y,z_out,0,0,u,0 to allow comparison
  /// against full Navier Stokes at n_nplot x n_plot points (2D)
  //======================================================================
  template<unsigned DIM>
  void WomersleyEquations<DIM>::output_3d(std::ostream& outfile,
                                          const unsigned& nplot,
                                          const double& z_out)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }
      outfile << z_out << " 0.0 0.0 ";
      outfile << interpolated_u_womersley(s);
      outfile << " 0.0 " << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }

  //======================================================================
  /// Output function:
  ///
  ///   x,y,u   or    x,y,z,u
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  template<unsigned DIM>
  void WomersleyEquations<DIM>::output(std::ostream& outfile,
                                       const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }
      outfile << interpolated_u_womersley(s) << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }


  //======================================================================
  /// C-style output function:
  ///
  ///   x,y,u   or    x,y,z,u
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  template<unsigned DIM>
  void WomersleyEquations<DIM>::output(FILE* file_pt, const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Tecplot header info
    fprintf(file_pt, "%s", tecplot_zone_string(nplot).c_str());

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      for (unsigned i = 0; i < DIM; i++)
      {
        fprintf(file_pt, "%g ", interpolated_x(s, i));
      }
      fprintf(file_pt, "%g \n", interpolated_u_womersley(s));
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(file_pt, nplot);
  }


  //======================================================================
  /// Output exact solution
  ///
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points.
  ///
  ///   x,y,u_exact    or    x,y,z,u_exact
  //======================================================================
  template<unsigned DIM>
  void WomersleyEquations<DIM>::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Exact solution Vector (here a scalar)
    Vector<double> exact_soln(1);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Get x position as Vector
      interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Output x,y,...,u_exact
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }
      outfile << exact_soln[0] << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }


  //======================================================================
  /// Output exact solution at time t
  ///
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points.
  ///
  ///   x,y,u_exact    or    x,y,z,u_exact
  //======================================================================
  template<unsigned DIM>
  void WomersleyEquations<DIM>::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    const double& time,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)

  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);


    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Exact solution Vector (here a scalar)
    Vector<double> exact_soln(1);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Get x position as Vector
      interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(time, x, exact_soln);

      // Output x,y,...,u_exact
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }
      outfile << exact_soln[0] << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }


  //======================================================================
  /// Validate against exact solution
  ///
  /// Solution is provided via function pointer.
  /// Plot error at a given number of plot points.
  ///
  //======================================================================
  template<unsigned DIM>
  void WomersleyEquations<DIM>::compute_error(
    std::ostream& outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    double& error,
    double& norm)
  {
    // Initialise
    error = 0.0;
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);

    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    Shape psi(n_node);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Tecplot header info
    outfile << "ZONE" << std::endl;

    // Exact solution Vector (here a scalar)
    Vector<double> exact_soln(1);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get x position as Vector
      interpolated_x(s, x);

      // Get FE function value
      double u_fe = interpolated_u_womersley(s);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Output x,y,...,error
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }
      outfile << exact_soln[0] << " " << exact_soln[0] - u_fe << std::endl;

      // Add to error and norm
      norm += exact_soln[0] * exact_soln[0] * W;
      error += (exact_soln[0] - u_fe) * (exact_soln[0] - u_fe) * W;
    }
  }


  //======================================================================
  /// Validate against exact solution at time t.
  ///
  /// Solution is provided via function pointer.
  /// Plot error at a given number of plot points.
  ///
  //======================================================================
  template<unsigned DIM>
  void WomersleyEquations<DIM>::compute_error(
    std::ostream& outfile,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
    const double& time,
    double& error,
    double& norm)

  {
    // Initialise
    error = 0.0;
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);


    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    Shape psi(n_node);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Tecplot
    outfile << "ZONE" << std::endl;

    // Exact solution Vector (here a scalar)
    Vector<double> exact_soln(1);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get x position as Vector
      interpolated_x(s, x);

      // Get FE function value
      double u_fe = interpolated_u_womersley(s);

      // Get exact solution at this point
      (*exact_soln_pt)(time, x, exact_soln);

      // Output x,y,...,error
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }
      outfile << exact_soln[0] << " " << exact_soln[0] - u_fe << std::endl;

      // Add to error and norm
      norm += exact_soln[0] * exact_soln[0] * W;
      error += (exact_soln[0] - u_fe) * (exact_soln[0] - u_fe) * W;
    }
  }


  //====================================================================
  // Force build of templates
  //====================================================================

  template class WomersleyEquations<1>;

  template class QWomersleyElement<1, 2>;
  template class QWomersleyElement<1, 3>;
  template class QWomersleyElement<1, 4>;

  template class WomersleyEquations<2>;

  template class QWomersleyElement<2, 2>;
  template class QWomersleyElement<2, 3>;
  template class QWomersleyElement<2, 4>;

  // template class QWomersleyElement<3,2>;
  // template class QWomersleyElement<3,3>;
  // template class QWomersleyElement<3,4>;

} // namespace oomph
