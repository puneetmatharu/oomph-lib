//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
// Driver for adaptive 2D rectangular driven cavity. Solved with black
// box adaptation, using Taylor Hood and Crouzeix Raviart elements.

// Generic oomph-lib header
#include "generic.h"

// Navier Stokes headers
#include "navier_stokes.h"

// The mesh
#include "meshes/simple_rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//==start_of_namespace===================================================
/// Namespace for physical parameters
//=======================================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re=100;
} // end_of_namespace



//==============================start_of_mesh======================
/// Refineable equivalent of the SimpleRectangularQuadMesh.
/// Refinement is performed by the QuadTree-based procedures
/// implemented in the RefineableQuadMesh base class.
//=================================================================
template<class ELEMENT>
class SimpleRefineableRectangularQuadMesh : 
 public virtual SimpleRectangularQuadMesh<ELEMENT>,  
 public RefineableQuadMesh<ELEMENT>
{ 

public: 

 ///  Pass number of elements in the horizontal 
 /// and vertical directions, and the corresponding dimensions.
 /// Timestepper defaults to Static.
 SimpleRefineableRectangularQuadMesh(const unsigned &Nx,
                                     const unsigned &Ny, 
                                     const double &Lx, const double &Ly,
                                     TimeStepper* time_stepper_pt=
                                     &Mesh::Default_TimeStepper) :
  SimpleRectangularQuadMesh<ELEMENT>(Nx,Ny,Lx,Ly,time_stepper_pt)
  {
   // Nodal positions etc. were created in constructor for
   // SimpleRectangularQuadMesh<...> --> We only need to set up 
   // adaptivity information: Associate finite elements with their 
   // QuadTrees and plant them in a QuadTreeForest:
   this->setup_quadtree_forest();
  } // end of constructor
 

 /// Destructor: Empty
 virtual ~SimpleRefineableRectangularQuadMesh() {}

}; // end of mesh



//==start_of_problem_class============================================
/// Driven cavity problem in rectangular domain, templated
/// by element type. 
//====================================================================
template<class ELEMENT>
class PRefineableDrivenCavityProblem : public Problem
{

public:

 /// Constructor
 PRefineableDrivenCavityProblem();

 /// Destructor: Empty
 ~PRefineableDrivenCavityProblem() {}

 /// Update the after solve (empty)
 void actions_after_newton_solve() {}

 /// Update the problem specs before solve. 
 /// (Re-)set velocity boundary conditions just to be on the safe side...
 void actions_before_newton_solve()
  { 
  // Setup tangential flow along boundary 0:
  unsigned ibound=0; 
  unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    // Tangential flow
    unsigned i=0;
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,1.0);
    // No penetration
    i=1;
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
   }
  
  // Overwrite with no flow along all other boundaries
  unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned ibound=1;ibound<num_bound;ibound++)
   {
    unsigned num_nod= mesh_pt()->nboundary_node(ibound);
    //cout << num_nod << " nodes on boundary " << ibound << endl;
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      for (unsigned i=0;i<2;i++)
       {
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
       }
     }
   }
  } // end_of_actions_before_newton_solve


 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableNavierStokesEquations<2>::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());
    
    // Pin redundant pressure dofs
   RefineableNavierStokesEquations<2>::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   // Now set the first pressure dof in the first element to 0.0
   fix_pressure(0,0,0.0);
  } // end_of_actions_after_adapt
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
private:

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to proper element and fix pressure
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure

}; // end_of_problem_class



//==start_of_constructor==================================================
/// Constructor for RefineableDrivenCavity problem 
///
//========================================================================
template<class ELEMENT>
PRefineableDrivenCavityProblem<ELEMENT>::PRefineableDrivenCavityProblem()
{ 

 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=4;

 // # of elements in y-direction
 unsigned n_y=4;
 
 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;

 // Build and assign mesh
 Problem::mesh_pt() = 
  new SimpleRefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 dynamic_cast<SimpleRefineableRectangularQuadMesh<ELEMENT>*>(mesh_pt())->
  spatial_error_estimator_pt()=error_estimator_pt;
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here: All boundaries are Dirichlet boundaries.
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over values (u and v velocities)
     for (unsigned i=0;i<2;i++)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
      }
    }
  } // end loop over boundaries

 //Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to Reynolds
 // number
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
  } // end loop over elements
 
 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 
 // Now set the first pressure dof in the first element to 0.0
 fix_pressure(0,0,0.0);
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} // end_of_constructor



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void PRefineableDrivenCavityProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=3; 


 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
 
} // end_of_doc_solution



//==start_of_main======================================================
/// Driver for RefineableDrivenCavity test problem 
//=====================================================================
int main(int argc, char** argv)
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 // Set output directory
 DocInfo doc_info;
 doc_info.set_directory("RESLT");
 
 

 // Solve problem with hp-refineable Crouzeix Raviart elements
 //-----------------------------------------------------------
 // Build problem
 PRefineableDrivenCavityProblem<PRefineableQCrouzeixRaviartElement<2> > problem;
 
 cout << "\n\n\nProblem self-test ";
 if (problem.self_test()==0) 
  {
   cout << "passed: Problem can be solved." << std::endl;
  }
 else 
  {
   throw OomphLibError("Self test failed",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

 if(CommandLineArgs::Argc==1)
  {
   problem.p_adapt();
   problem.newton_solve();
   doc_info.number()=1;
   problem.doc_solution(doc_info);
 
   problem.adapt();
   problem.newton_solve();
   doc_info.number()=2;
   problem.doc_solution(doc_info);
 
   problem.p_adapt();
   problem.newton_solve();
   doc_info.number()=3;
   problem.doc_solution(doc_info);
 
   problem.adapt();
   problem.newton_solve();
   doc_info.number()=4;
   problem.doc_solution(doc_info);
 
   problem.p_adapt();
   problem.newton_solve();
   doc_info.number()=5;
   problem.doc_solution(doc_info);
 
   problem.adapt();
   problem.newton_solve();
   doc_info.number()=6;
   problem.doc_solution(doc_info);
 
   problem.p_adapt();
   problem.newton_solve();
   doc_info.number()=7;
   problem.doc_solution(doc_info);
 
   problem.adapt();
   problem.newton_solve();
   doc_info.number()=8;
   problem.doc_solution(doc_info);
  }
 // Validation run
 else
  {
   problem.refine_uniformly();
   problem.refine_uniformly();
   problem.newton_solve();
   doc_info.number()=1;
   problem.doc_solution(doc_info);
 
   problem.p_adapt();
   problem.newton_solve();
   doc_info.number()=2;
   problem.doc_solution(doc_info);
 
   problem.adapt();
   problem.newton_solve();
   doc_info.number()=3;
   problem.doc_solution(doc_info);
 
   problem.p_adapt();
   problem.newton_solve();
   doc_info.number()=4;
   problem.doc_solution(doc_info);
  }
 
   // Count hanging nodes
 //cout << "Hanging nodes:" << endl;
 unsigned num_hang=0;
 for (unsigned n=0; n<problem.mesh_pt()->nnode(); n++)
  {
   if (problem.mesh_pt()->node_pt(n)->is_hanging())
    {
     /*
     cout << "  node " << n << " is hanging... at  ("
          << problem.mesh_pt()->node_pt(n)->x(0) << ", "
          << problem.mesh_pt()->node_pt(n)->x(1) << ")" << endl;
     HangInfo* hang_pt = problem.mesh_pt()->node_pt(n)->hanging_pt();
     cout << "    Nmaster = " << hang_pt->nmaster() << endl;
     double totweight = 0.0;
     for (unsigned nm=0; nm<hang_pt->nmaster(); nm++)
      {
       cout << "    master node:  x = (" << hang_pt->master_node_pt(nm)->x(0)
            << ", " << hang_pt->master_node_pt(nm)->x(1) << ")  w = "
            << hang_pt->master_weight(nm) << endl;
       totweight += hang_pt->master_weight(nm);
      }
     cout << "    Total weights = " << totweight << endl;
     */
     num_hang++;
    }
  }
 cout << "There were "<<num_hang<<" hanging nodes." << endl;
 
 // Step number
 doc_info.number()=0;

 //Output solution
 problem.doc_solution(doc_info);
 
} // end_of_main











