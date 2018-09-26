/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/



#include "SiconosKernel.hpp"


#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;


#define DEBUG_MESSAGES
#define DEBUG_STDOUT
#include "debug.h"



static int counter = 0;

class LinearElacticMaterial
{
  double _E;
  double _nu;
  double _rho;
  double _lambda;
  double _mu;

public:
  LinearElacticMaterial(){};

  
  LinearElacticMaterial(double E, double nu, double rho)
  {
    _E = E;
    _nu = nu;
    _rho = rho;

    _lambda = _E * _nu /( (1.0 +_nu) * (1.0 - 2* _nu));
    _mu = _E /(2.0*(1+_nu));
  };

  double mu()
  {
    return _mu;
  };
  double lambda()
  {
    return _lambda;
  };

};


class PunchLagrangianLinearTIDS : public LagrangianLinearTIDS
{

  /** The material of the elastic body */
  LinearElacticMaterial _mat;

  
  Mesh * _mesh;

  int _order;
  
  int _dim;

  FiniteElementCollection *_fec;
  
  FiniteElementSpace *_fespace;
  
  LinearForm * _b;

  BilinearForm *_a;

  Array<int> _ess_tdof_list;
  
  SP::SimpleMatrix create_from_mfem_matrix(SparseMatrix A)
  {
    int * Ai  = A.GetI();
    int * Aj  = A.GetJ();
    double * Ax  = A.GetData();
    int size = A.Size();

    int nnz = Ai[size];
    
    //SP::SimpleMatrix M(new SimpleMatrix(size,size,Siconos::SPARSE,nnz));
    SP::SimpleMatrix M(new SimpleMatrix(size,size));
   
    for (int row =0; row < size ; row++)
    {
      for (int k = Ai[row], end = Ai[row+1]; k < end; k++)
      {
        M->setValue(row, Aj[k], Ax[k]);
      }
    }

    //M->display();
    return M;
  };

  bool equal_mfem_matrix(SparseMatrix A, SparseMatrix B, double tol)
  {
    int * Ai  = A.GetI();
    int * Aj  = A.GetJ();
    double * Ax  = A.GetData();
    int Asize = A.Size();
    
    int * Bi  = B.GetI();
    int * Bj  = B.GetJ();
    double * Bx  = B.GetData();

    
    for (int row =0; row < Asize ; row++)
    {
      if (Ai[row] != Bi[row])
      {
        
        return false;
      }
      for (int k = Ai[row], end = Ai[row+1]; k < end; k++)
      {
        if (Aj[k] != Bj[k]) return false;
        if (fabs(Ax[k]- Bx[k]) > tol) return false;
      }
    }
    return true;
  };


  
public:
  PunchLagrangianLinearTIDS( const char *mesh_file, LinearElacticMaterial mat, int order =1)
  {
    _mat =mat;

    _order = order;

    // 2. Read the mesh from the given mesh file. We can handle triangular,
    //    quadrilateral, tetrahedral or hexahedral elements with the same code.
    _mesh = new Mesh(mesh_file, 1, 1);
    _dim = _mesh->Dimension();

    if (_mesh->bdr_attributes.Max() < 2)
    {
      cerr << "\nInput mesh should have at least two materials and "
           << "two boundary attributes! (See schematic in ex2.cpp)\n"
           << endl;
    }

       // 4. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
   //    largest number that gives a final mesh with no more than 5,000
   //    elements.
   // {
   //    int ref_levels =
   //       (int)floor(log(5000./mesh->GetNE())/log(2.)/dim);
   //    for (int l = 0; l < ref_levels; l++)
   //    {
   //       mesh->UniformRefinement();
   //    }
   // }
       // 5. Define a finite element space on the mesh. Here we use vector finite
   //    elements, i.e. dim copies of a scalar finite element space. The vector
   //    dimension is specified by the last argument of the FiniteElementSpace
   //    constructor. For NURBS meshes, we use the (degree elevated) NURBS space
   //    associated with the mesh nodes.

   if (_mesh->NURBSext)
   {
      _fec = NULL;
      _fespace = _mesh->GetNodes()->FESpace();
   }
   else
   {
      _fec = new H1_FECollection(order, _dim);
      _fespace = new FiniteElementSpace(_mesh, _fec, _dim);
   }
   // 13. For non-NURBS meshes, make the mesh curved based on the finite element
   //     space. This means that we define the mesh elements through a fespace
   //     based transformation of the reference element. This allows us to save
   //     the displaced mesh as a curved mesh when using high-order finite
   //     element displacement field. We assume that the initial mesh (read from
   //     the file) is not higher order curved mesh compared to the chosen FE
   //     space.
   if (!_mesh->NURBSext)
   {
      _mesh->SetNodalFESpace(_fespace);
   }

   _ndof = _fespace->GetTrueVSize();

   cout << "Number of finite element unknowns: " << _fespace->GetTrueVSize()
        << endl << "Assembling: " << flush;

   double position_init=0.0;
   double velocity_init=0.0;

   // -- Initial positions and velocities --
   SP::SiconosVector q0(new SiconosVector(_ndof,position_init));
   SP::SiconosVector v0(new SiconosVector(_ndof,velocity_init));

   _init(q0,v0);

  };

  void computeFEMForces()
  {
    // 7. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system. In this case, b_i equals the boundary integral
   //    of f*phi_i where f represents a "pull down" force on the Neumann part
   //    of the boundary and phi_i are the basis functions in the finite element
   //    fespace. The force is defined by the VectorArrayCoefficient object f,
   //    which is a vector of Coefficient objects. The fact that f is non-zero
   //    on boundary attribute 2 is indicated by the use of piece-wise constants
   //    coefficient for its last component.
    
   VectorArrayCoefficient f(_dim);
   for (int i = 0; i < _dim-1; i++)
   {
      f.Set(i, new ConstantCoefficient(0.0));
   }
   {
      Vector pull_force(_mesh->bdr_attributes.Max());
      pull_force = 0.0;
      pull_force(1) = -1.0;
      f.Set(_dim-1, new PWConstCoefficient(pull_force));
   }

   _b = new LinearForm(_fespace);
   _b->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(f));
   cout << "r.h.s. ... " << flush;
   _b->Assemble();

   DEBUG_EXPR(_b->Print(););

   _fExt.reset(new SiconosVector(_ndof));
   
   double * data=   _b->GetData();
   for (unsigned k =0; k < _ndof; k++)
   {
     _fExt->setValue(k,data[k]);
     
   }
   DEBUG_EXPR(_fExt->display(););
  };

  void computeFEMStiffness(){
    // 8. Define the solution vector x as a finite element grid function
    //    corresponding to fespace. Initialize x with initial guess of zero,
    //    which satisfies the boundary conditions.
   GridFunction x(_fespace);
   x = 0.0;
   
    // 9. Set up the bilinear form a(.,.) on the finite element space
    //    corresponding to the linear elasticity integrator with piece-wise
    //    constants coefficient lambda and mu.
   Vector lambda(_mesh->attributes.Max());
   lambda = 1.0;
   lambda(0) = lambda(1)*50;
   PWConstCoefficient lambda_func(lambda);
   Vector mu(_mesh->attributes.Max());
   mu = 1.0;
   mu(0) = mu(1)*50;
   PWConstCoefficient mu_func(mu);

   _a = new BilinearForm(_fespace);
   _a->AddDomainIntegrator(new ElasticityIntegrator(lambda_func,mu_func));

   // 10. Assemble the bilinear form and the corresponding linear system,
   //     applying any necessary transformations such as: eliminating boundary
   //     conditions, applying conforming constraints for non-conforming AMR,
   //     static condensation, etc.
   cout << "matrix ... " << flush;
   bool static_cond = false;

   if (static_cond) { _a->EnableStaticCondensation(); }
   _a->Assemble();
   _a->Finalize();

   SparseMatrix A;
   Vector B, X;
   A= _a->SpMat();

   _K = create_from_mfem_matrix(A);

   //Static test
   // computeFEMForces();
   // SparseMatrix A2;
   // _a->FormLinearSystem(_ess_tdof_list, x, *_b, A2, X, B);
   // // assert(equal_mfem_matrix(A, A2, 1e-15));
   // // cout << "done." << endl;
   // GSSmoother M_mfem(A2);
   // PCG(A2, M_mfem, B, X, 1, 500, 1e-8, 0.0);
   // std::cout << "\nB.Print()" <<std::endl;
   // B.Print();
   // std::cout << "\nX.Print()" <<std::endl;
   // X.Print();
   // std::cout << "\nb.Print()" <<std::endl;
   // _b->Print();
   // std::cout << "\nx.Print()" <<std::endl;
   // x.Print();
   // cout << "Size of linear system// : " << A.Height() << endl;

   // _a->RecoverFEMSolution(X, *_b, x);
  };
  
  void computeFEMMass(){

   BilinearForm *m = new BilinearForm(_fespace);
   m->AddDomainIntegrator(new VectorMassIntegrator());
   m->Assemble();
   // shift the eigenvalue corresponding to eliminated dofs to a large value
   // m->EliminateEssentialBCDiag(ess_bdr, numeric_limits<double>::min());
   m->Finalize();
   
   SparseMatrix M;
   M= m->SpMat();  
   //M.Print();

   // SparseMatrix M2;
   // m->FormLinearSystem(_ess_tdof_list, x, *b, M2, X, B);
   // assert(equal_mfem_matrix(M, M2, 1e-15));
   _mass = create_from_mfem_matrix(M);
   };

  void setDirichletBoundaryConditions(unsigned int boundary_id)
  {
    // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
    //    In this example, the boundary conditions are defined by marking only
    //    boundary attribute 1 from the mesh as essential and converting it to a
    //    list of true dofs.
   Array<int> ess_bdr(_mesh->bdr_attributes.Max());
   ess_bdr = 0;
   ess_bdr[0] = boundary_id;
   _fespace->GetEssentialTrueDofs(ess_bdr, _ess_tdof_list);
   std::cout  << "_ess_tdof_list.Print()" << std::endl;
   _ess_tdof_list.Print();
   std::cout  << "ess_bdr.Print()" << std::endl;
   ess_bdr.Print();

   int size = _ess_tdof_list.Size();
   std::cout  << "ess_tdof_list size() =" << size << std::endl;
   SP::IndexInt bdindex(new IndexInt(size));
   for(int k=0; k< size; k++)
   {
     (*bdindex)[k] = _ess_tdof_list[k];
   }
   SP::SiconosVector bdPrescribedVelocity(new SiconosVector(size));
   for(int k=0; k< size; k++) bdPrescribedVelocity->setValue(k,0.0);
   SP::BoundaryCondition bd (new BoundaryCondition(bdindex,bdPrescribedVelocity));
   setBoundaryConditions(bd);
  };

  
  void outputSolution()
  {
    GridFunction x(_fespace, _q[0]->getArray() );

    // Vector X;
    // _q[0]->display();
    // X.NewDataAndSize(_q[0]->getArray(), _ndof);
    // _a->RecoverFEMSolution(X, *_b, x);
    // std::cout << "\nX.Print()" <<std::endl;
    // X.Print();
    // std::cout << "\nx.Print()" <<std::endl;
    // x.Print();
    
    GridFunction *nodes = _mesh->GetNodes();
    //_mesh->Print();
    *nodes += x;
    x *= -1;


    
    ofstream mesh_ofs("displaced.mesh");
    mesh_ofs.precision(8);
    _mesh->Print(mesh_ofs);
    stringstream ss;
    ss << counter;
    string str = ss.str();
    ofstream sol_ofs("sol.gf" + str);
    sol_ofs.precision(8);
    x.Save(sol_ofs);
    counter ++;
    bool visualization= true;
    if (visualization)
    {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << *_mesh << x << flush;
    }
    
  };

  


};

TYPEDEF_SPTR(PunchLagrangianLinearTIDS)
