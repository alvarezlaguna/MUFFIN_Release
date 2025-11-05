//
//  TwoFluidIsothermal1DSourceTermMPI.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/01/2021.
//  Copyright © 2021 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef TwoFluidIsothermal1DSourceTermMPI_hpp
#define TwoFluidIsothermal1DSourceTermMPI_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include "SourceTerm.hpp"
#include "../LinearSolver/LinearSolver.hpp"
#include "../MeshData/Cell1D.hpp"
#include "../MeshData/MeshData.hpp"
//#include <mpi4py/mpi4py.h>
#include <petsc.h>
#include <petscksp.h>

#include <mpi.h>
#include <petscmat.h>
#include <petscvec.h>

using namespace Parameters;
using namespace std;

typedef vector< vector<double> > Matrix;

class TwoFluidIsothermal1DSourceTermMPI : public SourceTerm
{
public:
    TwoFluidIsothermal1DSourceTermMPI(string name) : SourceTerm(name), m_B(NBCELLS, 0.)
    {
        /* Taken from ex23.c in PETSc website */
        PetscErrorCode ierr;
        cout<<"ierr = "<<ierr<<"\n";
        PetscInitialize(NULL, NULL, NULL, NULL);
        int    rank, size;
        MPI_Comm_size(MPI_COMM_WORLD,&size);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        PetscPrintf(MPI_COMM_WORLD,"Number of processors = %d, rank = %d\n", size, rank);
        MatCreate(PETSC_COMM_WORLD,&m_A_matrix);
        int nLocal, nGlobal;
        nLocal = NBCELLS_MPI[rank];
        nGlobal = 0;
        for(int iP = 0; iP < size; ++iP){
            nGlobal += NBCELLS_MPI[iP];
        }

        MatSetSizes(m_A_matrix,nLocal,nLocal,nGlobal,nGlobal);
        MatSetFromOptions(m_A_matrix);
        MatSetUp(m_A_matrix);
        
        Vec            x, b, u;          /* approx solution, RHS, exact solution */
        PetscInt       i, col[3], its, rstart, rend;
        PetscScalar    neg_one = -1.0,one = 1.0,value[3];
        // These vars should go as members
        PetscReal      norm,tol=1.e-11;  /* norm of solution error */
        
        PetscPrintf(MPI_COMM_WORLD,"Create vectors\n");
        ierr = VecCreate(PETSC_COMM_WORLD,&m_x_vec);
        ierr = VecSetSizes(m_x_vec,nLocal,nGlobal);
        ierr = VecSetFromOptions(m_x_vec);
        ierr = VecDuplicate(m_x_vec,&m_b_vec);
        ierr = VecDuplicate(m_x_vec,&u);
        
        ierr = VecGetOwnershipRange(m_x_vec,&rstart,&rend);
        ierr = VecGetLocalSize(m_x_vec,&nLocal);

        PetscPrintf(MPI_COMM_WORLD,"After create vectors\n");

        /*
         Assemble matrix.
         
         The linear system is distributed across the processors by
         chunks of contiguous rows, which correspond to contiguous
         sections of the mesh on which the problem is discretized.
         For matrix assembly, each processor contributes entries for
         the part that it owns locally.
         */
        
        
        if (!rstart) {
            rstart = 1;
            i      = 0; col[0] = 0; col[1] = 1; value[0] = -2.0; value[1] = 1.0;
            ierr   = MatSetValues(m_A_matrix,1,&i,2,col,value,INSERT_VALUES);
        }
        if (rend == nGlobal) {
            rend = nGlobal-1;
            i    = nGlobal-1; col[0] = nGlobal-2; col[1] = nGlobal-1; value[0] = 1.0; value[1] = -2.0;
            ierr = MatSetValues(m_A_matrix,1,&i,2,col,value,INSERT_VALUES);
        }
        
        /* Set entries corresponding to the mesh interior */
        value[0] = 1.0; value[1] = -2.0; value[2] = 1.0;
        for (i=rstart; i<rend; i++) {
            col[0] = i-1; col[1] = i; col[2] = i+1;
            ierr   = MatSetValues(m_A_matrix,1,&i,3,col,value,INSERT_VALUES);
        }
        
        /* Assemble the matrix */
        ierr = MatAssemblyBegin(m_A_matrix,MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(m_A_matrix,MAT_FINAL_ASSEMBLY);
        PetscPrintf(MPI_COMM_WORLD,"After ensembled matix\n");

        // Visualized the matrix
        //MatView(m_A_matrix,PETSC_VIEWER_STDOUT_WORLD);

        
        /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Create the linear solver and set various options
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
        /*
         Create linear solver context
         */
        ierr = KSPCreate(PETSC_COMM_WORLD,&m_ksp);
        
        /*
         Set operators. Here the matrix that defines the linear system
         also serves as the preconditioning matrix.
         */
        ierr = KSPSetOperators(m_ksp,m_A_matrix,m_A_matrix);
        
        /*
         Set linear solver defaults for this problem (optional).
         - By extracting the KSP and PC contexts from the KSP context,
         we can then directly call any KSP and PC routines to set
         various options.
         - The following four statements are optional; all of these
         parameters could alternatively be specified at runtime via
         KSPSetFromOptions();
         */
        ierr = KSPGetPC(m_ksp,&m_pc);
        ierr = PCSetType(m_pc,PCGAMG); // Other types of preconditioner PCSOR, PCJACOBI
        ierr = KSPSetTolerances(m_ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
        // Default options
        //ierr = PCSetType(m_pc,PCJACOBI);
        //ierr = KSPSetTolerances(m_ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
        
        /*
         Set runtime options, e.g.,
         -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
         These options will override those specified above as long as
         KSPSetFromOptions() is called _after_ any other customization
         routines.
         */
        //ierr = KSPSetFromOptions(ksp);
        

        
        /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Check solution and clean up
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
        /*
         Check the error
         */
//        ierr = VecAXPY(x,neg_one,u);
//        ierr = VecNorm(x,NORM_2,&norm);
//        ierr = KSPGetIterationNumber(m_ksp,&its);
//        ierr = MatMult(m_A_matrix,m_x_vec,m_b_vec);
//        VecView(m_b_vec,PETSC_VIEWER_STDOUT_WORLD);
//        VecView(m_x_vec,PETSC_VIEWER_STDOUT_WORLD);
//        if (norm > tol) {
//            ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",norm,its);
//        }

    }
    ~TwoFluidIsothermal1DSourceTermMPI(){PetscFinalize();}
    void solveLinearSystem(Mat& m_A,Vec& m_x,Vec& m_B);
    void setup();
    virtual void computeSource();
    double computeIonizationConstant();

private:
    /// TODO: write here the matrix and the vectors for the linear system
    unique_ptr<LinearSolver> m_linearSolver;
    //Matrix m_A;         // Matrix of the linear solver is a PETSc matrix
    Mat m_A_matrix;     // Matrix for Petsc
    Vec m_x_vec;        // Solution of the linear solver
    Vec m_b_vec;        // RHD of the linear solver
    KSP m_ksp;          // Krilov method
    PC  m_pc;           // Preconditioner
    vector<double> m_x; // Solution of the linear solver
    vector<double> m_B; // rhs of the linear solver
    
};

#endif /* TwoFluidIsothermal1DSourceTermMPI_hpp */
