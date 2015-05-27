#include <stdio.h>
#include <mpi.h>
#include <sstream>
#include <vector>
#include "sparse_matrix.h"
#include "mpimatrix.h"
#include "generator.h"
#include "dense_matrix.h"
#include <ctime>

using namespace std;

#define MAIN_DEBUG 1

void printHelp(int rank)
{
    if (rank != 0) return;
    printf("This program is using MPI to do operations on sparse matrices\n");
    printf("Usage:\n");
    printf("./program <matrix_file_path> <matrix_file_path> <operation>\n");
}

void printLU(int rank, sparse_matrix &L, sparse_matrix &U, const char* name)
{
    if (rank != 0) return;
    printf("%s L:\n", name);
    L.printDense();
    printf("%s U:\n", name);
    U.printDense();
}

void printM(int rank, sparse_matrix &M, const char* name)
{
    if (rank != 0) return;
    printf("%s:\n", name);
    M.printDense();
}

void genRandomM(int rank, int size)
{
    if (rank != 0) return;
    auto randM = Generator(rank, size).GenerateRandomMatrix(10, 10, 25, column_wise);
    printf("Random sparse matrix: \n");
    randM.printDense();
    auto densM = dense_matrix(randM);
    printf("dense matrix: \n");
    densM.printDense();
}

int main (int argc, char *argv[])
{
    int rank, size;
    std::clock_t start;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    MpiMatrixHelper mpi_helper(rank, size);

    auto m1 = mpi_helper.load("big", sparse, column_wise, 1);
    auto m2 = mpi_helper.load("big", sparse, row_wise, 1);

    sparse_matrix mult_result = mpi_helper.mul(m1, m2);
    sparse_matrix add_result = mpi_helper.add(m1, m2);

    sparse_matrix L, U, CL, CU;

    mpi_helper.LU(m1, CL, CU);
    auto CLU = mpi_helper.mul(CL, CU);
    if(rank == 0)
    {
        if (CLU == m1) printf("LU SUCCESS!!\n");
        else
        {
            printf("LU FAILURE!! Result:\n");
            CLU.printDense();
        }
    }

    mpi_helper.ILU(m1, L, U);
    if(rank == 0) printf("ILU DONE!\n");

    auto L_inv = mpi_helper.Inverse(L);
    auto I = mpi_helper.mul(L, L_inv);
    if(rank == 0)
    {
        if (I == sparse_matrix::identity(L.getWidth()))
            printf("INVERSE SUCCESS!!\n");
        else printf("INVERSE FAILURE!!\n");
    }

    // CG non-mpi test
    int n = m1.getHeight();
    if(rank == 0) printf("n = %d\n", n);
    sparse_vector b(n, column_wise);
    for(int i=0; i<n; i++)
        b[i] = i;
    start = std::clock();
    auto cg_result = mpi_helper.CG_ILU_PRECONDITIONED_WITH_ILU(m1, b);
    auto time = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    printf("TIME = %f\n", time);

    // Check
    if(rank == 0)
    {
        sparse_vector wyn = m1 * cg_result;

        if (wyn == b) printf("CG SUCCESS!!\n");
        else
        {
            printf("CG FAILURE!! Result is:\n");

        }
        wyn.print();
    }

    MPI_Finalize();
    return 0;
}
