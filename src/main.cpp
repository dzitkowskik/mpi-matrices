#include <stdio.h>
#include <mpi.h>
#include <sstream>
#include <vector>
#include "sparse_matrix.h"
#include "mpimatrix.h"
#include "generator.h"
#include "dense_matrix.h"
#include "cg.h"

using namespace std;

void printHelp()
{
    printf("This program is using MPI to do operations on sparse matrices\n");
    printf("Usage:\n");
    printf("./program <matrix_file_path> <matrix_file_path> <operation>\n");
}

void printLU(MpiMatrix &L, MpiMatrix &U, const char* name)
{
    printf("%s L:\n", name);
    L.print();
    printf("%s U:\n", name);
    U.print();
}

void printM(MpiMatrix &M, const char* name)
{
    printf("%s:\n", name);
    M.print();
}

void genRandomM(int rank, int size)
{
    auto randM = Generator(rank, size).GenerateRandomMatrix(10, 10, 25, column_wise);
    printf("Random sparse matrix: \n");
    randM.print();
    auto densM = dense_matrix(randM.matrix);
    printf("dense matrix: \n");
    densM.printDense();
}

int main (int argc, char *argv[])
{
    int rank, size;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    MpiMatrix m1, m2;
    if (rank == 0)
    {
//        genRandomM(rank, size);
        printHelp();
    }

    m1 = MpiMatrix::load("/tmp/mat", rank, size, sparse);
    m2 = MpiMatrix::load("/tmp/mat", rank, size, sparse);

    MpiMatrix mult_result = m1*m2;
    MpiMatrix add_result = m1+m2;

    MpiMatrix L, U, CL, CU;
    m1.LU(CL, CU);
    m1.ILU(L, U);

    // CG non-mpi test
    if (rank == 0)
    {
        m1.matrix.transpose();
        sparse_vector x(4, column_wise);
        x[0] = x[1] = x[2] = x[3] = 0;
        sparse_vector b(4, column_wise);
        b[0] = 1; b[1] = 2; b[2] = 3; b[3] = 4;
        cg(m1.matrix, x, b, CL.matrix, CU.matrix);
        printf("RESULT X =\n");
        x.print();
    }

    if (rank == 0)
    {
        printM(m1, "Matrix 1");
        printM(m2, "Matrix 2");
        printM(mult_result, "Multiplication");
        printM(add_result, "Sum");
        printLU(L, U, "ILU ");
        printLU(CL, CU, "LU ");
    }
    MPI_Finalize();
    return 0;
}
