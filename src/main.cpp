#include <stdio.h>
#include <mpi.h>
#include <sstream>
#include <vector>
#include "sparse_matrix.h"
#include "mpimatrix.h"
#include "generator.h"
#include "dense_matrix.h"

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
