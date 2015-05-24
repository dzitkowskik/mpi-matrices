#include "stdio.h"
#include "mpi.h"
#include "../mpimatrix.h"
#include "../generator.h"
#include "../dense_matrix.h"
#include <ctime>
#include <unistd.h>

#define RANDOM_TESTS_COUNT 2
#define MATRIX_SIZE 400

#define TEST_ADD 1
#define TEST_MUL 1
#define TEST_LU 0

bool test_multiplication(int rank, int size, double &mpi_duration, double &normal_duration)
{
    std::clock_t start;
    bool test_result = false;
    MpiMatrixHelper helper(rank, size);
    mpi_duration = 0;
    normal_duration = 0;

    for(int i = 0; i < RANDOM_TESTS_COUNT; i++)
    {
      Generator gen(rank, size);
      auto test_matrix_1 = gen.GenerateRandomMatrix(MATRIX_SIZE, MATRIX_SIZE, 2*MATRIX_SIZE, column_wise);
      auto test_matrix_2 = gen.GenerateRandomMatrix(MATRIX_SIZE, MATRIX_SIZE, 2*MATRIX_SIZE, row_wise);
      
      start = std::clock();
      auto sparse_result = helper.mul(test_matrix_1, test_matrix_2);
      mpi_duration += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

      if (rank == 0)
      {
        dense_matrix actual(sparse_result);
        dense_matrix dense_m_1(test_matrix_1);
        dense_matrix dense_m_2(test_matrix_2);

        start = std::clock();
        auto expected = dense_m_1 * dense_m_2;
        normal_duration += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;  

        if(MATRIX_SIZE < 7)
        {
            printf("A:\n");
            dense_m_1.printDense();

            printf("B:\n");
            dense_m_2.printDense();

            printf("Result:\n");
            actual.printDense();

            printf("Expected:\n");
            expected.printDense();
        }

        test_result = (expected == actual);
      }

      MPI_Bcast(&test_result, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
      //printf("I am rank %d and test result is %d\n", rank, test_result);
      if(!test_result) return false;
    }

    mpi_duration /= RANDOM_TESTS_COUNT;
    normal_duration /= RANDOM_TESTS_COUNT;
    return true;
}

bool test_addition(int rank, int size, double &mpi_duration, double &normal_duration)
{
    std::clock_t start;
    bool test_result = false;
    MpiMatrixHelper helper(rank, size);
    mpi_duration = 0;
    normal_duration = 0;

    for(int i = 0; i < RANDOM_TESTS_COUNT; i++)
    {
      Generator gen(rank, size);
      auto test_matrix_1 = gen.GenerateRandomMatrix(MATRIX_SIZE, MATRIX_SIZE, 2*MATRIX_SIZE, column_wise);
      auto test_matrix_2 = gen.GenerateRandomMatrix(MATRIX_SIZE, MATRIX_SIZE, 2*MATRIX_SIZE, column_wise);
      start = std::clock();
      auto sparse_result = helper.add(test_matrix_1, test_matrix_2);
      mpi_duration += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      
      if (rank == 0)
      {
        dense_matrix actual(sparse_result);
        dense_matrix dense_m_1(test_matrix_1);
        dense_matrix dense_m_2(test_matrix_2);
        
        start = std::clock();
        auto expected = dense_m_1 + dense_m_2;
        normal_duration += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

        if(MATRIX_SIZE < 7)
        {
            printf("A:\n");
            dense_m_1.printDense();

            printf("B:\n");
            dense_m_2.printDense();

            printf("Result:\n");
            actual.printDense();

            printf("Expected:\n");
            expected.printDense();
        }

        test_result = (expected == actual);
      }

      MPI_Bcast(&test_result, 1, MPI_INT, 0, MPI_COMM_WORLD);
      //printf("I am rank %d and test result is %d\n", rank, test_result);
      if(!test_result) return false;
    }
    mpi_duration /= RANDOM_TESTS_COUNT;
    normal_duration /= RANDOM_TESTS_COUNT;
    return true;
}

int main(int argc, char** argv)
{
    int rank, size;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    double mpi_duration=0, normal_duration=0;

    if(rank == 0)
      srand(time(NULL));

    if(TEST_MUL)
    if(test_multiplication(rank, size, mpi_duration, normal_duration))
    {
      if(rank == 0)
        printf("test_multiplication [SUCCESS] | time: mpi=%f, normal=%f\n", mpi_duration, normal_duration);
    }
    else
    {
      if(rank == 0)
        printf("test_multiplication [FAIL]\n");
    }

    if(TEST_ADD)
    if(test_addition(rank, size, mpi_duration, normal_duration))
    {
        if(rank == 0)
            printf("test_addition [SUCCESS] | time mpi=%f, normal=%f\n", mpi_duration, normal_duration);
    }
    else
    {
        if(rank == 0)
            printf("test_addition [FAIL]\n");
    }

    MPI_Finalize();
    return 0;
}
