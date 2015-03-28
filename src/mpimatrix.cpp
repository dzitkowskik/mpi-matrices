#include <cstdlib>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include "mpimatrix.h"

using namespace std;

void MpiMatrix::init()
{
    createSparseElemDatatype();
}

void MpiMatrix::loadFromFile(const char* path, direction dir)
{
    if (rank == 0)
        matrix = sparse_matrix::fromFile(path, dir);
}

void MpiMatrix::print()
{
    matrix.printSparse();
}

void MpiMatrix::createSparseElemDatatype()
{
    const int nitems=2;
    int          blocklengths[2] = {2,1};
    MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};
    MPI_Aint     offsets[2];
    MPI_Aint extent;
    MPI_Type_extent(MPI_INT, &extent);
    offsets[0] = 0;
    offsets[1] = 2 * extent;

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &sparse_elem_type);
    MPI_Type_commit(&sparse_elem_type);
}

void MpiMatrix::sendMatrix(int node, sparse_matrix matrix)
{
    int size = matrix.raw_data.size();
    int width = matrix.width;
    int height = matrix.height;

    MPI_Send(&size, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
    MPI_Send(&matrix.raw_data.front(), size, sparse_elem_type, node, 0, MPI_COMM_WORLD);

    MPI_Send(&width, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
    MPI_Send(&height, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
}

sparse_matrix MpiMatrix::receiveMatrix(int node, direction dir)
{
    int size;
    MPI_Recv(&size, 1, MPI_INT, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    vector<sparse_elem> data(size);
    MPI_Recv(&data[0], size, sparse_elem_type, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    int width, height;
    MPI_Recv(&width, 1, MPI_INT, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&height, 1, MPI_INT, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return sparse_matrix(data, width, height, dir);
}

MpiMatrix MpiMatrix::operator+(const MpiMatrix &m)
{
    int done;
    sparse_matrix result;

    if (rank == 0)
    {
         if (matrix.width != m.matrix.height)
            throw std::runtime_error("Dimensions of matrices do not match");

        if (matrix.width < processors_cnt || processors_cnt == 1)
        {
            // Do sequential multiplication
            result = sparse_matrix(
              sparse_matrix::addMatrices(matrix.raw_data, m.matrix.raw_data),
              matrix.width,
              matrix.height,
              row_wise);
            done = 1;
        }

        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            // Do parallel multiplication using MPI
            vector<sparse_matrix> matrices1 = matrix.splitToN(processors_cnt-1);
            vector<sparse_matrix> matrices2 = m.matrix.splitToN(processors_cnt-1);

            for(int i=1; i < processors_cnt; i++)
            {
                sendMatrix(i, matrices1[i-1]);
                sendMatrix(i, matrices2[i-1]);
            }

            vector<sparse_elem> elements;

            for(int i=1; i < processors_cnt; i++)
            {
                sparse_matrix part_result = receiveMatrix(i, column_wise);
                elements.insert(elements.begin(), part_result.raw_data.begin(), part_result.raw_data.end());
            }

            result = sparse_matrix(elements, matrix.width, matrix.height, column_wise);
        }
    }

    if (rank != 0)
    {
        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if(!done)
        {
            sparse_matrix matrix1 = receiveMatrix(0, column_wise);
            sparse_matrix matrix2 = receiveMatrix(0, column_wise);
            
            //printf("I am %d and received:\n", rank);
            //printf("matrix1 (w=%d, h=%d): \n", matrix1.width, matrix1.height);
            //matrix1.printSparse();
            //printf("matrix2 (w=%d, h=%d): \n", matrix2.width, matrix2.height);
            //matrix2.printSparse();


            result = sparse_matrix(sparse_matrix::addMatrices(matrix1.raw_data, matrix2.raw_data),
              0,
              0,
              column_wise);
            sendMatrix(0, result);
        }
    }

    return MpiMatrix(rank, processors_cnt, result);
}

MpiMatrix MpiMatrix::operator*(const MpiMatrix &m)
{
    int done;
    sparse_matrix result;

    if (rank == 0)
    {
        if (matrix.width != m.matrix.height)
            throw std::runtime_error("Dimensions of matrices do not match");

        if (matrix.width < processors_cnt || processors_cnt == 1)
        {
            // Do sequential multiplication
            result = sparse_matrix(
              sparse_matrix::multiplyMatrices(matrix, m.matrix),
              m.matrix.width,
              matrix.height,
              row_wise);
            done = 1;
        }

        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            // Do parallel multiplication using MPI
            vector<sparse_matrix> matrices_col = matrix.splitToN(processors_cnt-1);
            vector<sparse_matrix> matrices_row = m.matrix.splitToN(processors_cnt-1);

            for(int i=1; i < processors_cnt; i++)
            {
                sendMatrix(i, matrices_col[i-1]);
                sendMatrix(i, matrices_row[i-1]);
            }

            vector<sparse_matrix> result_matrices;
            vector<sparse_elem> acc;

            for(int i=1; i < processors_cnt; i++)
            {
                result_matrices.push_back(receiveMatrix(i, column_wise));
                acc = sparse_matrix::addMatrices(acc, result_matrices[i-1].raw_data);
            }

            result = sparse_matrix(acc, m.matrix.width, matrix.height, column_wise);
        }
    }

    if (rank != 0)
    {
        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if(!done)
        {
            sparse_matrix col_matrix = receiveMatrix(0, column_wise);
            sparse_matrix row_matrix = receiveMatrix(0, row_wise);

            //printf("I am %d and received:\n", rank);
            //printf("col_matrix (w=%d, h=%d): \n", col_matrix.width, col_matrix.height);
            //col_matrix.printSparse();
            //printf("row_matrix (w=%d, h=%d): \n", row_matrix.width, row_matrix.height);
            //row_matrix.printSparse();

            result = sparse_matrix(
              sparse_matrix::multiplyMatrices(col_matrix, row_matrix),
              0,
              0,
              row_wise);
            sendMatrix(0, result);
        }
    }

    return MpiMatrix(rank, processors_cnt, result);
}
