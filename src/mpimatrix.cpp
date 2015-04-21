#include <cstdlib>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include "mpimatrix.h"

using namespace std;

MpiMatrix::MpiMatrix()
{ init(); }

MpiMatrix::MpiMatrix(int rank, int proc_cnt)
        : rank(rank), processors_cnt(proc_cnt)
{ init(); }

MpiMatrix::MpiMatrix(int rank, int proc_cnt, sparse_matrix sp)
        : rank(rank), processors_cnt(proc_cnt), matrix(sp)
{ init(); }

MpiMatrix::MpiMatrix(const MpiMatrix &m)
{
    matrix = m.matrix;
    rank = m.rank;
    processors_cnt = m.processors_cnt;
    sparse_elem_type = m.sparse_elem_type;
}

MpiMatrix::~MpiMatrix()
{ }

void MpiMatrix::createSparseElemDatatype()
{
    const int nitems = 2;
    int blocklengths[2] = {2, 1};
    MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};
    MPI_Aint offsets[2];
    MPI_Aint extent;
    MPI_Type_extent(MPI_INT, &extent);
    offsets[0] = 0;
    offsets[1] = 2 * extent;

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &sparse_elem_type);
    MPI_Type_commit(&sparse_elem_type);
}

void MpiMatrix::init()
{
    createSparseElemDatatype();
}

MpiMatrix MpiMatrix::load(const char *path, int rank, int proc_cnt, MatrixType type)
{
    MpiMatrix result(rank, proc_cnt);
    if (type == sparse)
        result.loadSparse(path, column_wise);
    else
        result.loadDense(path, column_wise);
    return result;
}

void MpiMatrix::loadDense(const char *path, direction dir)
{
    if (rank == 0)
        matrix = sparse_matrix::fromDenseFile(path, dir);
}

void MpiMatrix::loadSparse(const char *path, direction dir)
{
    if (rank == 0)
        matrix = sparse_matrix::fromSparseFile(path, dir);
}

void MpiMatrix::print()
{
    matrix.printSparse();
}

void MpiMatrix::sendVector(int node, sparse_vector vector)
{
    direction dir = vector.getDir();
    auto raw_data = vector.getElements(dir);
    int size = raw_data.size();

    MPI_Send(&size, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
    MPI_Send(&dir, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
    MPI_Send(&raw_data.front(), size, sparse_elem_type, node, 0, MPI_COMM_WORLD);
}

sparse_vector MpiMatrix::receiveVector(int node)
{
    int size;
    direction dir;
    MPI_Recv(&size, 1, MPI_INT, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&dir, 1, MPI_INT, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    vector<sparse_matrix_elem> data(size);
    MPI_Recv(&data[0], size, sparse_elem_type, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return sparse_vector(size, dir, data);
}

void MpiMatrix::sendMatrix(int node, sparse_matrix matrix)
{
    auto raw_data = matrix.getRawData();
    int size = raw_data.size();
    int width = matrix.getWidth();
    int height = matrix.getHeight();

    MPI_Send(&size, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
    MPI_Send(&raw_data.front(), size, sparse_elem_type, node, 0, MPI_COMM_WORLD);

    MPI_Send(&width, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
    MPI_Send(&height, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
}

sparse_matrix MpiMatrix::receiveMatrix(int node, direction dir)
{
    int size;
    MPI_Recv(&size, 1, MPI_INT, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    vector<sparse_matrix_elem> data(size);
    MPI_Recv(&data[0], size, sparse_elem_type, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    int width, height;
    MPI_Recv(&width, 1, MPI_INT, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&height, 1, MPI_INT, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return sparse_matrix(data, width, height, dir);
}

MpiMatrix MpiMatrix::operator+(const MpiMatrix &m)
{
    int done = 0;
    sparse_matrix result;

    if (rank == 0)
    {
        if (matrix.getWidth() != m.matrix.getHeight())
            throw std::runtime_error("Dimensions of matrices do not match");

        if (matrix.getWidth() < processors_cnt || processors_cnt == 1)
        {
            // Do sequential addition
            result = matrix + m.matrix;
            done = 1;
        }

        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            // Do parallel multiplication using MPI
            vector<sparse_matrix> matrices1 = matrix.splitToN(processors_cnt - 1);
            vector<sparse_matrix> matrices2 = m.matrix.splitToN(processors_cnt - 1);

            for (int i = 1; i < processors_cnt; i++)
            {
                sendMatrix(i, matrices1[i - 1]);
                sendMatrix(i, matrices2[i - 1]);
            }

            vector<sparse_matrix_elem> elements;

            for (int i = 1; i < processors_cnt; i++)
            {
                auto part_result = receiveMatrix(i, column_wise).getRawData();
                elements.insert(
                        elements.begin(),
                        part_result.begin(),
                        part_result.end());
            }

            result = sparse_matrix(elements, matrix.getWidth(), matrix.getHeight(), column_wise);
        }
    }

    if (rank != 0)
    {
        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            sparse_matrix matrix1 = receiveMatrix(0, column_wise);
            sparse_matrix matrix2 = receiveMatrix(0, column_wise);

//			printf("I am %d and received:\n", rank);
//			printf("matrix1 (w=%d, h=%d): \n", matrix1.getWidth(), matrix1.getHeight());
//			matrix1.printSparse();
//			printf("matrix2 (w=%d, h=%d): \n", matrix2.getWidth(), matrix2.getHeight());
//			matrix2.printSparse();

            result = matrix1 + matrix2;
            sendMatrix(0, result);
        }
    }

    return MpiMatrix(rank, processors_cnt, result);
}

MpiMatrix MpiMatrix::operator*(const MpiMatrix &m)
{
    int done = 0;
    sparse_matrix result;

    if (rank == 0)
    {
        if (matrix.getWidth() != m.matrix.getHeight())
            throw std::runtime_error("Dimensions of matrices do not match");

        if (matrix.getWidth() < processors_cnt || processors_cnt == 1)
        {
            // Do sequential multiplication
            result = matrix * m.matrix;
            done = 1;
        }

        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            // Do parallel multiplication using MPI
            vector<sparse_matrix> matrices_col = matrix.splitToN(processors_cnt - 1);
            vector<sparse_matrix> matrices_row = m.matrix.splitToN(processors_cnt - 1);

            for (int i = 1; i < processors_cnt; i++)
            {
                sendMatrix(i, matrices_col[i - 1]);
                sendMatrix(i, matrices_row[i - 1]);
            }

            result = sparse_matrix(m.matrix.getWidth(), matrix.getHeight(), column_wise);
            for (int i = 1; i < processors_cnt; i++)
                result = result + receiveMatrix(i, column_wise);
        }
    }

    if (rank != 0)
    {
        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            sparse_matrix col_matrix = receiveMatrix(0, column_wise);
            sparse_matrix row_matrix = receiveMatrix(0, row_wise);


//			printf("I am %d and received:\n", rank);
//			printf("col_matrix (w=%d, h=%d): \n", col_matrix.getWidth(), col_matrix.getHeight());
//			col_matrix.printSparse();
//			printf("row_matrix (w=%d, h=%d): \n", row_matrix.getWidth(), row_matrix.getHeight());
//			row_matrix.printSparse();

            result = col_matrix * row_matrix;
            sendMatrix(0, result);
        }
    }

    return MpiMatrix(rank, processors_cnt, result);
}

void MpiMatrix::LU(MpiMatrix &L, MpiMatrix &U)
{
    int done = 0;
    int width = matrix.getWidth();
    int height = matrix.getHeight();

    sparse_matrix local(matrix);
    // We want to have column wise sparse matrix
    if (local.getDir() == row_wise) local.transpose();

    if (rank == 0)
    {
        if (width < processors_cnt || processors_cnt == 1)
        {
            // Do sequential LU
            for (int k = 0; k < width; k++)
            {
                for (int i = k + 1; i < height; i++)
                    local[k][i] /= local[k][k];
                for (int i = k + 1; i < width; i++)
                    for (int j = k + 1; j < height; j++)
                        local[i][j] = local[i][j] - local[i][k] * local[k][j];
            }
            done = 1;
        }

        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            // Split the matrix to submatrices
            vector<sparse_matrix> matrices = matrix.splitToN(processors_cnt - 1);

            // Send submatrices to processors with their positions and with of original matrix
            int pos = 0;
            for (int i = 1; i < processors_cnt; i++)
            {
                sendMatrix(i, matrices[i - 1]);
                MPI_Send(&pos, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&width, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                pos += matrices[i - 1].getWidth();
            }

            // Receive result from the last processor and store it as local
            local = receiveMatrix(processors_cnt - 1, column_wise);
        }
    }
    if (rank > 0)
    {
        // Receive matrix and position
        int pos = 0, n = 0;
        local = receiveMatrix(0, column_wise);
        height = local.getHeight();
        MPI_Recv(&pos, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int min = pos;
        int len = n / (processors_cnt - 1);
        int max = rank == processors_cnt - 1 ? n - 1 : len - 1;

//		printf("I am proc nr %d and got n = %d and pos = %d\n", rank, n, pos);
//		printf("have min = %d, max = %d and matrix:\n", min, max);
//		local.printSparse();
//		printf("Width = %d, height = %d\n", width, height);
//		printf("\n");

        // Compute
        for (int k = 0; k < n; k++)
        {
            if (k <= max)
            {
                if (k >= min)
                {
                    for (int i = k + 1; i < height; i++)
                        local[k][i] /= local[k][k];
                    for (int i = rank + 1; i < processors_cnt; i++)
                        sendVector(i, local[k]);
                }
                else local[k] = receiveVector(MPI_ANY_SOURCE);
            }
            for (int i = (((k + 1) > min) ? (k + 1) : min); i <= max; i++)
                for (int j = k + 1; j < n; j++)
                    local[i][j] = local[i][j] - local[i][k] * local[k][j];
        }

        // Last processor sends result
        if (rank == processors_cnt - 1)
            sendMatrix(0, local);
    }
    if (rank == 0)
    {
        // Retrieve L and U matrices
        L = MpiMatrix(rank, processors_cnt, local.getL());
        U = MpiMatrix(rank, processors_cnt, local.getU());
    }
}

bool MpiMatrix::operator==(MpiMatrix const &m)
{
    return matrix == m.matrix;
}
