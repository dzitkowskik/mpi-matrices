#include <stdexcept>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include "mpimatrix.h"

using namespace std;

MpiMatrixHelper::MpiMatrixHelper()
{ init(); }

MpiMatrixHelper::MpiMatrixHelper(int rank, int proc_cnt)
        : rank(rank), processors_cnt(proc_cnt)
{ init(); }

MpiMatrixHelper::MpiMatrixHelper(const MpiMatrixHelper &m)
{
    rank = m.rank;
    processors_cnt = m.processors_cnt;
    sparse_elem_type = m.sparse_elem_type;
}

MpiMatrixHelper::~MpiMatrixHelper()
{ }

void MpiMatrixHelper::createSparseElemDatatype()
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

void MpiMatrixHelper::init()
{
    createSparseElemDatatype();
}

sparse_matrix MpiMatrixHelper::load(const char *path, MatrixType type, direction dir, int offset)
{
    if (rank != 0) return sparse_matrix();
    if (type == sparse) return sparse_matrix::fromSparseFile(path, dir, offset);
    else if (type == dense) return sparse_matrix::fromDenseFile(path, dir);
    else throw std::runtime_error("wrong data matrix file type");
}

void MpiMatrixHelper::sendVector(int node, sparse_vector vector)
{
    direction dir = vector.getDir();
    auto raw_data = vector.getElements(dir);
    int size = raw_data.size();

    MPI_Send(&size, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
    MPI_Send(&dir, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
    MPI_Send(&raw_data.front(), size, sparse_elem_type, node, 0, MPI_COMM_WORLD);
}

sparse_vector MpiMatrixHelper::receiveVector(int node)
{
    int size;
    direction dir;
    MPI_Recv(&size, 1, MPI_INT, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&dir, 1, MPI_INT, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    vector<sparse_matrix_elem> data(size);
    MPI_Recv(&data[0], size, sparse_elem_type, node, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return sparse_vector(size, dir, data);
}

void MpiMatrixHelper::sendMatrix(int node, sparse_matrix matrix)
{
    auto raw_data = matrix.getRawData();
    int size = static_cast<int>(raw_data.size());
    int width = matrix.getWidth();
    int height = matrix.getHeight();

    MPI_Send(&size, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
    MPI_Send(&raw_data.front(), size, sparse_elem_type, node, 0, MPI_COMM_WORLD);

    MPI_Send(&width, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
    MPI_Send(&height, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
}

sparse_matrix MpiMatrixHelper::receiveMatrix(int node, direction dir)
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
