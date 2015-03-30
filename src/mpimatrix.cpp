#include <cstdlib>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include "mpimatrix.h"

using namespace std;

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

void MpiMatrix::loadFromFile(const char *path, direction dir)
{
	if (rank == 0)
		matrix = sparse_matrix::fromFile(path, dir);
}

void MpiMatrix::print()
{
	matrix.printSparse();
}

void MpiMatrix::sendMatrix(int node, sparse_matrix matrix)
{
	int size = matrix.numberOfElements();
	int width = matrix.getWidth();
	int height = matrix.getHeight();

	MPI_Send(&size, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
	MPI_Send(&matrix.getRawData().front(), size, sparse_elem_type, node, 0, MPI_COMM_WORLD);

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
				sparse_matrix part_result = receiveMatrix(i, column_wise);
				elements.insert(
						elements.begin(),
						part_result.getRawData().begin(),
						part_result.getRawData().end());
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
	int size = matrix.getWidth();
	if (rank == 0)
	{
		if (matrix.getWidth() < processors_cnt || processors_cnt == 1)
		{
			// Do sequential LU
			for(int k=0; k<size; k++)
			{
				//matrix.divCol(k, matrix.get(k,k));
				for(int i=k+1; i < size; i++)
				{
					for(int j=k+1; j < matrix.getHeight(); j++)
					{

					}
				}
			}
			done = 1;
		}

		// Send submatrices to processors
	}
	if (rank > 0)
	{

	}
	if (rank == 0)
	{

	}
}
