//
// Created by Karol Dzitkowski on 05.05.15.
//

#include "mpimatrix.h"
#include <stdexcept>

sparse_matrix MpiMatrixHelper::add(const sparse_matrix &a, const sparse_matrix &b)
{
    int done = 0;
    sparse_matrix result(a.getWidth(), a.getHeight(), column_wise);

    if (rank == 0)
    {
        if (a.getWidth() != b.getHeight())
            throw std::runtime_error("Dimensions of matrices do not match");

        if (a.getWidth() < processors_cnt || processors_cnt == 1)
        {
            // Do sequential addition
            result = a + b;
            done = 1;
        }

        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            // Do parallel multiplication using MPI
            auto matrices1 = a.splitToN(processors_cnt - 1);
            auto matrices2 = b.splitToN(processors_cnt - 1);

            for (int i = 1; i < processors_cnt; i++)
            {
                sendMatrix(i, matrices1[i - 1].first);
                sendMatrix(i, matrices2[i - 1].first);
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

            result.fill(elements);
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

    return result;
}

void MpiMatrixHelper::addto(sparse_matrix &to, const sparse_matrix &what)
{
    int done = 0;
    if (rank == 0)
    {
        if (to.getWidth() != what.getWidth() || to.getHeight() != what.getHeight())
            throw std::runtime_error("Dimensions of matrices do not match");

        if (to.getWidth() < processors_cnt || processors_cnt == 1)
        {
            to += what;
            done = 1;
        }

        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            // Do parallel multiplication using MPI
            auto matrices1 = to.splitToN(processors_cnt - 1);
            auto matrices2 = what.splitToN(processors_cnt - 1);

            for (int i = 1; i < processors_cnt; i++)
            {
                sendMatrix(i, matrices1[i - 1].first);
                sendMatrix(i, matrices2[i - 1].first);
            }

            vector<sparse_matrix_elem> elements;
            for (int i = 1; i < processors_cnt; i++)
            {
                auto part_result = receiveMatrix(i, column_wise).getRawData();
                elements.insert(elements.begin(), part_result.begin(), part_result.end());
            }

            to.fill(elements);
        }
    }
    else
    {
        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            sparse_matrix matrix1 = receiveMatrix(0, column_wise);
            sparse_matrix matrix2 = receiveMatrix(0, column_wise);

            matrix1 += matrix2;
            sendMatrix(0, matrix1);
        }
    }
}

sparse_matrix MpiMatrixHelper::sub(const sparse_matrix &a, const sparse_matrix &b)
{
    int done = 0;
    sparse_matrix result(a.getWidth(), a.getHeight(), column_wise);

    if (rank == 0)
    {
        if (a.getWidth() != b.getHeight())
            throw std::runtime_error("Dimensions of matrices do not match");

        if (a.getWidth() < processors_cnt || processors_cnt == 1)
        {
            // Do sequential addition
            result = a - b;
            done = 1;
        }

        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            // Do parallel multiplication using MPI
            auto matrices1 = a.splitToN(processors_cnt - 1);
            auto matrices2 = b.splitToN(processors_cnt - 1);

            for (int i = 1; i < processors_cnt; i++)
            {
                sendMatrix(i, matrices1[i - 1].first);
                sendMatrix(i, matrices2[i - 1].first);
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

            result.fill(elements);
        }
    }

    if (rank != 0)
    {
        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            sparse_matrix matrix1 = receiveMatrix(0, column_wise);
            sparse_matrix matrix2 = receiveMatrix(0, column_wise);

            result = matrix1 - matrix2;
            sendMatrix(0, result);
        }
    }

    return result;
}

void MpiMatrixHelper::subto(sparse_matrix &to, const sparse_matrix &what)
{
    int done = 0;
    if (rank == 0)
    {
        if (to.getWidth() != what.getWidth() || to.getHeight() != what.getHeight())
            throw std::runtime_error("Dimensions of matrices do not match");

        if (to.getWidth() < processors_cnt || processors_cnt == 1)
        {
            to -= what;
            done = 1;
        }

        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            // Do parallel multiplication using MPI
            auto matrices1 = to.splitToN(processors_cnt - 1);
            auto matrices2 = what.splitToN(processors_cnt - 1);

            for (int i = 1; i < processors_cnt; i++)
            {
                sendMatrix(i, matrices1[i - 1].first);
                sendMatrix(i, matrices2[i - 1].first);
            }

            vector<sparse_matrix_elem> elements;
            for (int i = 1; i < processors_cnt; i++)
            {
                auto part_result = receiveMatrix(i, column_wise).getRawData();
                elements.insert(elements.begin(), part_result.begin(), part_result.end());
            }

            to.fill(elements);
        }
    }
    else
    {
        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            sparse_matrix matrix1 = receiveMatrix(0, column_wise);
            sparse_matrix matrix2 = receiveMatrix(0, column_wise);

//            printf("Received matrix 1 (%dx%d):\n", matrix1.getWidth(), matrix1.getHeight());
//            matrix1.printSparse();
//            printf("Received matrix 2 (%dx%d):\n", matrix2.getWidth(), matrix2.getHeight());
//            matrix2.printSparse();

            matrix1 -= matrix2;
            sendMatrix(0, matrix1);
        }
    }
}

sparse_matrix MpiMatrixHelper::mul(const sparse_matrix &a, const sparse_matrix &b)
{
    int done = 0;
    sparse_matrix result(b.getWidth(), a.getHeight(), column_wise);

    if (rank == 0)
    {
        if (a.getWidth() != b.getHeight())
            throw std::runtime_error("Dimensions of matrices do not match");

        if (a.getWidth() < processors_cnt || processors_cnt == 1)
        {
            // Do sequential multiplication
            result = a * b;
            done = 1;
        }

        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            // Do parallel multiplication using MPI
            auto matrices_col = a.splitToN(processors_cnt - 1);
            auto matrices_row = b.splitToN(processors_cnt - 1);

            for (int i = 1; i < processors_cnt; i++)
            {
                sendMatrix(i, matrices_col[i - 1].first);
                sendMatrix(i, matrices_row[i - 1].first);
            }

            for (int i = 1; i < processors_cnt; i++)
                result += receiveMatrix(i, column_wise);
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

    return result;
}

sparse_vector MpiMatrixHelper::mul(const sparse_matrix &A, const sparse_vector &x)
{
    int done = 0;
    sparse_vector result(x.size(), column_wise);

    if (rank == 0)
    {
        if (A.getWidth() != x.size())
            throw std::runtime_error("Dimensions of matrices do not match");

        if (A.getWidth() < processors_cnt || processors_cnt == 1)
        {
            // Do sequential multiplication
            result = A * x;
            done = 1;
        }

        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            // Do parallel multiplication using MPI
            auto matrices_col = A.splitToN(processors_cnt - 1);

            for (int i = 1; i < processors_cnt; i++)
            {
                sendVector(i, x);
                sendMatrix(i, matrices_col[i - 1].first);
            }

            for (int i = 1; i < processors_cnt; i++)
                result += receiveVector(i);
        }
    }

    if (rank != 0)
    {
        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            sparse_vector vec = receiveVector(0);
            sparse_matrix col_matrix = receiveMatrix(0, column_wise);

            result = col_matrix * vec;
            sendVector(0, result);
        }
    }

    return result;
}
