//
// Created by Karol Dzitkowski on 26.05.15.
//

#include "mpimatrix.h"


// Solves A * x = b where A is lower triangular matrix
sparse_vector MpiMatrixHelper::solveTrian(const sparse_matrix &A, const sparse_vector &b)
{
    int n = b.size();
    sparse_vector x(n);

    // Solve L
    for (int r = 0; r < n; ++r)
    {
        auto sum = 0.0f;
        for (int c = 0; c < r; ++c)
            sum += A[c][r] * x[c];
        x[r] = (b[r] - sum) / A[r][r];
    }

    return x;
}

// Solves multiple linear equations: A * X = B where A is lower triangular
// matrix and X and B are matrices in parallel
sparse_matrix MpiMatrixHelper::SolveManyTrian(const sparse_matrix &A, const sparse_matrix &B)
{
    int done = 0;
    sparse_matrix result(B.getWidth(), B.getHeight(), column_wise);

    if(rank == 0)
    {
        if (B.getDir() != column_wise) throw new std::runtime_error("B must be column matrix");
        int n = B.getWidth();
        if(n < processors_cnt || processors_cnt == 1)
        {
            // Solve sequentially
            for(int i=0; i<B.getWidth(); i++)
            {
                auto solution = solveTrian(A, B[i]);
                result[i] = solution;
            }

            done = 1;
        }

        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            int pos = 0;
            auto b_parts = B.splitToN(processors_cnt - 1);
            for (int i = 1; i < processors_cnt; i++)
            {
                sendMatrix(i, A);
                sendMatrix(i, b_parts[i - 1].first);
                MPI_Send(&pos, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                pos += b_parts[i - 1].second;
                MPI_Send(&pos, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }

            for (int i = 1; i < processors_cnt; i++)
                result += receiveMatrix(i, column_wise);
        }

    }
    if(rank > 0)
    {
        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            auto left = receiveMatrix(0, column_wise);
            auto right = receiveMatrix(0, column_wise);
            int from, to;
            MPI_Recv(&from, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&to, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            sparse_matrix part_result(right.getWidth(), right.getHeight(), column_wise);

//            printf("I am %d from = %d to = %d\n", rank, from, to);

            for(int i=from; i<to; i++)
            {
                auto solution = solveTrian(left, right[i]);
                part_result[i] = solution;
            }

            sendMatrix(0, part_result);
        }
    }

    return result;
}

sparse_matrix MpiMatrixHelper::Inverse(const sparse_matrix &A)
{
    // A * A^(-1) = I
    auto I = sparse_matrix::identity(A.getWidth());
    return SolveManyTrian(A, I);
}
