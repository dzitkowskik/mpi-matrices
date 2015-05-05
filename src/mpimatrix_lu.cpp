//
// Created by ghash on 05.05.15.
//

#include "mpimatrix.h"

void MpiMatrixHelper::LU(const sparse_matrix &A, sparse_matrix &L, sparse_matrix &U)
{
    int done = 0;
    int width = A.getWidth();
    int height = A.getHeight();

    sparse_matrix local(A);
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
            auto matrices = A.splitToN(processors_cnt - 1);

            // Send submatrices to processors with their positions and with of original matrix
            int pos = 0;
            for (int i = 1; i < processors_cnt; i++)
            {
                sendMatrix(i, matrices[i - 1].first);
                MPI_Send(&pos, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&width, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                pos += matrices[i - 1].second;
            }

            // Receive result from the last processor and store it as local
            local = receiveMatrix(processors_cnt - 1, column_wise);
        }
    }
    if (rank > 0)
    {
        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (done) return;

        // Receive matrix and position
        int pos = 0, n = 0;
        local = receiveMatrix(0, column_wise);
        height = local.getHeight();
        MPI_Recv(&pos, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int min = pos;
        int len = n / (processors_cnt - 1);
        int max = rank == processors_cnt - 1 ? n - 1 : len - 1;

//        printf("I am proc nr %d and got n = %d and pos = %d\n", rank, n, pos);
//        printf("have min = %d, max = %d and matrix:\n", min, max);
//        local.printSparse();
//        printf("Width = %d, height = %d\n", width, height);
//        printf("\n");

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
        L = local.getL();
        U = local.getU();
    }
}

void MpiMatrixHelper::ILU(const sparse_matrix &A, sparse_matrix &L, sparse_matrix &U)
{
    int done = 0;
    int width = A.getWidth();
    int height = A.getHeight();

    sparse_matrix local(A);

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
                        if (local[i][j] != 0)
                            local[i][j] = local[i][j] - local[i][k] * local[k][j];
            }
            done = 1;
        }

        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!done)
        {
            // Split the matrix to submatrices
            auto matrices = A.splitToN(processors_cnt - 1);

            // Send submatrices to processors with their positions and with of original matrix
            int pos = 0;
            for (int i = 1; i < processors_cnt; i++)
            {
                sendMatrix(i, matrices[i - 1].first);
                MPI_Send(&pos, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&width, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                pos += matrices[i - 1].second;
            }

            // Receive result from the last processor and store it as local
            local = receiveMatrix(processors_cnt - 1, column_wise);
        }
    }
    if (rank > 0)
    {
        MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (done) return;

        // Receive matrix and position
        int pos = 0, n = 0;
        local = receiveMatrix(0, column_wise);
        height = local.getHeight();
        MPI_Recv(&pos, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int min = pos;
        int len = n / (processors_cnt - 1);
        int max = rank == processors_cnt - 1 ? n - 1 : len - 1;

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
                    if (local[i][j] != 0)
                        local[i][j] = local[i][j] - local[i][k] * local[k][j];
        }

        // Last processor sends result
        if (rank == processors_cnt - 1)
            sendMatrix(0, local);
    }
    if (rank == 0)
    {
        // Retrieve L and U matrices
        L = local.getL();
        U = local.getU();
    }
}