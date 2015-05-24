//
// Created by Karol Dzitkowski on 05.05.15.
//

#include "mpimatrix.h"

#define CG_EPS 0.001
#define CG_MAX_ITERS 1000

sparse_vector solveILU(const sparse_matrix &L, const sparse_matrix &U, const sparse_vector &b)
{
    int n = b.size();
    double *tmp = new double[n];
    sparse_vector x(n);


    // Solve L
    for (int r = 0; r < n; ++r)
    {
        auto sum = 0.0f;
        for (int c = 0; c < r; ++c)
            sum += L[c][r] * tmp[c];
        tmp[r] = (b[r] - sum) / L[r][r];
    }

    // Solve U
    for (int r = n - 1; r >= 0; --r)
    {
        auto sum = 0.0f;
        for (int c = r + 1; c < n; ++c)
            sum += U[c][r] * x[c];
        x[r] = (tmp[r] - sum) / U[r][r];
    }

    delete[] tmp;
    return x;
}

sparse_vector MpiMatrixHelper::CG(const sparse_matrix &A, const sparse_vector &b)
{
    sparse_vector x(b.size(), column_wise);
    sparse_matrix L, U;
    LU(A, L, U);

    double alpha, beta, rho0, rho1, norm_b, residual;
    sparse_vector p, z, q, r;

    if (rank == 0)
    {
        norm_b = b.l2_norm();
        r = b; // - mul(A, x); // but x is zero vector in this case
        if (norm_b == 0.0) norm_b = 1.0;
        residual = r.l2_norm() / norm_b;
    }

    int iters = 0;
    for(int iteration = 1; iteration <= CG_MAX_ITERS; iteration++)
    {
        iters++;
        MPI_Bcast(&residual, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if(residual <= CG_EPS) break;

        if(rank == 0)
        {
            z = solveILU(L, U, r);
            rho0 = r.dot(z);

            if (iteration == 1) p = z;
            else
            {
                beta = rho0 / rho1;
                p = z + p * beta;
            }
        }

        q = mul(A, p);

        if(rank == 0)
        {
            alpha = rho0 / p.dot(q);
            x += p * alpha;
            r -= q * alpha;

            rho1 = rho0;
            residual = r.l2_norm() / norm_b;
        }
    }

    if(rank == 0)
        printf("\n\n ITER CNT = %d", iters);

    return x;
}