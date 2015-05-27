//
// Created by Karol Dzitkowski on 05.05.15.
//

#include "mpimatrix.h"

#define CG_EPS 1e-6
#define CG_MAX_ITERS 500

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

    // CONJUGATE GRADIENT
    double alpha, beta, rho0, rho1, norm_b, residual;
    sparse_vector p, q, r;
    int iteration;

    if (rank == 0)
    {
        norm_b = b.l2_norm();
        r = p = b;
        if (norm_b == 0.0) norm_b = 1.0;
        residual = r.l2_norm() / norm_b;
        rho1 = r.dot(r);
    }

    for(iteration = 1; iteration <= CG_MAX_ITERS; iteration++)
    {
        MPI_Bcast(&residual, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if(residual <= CG_EPS) break;

        q = mul(A, p);
        if(rank == 0)
        {
            rho0 = r.dot(r);
            alpha = rho0 / p.dot(q);
            x += p * alpha;
            r -= q * alpha;
            beta = rho0 / rho1;
            p = r + p * beta;
            residual = r.l2_norm() / norm_b;
            rho1 = rho0;
        }

        //printf("Residual: %f\n", residual);
    }

    if(rank == 0) printf("ITER CNT = %d\n", iteration-1);
    return x;
}

sparse_vector MpiMatrixHelper:: PRECONDITIONED_2(const sparse_matrix &A, const sparse_vector &b, const sparse_matrix &M_inv)
{
    sparse_vector x(b.size(), column_wise);
    for(int i=0; i<b.size(); i++)
        x[i] = 1;

    double alpha, beta, rho0, rho1, norm_b, residual;
    sparse_vector p, z, q, r;
    int iteration;

    if (rank == 0)
    {
        norm_b = b.l2_norm();
        r = b - mul(A, x); // but x is zero vector in this case
        if (norm_b == 0.0) norm_b = 1.0;
        residual = r.l2_norm() / norm_b;
    }

    p = mul(M_inv, r);
    z = p; // z0

    for(iteration = 1; iteration <= CG_MAX_ITERS; iteration++)
    {
        q = mul(A, p);
        if(rank == 0)
        {
            rho0 = r.dot(z);
            alpha = rho0 / p.dot(q);
            x += p * alpha;
            r -= q * alpha;
            residual = r.l2_norm() / norm_b;
        }
        z = mul(M_inv, r);
        if(rank == 0)
        {
            rho1 = r.dot(z);
            beta = rho1 / rho0;
            p = z + p * beta;
        }
        MPI_Bcast(&residual, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if(residual <= CG_EPS) break;
        printf("Residual: %f\n", residual);
    }

    if(rank == 0)
        printf("ITER CNT = %d\n", iteration);

    return x;
}

sparse_vector MpiMatrixHelper::CG_ILU(const sparse_matrix &A, const sparse_vector &b)
{
    sparse_vector x(b.size(), column_wise);
    sparse_matrix L, U;
    ILU(A, L, U);

    double alpha, beta, rho0, rho1, norm_b, residual;
    sparse_vector p, z, q, r;
    int iteration;

    if (rank == 0)
    {
        norm_b = b.l2_norm();
        r = b; // - mul(A, x); // but x is zero vector in this case
        if (norm_b == 0.0) norm_b = 1.0;
        residual = r.l2_norm() / norm_b;
    }

    for(iteration = 1; iteration <= CG_MAX_ITERS; iteration++)
    {
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
        MPI_Bcast(&residual, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if(residual <= CG_EPS) break;
        printf("IT: %d, Residual = %f\n", iteration, residual);
    }

    if(rank == 0)
        printf("ITER CNT = %d\n", iteration);

    return x;
}

sparse_vector MpiMatrixHelper::CG_ILU_PRECONDITIONED(const sparse_matrix &A, const sparse_vector &b)
{
    sparse_matrix L, U;
    ILU(A, L, U);
    auto L_inv = Inverse(L);
    auto L_inv_T(L_inv);
    if(rank == 0) L_inv_T.transpose();
    auto left = mul(L_inv, A);
    auto A_prec = mul(left, L_inv_T);
    auto b_prec = mul(L_inv, b);

    A_prec.clean();
    b_prec.clean();

    auto x = CG(A_prec, b_prec);
    return mul(L_inv_T, x);
}

sparse_vector MpiMatrixHelper::CG_ILU_PRECONDITIONED_WITH_ILU(const sparse_matrix &A, const sparse_vector &b){
    sparse_matrix L, U;
    ILU(A, L, U);
    auto L_inv = Inverse(L);
    auto L_inv_T(L_inv);
    if(rank == 0) L_inv_T.transpose();
    auto left = mul(L_inv, A);
    auto A_prec = mul(left, L_inv_T);
    auto b_prec = mul(L_inv, b);

    A_prec.clean();
    b_prec.clean();

    auto x = CG_ILU(A_prec, b_prec);
    return mul(L_inv_T, x);
}

sparse_vector MpiMatrixHelper::CG_ILU_PRECONDITIONED_COPY(const sparse_matrix &A, const sparse_vector &b)
{
    sparse_matrix L, U;
    ILU(A, L, U);
    auto L_inv = Inverse(L);
    auto L_inv_T(L_inv);
    if(rank == 0) L_inv_T.transpose();
    auto A_prec = mul(mul(L_inv, A), L_inv_T);
    auto b_prec = mul(L_inv, b);
    auto x = CG(A_prec, b_prec);
    return mul(L_inv_T, x);
}

sparse_vector MpiMatrixHelper::CG_ILU_PRECONDITIONED_2(const sparse_matrix &A, const sparse_vector &b)
{
    sparse_matrix L, U;
    ILU(A, L, U);

    auto L_T(L);
    L_T.transpose();
    auto M = mul(L, L_T);
    auto M_inv = Inverse(M);

    sparse_vector x(b.size(), column_wise);
    for(int i=0; i<b.size(); i++)
        x[i] = 1;

    double alpha, beta, rho0, rho1, norm_b, residual;
    sparse_vector p, z, q, r;
    int iteration;

    if (rank == 0)
    {
        norm_b = b.l2_norm();
        r = b - mul(A, x); // but x is zero vector in this case
        if (norm_b == 0.0) norm_b = 1.0;
        residual = r.l2_norm() / norm_b;
    }

    p = mul(M_inv, r);
    z = p; // z0

    for(iteration = 1; iteration <= CG_MAX_ITERS; iteration++)
    {
        q = mul(A, p);
        if(rank == 0)
        {
            rho0 = r.dot(z);
            alpha = rho0 / p.dot(q);
            x += p * alpha;
            r -= q * alpha;
            residual = r.l2_norm() / norm_b;
        }
        z = mul(M_inv, r);
        if(rank == 0)
        {
            rho1 = r.dot(z);
            beta = rho1 / rho0;
            p = z + p * beta;
        }
        MPI_Bcast(&residual, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if(residual <= CG_EPS) break;
    }

    if(rank == 0)
        printf("ITER CNT = %d\n", iteration);

    return x;
}