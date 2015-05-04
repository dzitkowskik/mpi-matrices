//
// Created by ghash on 04.05.15.
//

#include "cg.h"

#define CG_EPS 0.001
#define CG_MAX_ITERS 10

sparse_vector solveILU(const sparse_matrix &L, const sparse_matrix &U, const sparse_vector &b)
{
    int n = b.size();
    double *tmp = new double[n];
    sparse_vector x(n);
    auto sum = 0.0f;

    // Solve L
    for (int r = 0; r < n; ++r)
    {
        for (int c = 0; c < r; ++c)
            sum += L[c][r] * tmp[c];
        tmp[r] = (b[r] - sum) / L[r][r];
    }

    // Solve U
    for (int r = n - 1; r >= 0; --r)
    {
        for (int c = r + 1; c < n; ++c)
            sum += U[c][r] * x[c];
        x[r] = (tmp[r] - sum) / U[r][r];
    }

    delete[] tmp;
    return x;
}

int cg(const sparse_matrix &A, sparse_vector &x, sparse_vector &b, sparse_matrix &L, sparse_matrix &U)
{
    int iter = 1;
    double alpha, beta, rho0, rho1;
    sparse_vector p, z, q;

    double normb = b.l2_norm();
    sparse_vector r = b - (A * x);

    if (normb == 0.0) normb = 1.0;
    double residual = r.l2_norm() / normb;

    for(; iter <= CG_MAX_ITERS; iter++)
    {
        if(residual <= CG_EPS) break;

        z = solveILU(L, U, r);
        printf("r=\n");
        r.print();
        printf("z=\n");
        z.print();
        rho0 = r.dot(z);

        if(iter==1) p = z;
        else
        {
            beta = rho0 / rho1;
            p = z + p * beta;
        }

        q = A * p;
        alpha = rho0 / p.dot(q);
        x += p * alpha;
        r -= q * alpha;

        rho1 = rho0;
        residual = r.l2_norm() / normb;
    }

    return iter;
}