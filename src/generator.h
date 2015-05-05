#include "mpimatrix.h"


class Generator
{
    int rank;
    int proc_cnt;

public:
    Generator(int rank, int proc_cnt) : rank(rank), proc_cnt(proc_cnt) {};
    ~Generator() {};

    sparse_matrix GenerateRandomMatrix(int width, int height, int num, direction dir);
};
