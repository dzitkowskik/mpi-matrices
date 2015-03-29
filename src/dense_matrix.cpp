#include "dense_matrix.h"
#include <stdio.h>
#include <assert.h>

dense_matrix dense_matrix::operator+(const dense_matrix &m)
{
    dense_matrix result(width, height);
    for(int i = 0; i < width; i++)
        for(int j = 0; j < height; j++)
            result.data[i][j] = data[i][j] + m.data[i][j];
    return result;
}

dense_matrix dense_matrix::operator*(const dense_matrix &m)
{
    dense_matrix result(m.width, height);
    for(int i = 0; i < m.width; i++)
        for(int j = 0; j < height; j++)
            for(int k = 0; k < width; k++)
                result.data[i][j] += data[k][j] * m.data[i][k];
    return result;
}

bool AlmostEqual2sComplement(float A, float B, int maxUlps = 4)
{
    // Make sure maxUlps is non-negative and small enough that the
    // default NAN won't compare as equal to anything.
    assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
    int aInt = *(int*)&A;
    // Make aInt lexicographically ordered as a twos-complement int
    if (aInt < 0)
        aInt = 0x80000000 - aInt;
    // Make bInt lexicographically ordered as a twos-complement int
    int bInt = *(int*)&B;
    if (bInt < 0)
        bInt = 0x80000000 - bInt;
    int intDiff = abs(aInt - bInt);
    if (intDiff <= maxUlps)
        return true;
    return false;
}

bool operator== (dense_matrix &m1, dense_matrix &m2)
{
    for(int i = 0; i < m1.width; i++)
      for(int j = 0; j < m1.height; j++)
        if(!AlmostEqual2sComplement(m1.data[i][j], m2.data[i][j]))
          return false;
    return true;
}

bool operator!= (dense_matrix &m1, dense_matrix &m2)
{
    return !(m1 == m2);
}

void dense_matrix::printDense()
{
    for(int j = 0; j < height; j++)
    {
        for(int i = 0; i < width; i++)
        {
          printf("\t%.2f", data[i][j]);
        }
        printf("\n");
    }
}
