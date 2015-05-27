//
// Created by Karol Dzitkowski on 04.05.15.
//

#include <stdio.h>
#include "sparse_matrix.h"

sparse_matrix sparse_matrix::operator+(const sparse_matrix &m) const
{
    if(dir != m.getDir()) throw new std::runtime_error("adding matrices with different directions");
    int w = width > m.getWidth() ? width : m.getWidth();
    int h = height > m.getHeight() ? height : m.getHeight();
    sparse_matrix result(w, h, dir);
    for(int i=0; i<data.size(); i++)
        result[i] = data[i] + m[i];
    return result;
}

sparse_matrix sparse_matrix::operator-(const sparse_matrix &m) const
{
    if(dir != m.getDir()) throw new std::runtime_error("adding matrices with different directions");
    int w = width > m.getWidth() ? width : m.getWidth();
    int h = height > m.getHeight() ? height : m.getHeight();
    sparse_matrix result(w, h, dir);
    for(int i=0; i<data.size(); i++)
        result[i] = data[i] - m[i];
    return result;
}

sparse_matrix sparse_matrix::operator*(const sparse_matrix &m) const
{
    sparse_matrix result(m.width, height, dir);

    for (int i = 0; i < data.size(); i++)
    {
        vector<sparse_matrix_elem> part_result = data[i] * m.data[i];
        for(vector<sparse_matrix_elem>::iterator it = part_result.begin(); it != part_result.end(); it++)
        {
            result[it->col][it->row] += it->value;
        }
    }
//    result.clean();
    return result;
}

sparse_matrix &sparse_matrix::operator+=(const sparse_matrix &m)
{
    if(dir != m.getDir()) throw new std::runtime_error("adding matrices with different directions");
    for (int i = 0; i < data.size(); i++)
            data[i] += m[i];
    return *this;
}

sparse_matrix &sparse_matrix::operator-=(const sparse_matrix &m)
{
    for (int i = 0; i < width; i++)
        for (int j = 0; j < height; j++)
            data[i][j] -= m[i][j];
    this->clean();
    return *this;
}

sparse_vector sparse_matrix::operator*(const sparse_vector &v) const
{
    auto n = v.size();
    sparse_vector result(n, v.getDir());

    if (dir == row_wise)
    {
        for (int i = 0; i < n; i++)
            result[i] = data[i].dot(v);
    } else {
        auto local = *this;
        local.toggleDir();
        for (int i = 0; i < n; i++)
            result[i] = local[i].dot(v);
    }

    return result;
}

sparse_vector &sparse_matrix::operator[](size_t el)
{ return data[el]; }

const sparse_vector &sparse_matrix::operator[](size_t el) const
{ return data[el]; }

bool sparse_matrix::operator==(const sparse_matrix &m)
{
    try
    {
        for (int i = 0; i < data.size(); i++)
            if (data[i] != m[i]) return false;
        return true;
    } catch (...)
    {
        return false;
    }
}

bool sparse_matrix::operator!=(const sparse_matrix &m)
{
    return !(*this == m);
}
