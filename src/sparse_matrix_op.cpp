//
// Created by Karol Dzitkowski on 04.05.15.
//

#include "sparse_matrix.h"

void addToMap(map<pair<int, int>, double> *value_map, const vector<sparse_matrix_elem> &v)
{
    for (auto it = v.cbegin(); it != v.cend(); it++)
    {
        auto key = make_pair(it->col, it->row);
        if (value_map->count(key) > 0)
            (*value_map)[key] += it->value;
        else value_map->insert(pair<pair<int, int>, double>(key, it->value));
    }
}

void subToMap(map<pair<int, int>, double> *value_map, const vector<sparse_matrix_elem> &v)
{
    for (auto it = v.cbegin(); it != v.cend(); it++)
    {
        auto key = make_pair(it->col, it->row);
        if (value_map->count(key) > 0)
            (*value_map)[key] -= it->value;
        else value_map->insert(pair<pair<int, int>, double>(key, it->value));
    }
}

sparse_matrix sparse_matrix::operator+(const sparse_matrix &m)
{
    vector<sparse_matrix_elem> elements;
    try
    {
        map<pair<int, int>, double> value_map;
        addToMap(&value_map, getRawData());
        addToMap(&value_map, m.getRawData());
        for (map<pair<int, int>, double>::iterator it = value_map.begin(); it != value_map.end(); it++)
        {
            elements.push_back(sparse_matrix_elem{get<0>(it->first), get<1>(it->first), it->second});
        }
    }
    catch (exception e)
    {
        printf("sparse_matrix operator+ : %s", e.what());
    }
    return sparse_matrix(elements, width, height, dir);
}

sparse_matrix sparse_matrix::operator-(const sparse_matrix &m)
{
    vector<sparse_matrix_elem> elements;
    try
    {
        map<pair<int, int>, double> value_map;
        subToMap(&value_map, getRawData());
        subToMap(&value_map, m.getRawData());
        for (map<pair<int, int>, double>::iterator it = value_map.begin(); it != value_map.end(); it++)
        {
            elements.push_back(sparse_matrix_elem{get<0>(it->first), get<1>(it->first), it->second});
        }
    }
    catch (exception e)
    {
        printf("sparse_matrix operator+ : %s", e.what());
    }
    return sparse_matrix(elements, width, height, dir);
}

sparse_matrix sparse_matrix::operator*(const sparse_matrix &m)
{
    sparse_matrix result(m.width, height, dir);

    // TODO: FIX THIS - IT IS VERY SLOW!!!

    for (int i = 0; i < width; i++)
    {
        sparse_matrix tmp = sparse_matrix(data[i] * m.data[i], m.width, height, dir);
        result = result + tmp;
    }

    return result;
}

sparse_matrix &sparse_matrix::operator+=(const sparse_matrix &m)
{
    for (int i = 0; i < width; i++)
        for (int j = 0; j < height; j++)
            if (data[i][j] != 0) data[i][j] += m[i][j];
    this->clean();
    return *this;
}

sparse_matrix &sparse_matrix::operator-=(const sparse_matrix &m)
{
    for (int i = 0; i < width; i++)
        for (int j = 0; j < height; j++)
            if (data[i][j] != 0) data[i][j] -= m[i][j];
    this->clean();
    return *this;
}

sparse_vector sparse_matrix::operator*(const sparse_vector &v) const
{
    auto n = v.size();
    sparse_vector result(n);

    for(int i = 0; i < n; i++)
        result[i] = data[i].sum() * v[i];

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

