//
// Created by Karol Dzitkowski on 30.03.15.
//

#include "sparse_vector.h"

sparse_vector sparse_vector::operator+(const sparse_vector &m)
{
	auto result = sparse_vector(m);
	for (auto it = data.begin(); it != data.end(); it++)
		result.add(*it);
	result.clean();
	return result;
}

sparse_vector sparse_vector::operator-(const sparse_vector &m)
{
	auto result = sparse_vector(m);
	for (auto it = data.begin(); it != data.end(); it++)
		result.sub(*it);
	result.clean();
	return result;
}

std::vector<sparse_matrix_elem> sparse_vector::operator*(const sparse_vector &m)
{
	auto result = std::vector<sparse_matrix_elem>{};
	auto value = double{ 0.0f };
	if (dir == row_wise && m.dir == column_wise) // result will be one double
	{
		for (int i = 0; i < length; i++)
			value += this->get(i) * m.get(i);
		result.push_back(sparse_matrix_elem{0, 0, value});
	}
	else if (dir == m.dir) // result will be vector
	{
		for (int i = 0; i < length; i++)
		{
			value = get(i) * m.get(i);
			if (value == 0) continue;
			if (dir == column_wise)
				result.push_back(sparse_matrix_elem{0, i, value});
			else result.push_back(sparse_matrix_elem{i, 0, value});
		}
	}
	else // result will be matrix
	{
		for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < m.length; j++)
			{
				value = get(i) * m.get(j);
				if (value != 0)
					result.push_back(sparse_matrix_elem{j, i, value});
			}
		}
	}
	return result;
}

sparse_vector sparse_vector::operator*(const double &m)
{
	auto result = sparse_vector(*this);
	if (m == 0) result.clear();
	for (int i = 0; i < result.size(); i++)
		result.mul(i, m);
	return result;
}

sparse_vector sparse_vector::operator+(const double &m)
{
	auto result = sparse_vector(*this);
	if (m == 0) return result;
	for (int i = 0; i < result.size(); i++)
		result.add(i, m);
	result.clean();
	return result;
}

sparse_vector sparse_vector::operator-(const double &m)
{
	auto result = sparse_vector(*this);
	if (m == 0) return result;
	for (int i = 0; i < result.size(); i++)
		result.sub(i, m);
	result.clean();
	return result;
}

double& sparse_vector::operator[](int nIndex)
{
	if (nIndex >= length || nIndex < 0)
		throw std::runtime_error("index out of bounds");
	if (data.count(nIndex) <= 0)
		data.insert(std::make_pair(nIndex,0));
	return data[nIndex];
}

const double sparse_vector::operator[](int el) const
{
	return get(el);
}

sparse_vector& sparse_vector::operator*=(const double &v)
{
	if (v == 0) data.clear();
	for (int i = 0; i < size(); i++)
		mul(i, v);
	return *this;
}

sparse_vector &sparse_vector::operator+=(const double &v)
{
	if (v == 0) return *this;
	for (int i = 0; i < size(); i++)
		add(i, v);
	clean();
	return *this;
}

sparse_vector &sparse_vector::operator-=(const double &v)
{
	if (v == 0) return *this;
	for (int i = 0; i < size(); i++)
		sub(i, v);
	clean();
	return *this;
}
