//
// Created by Karol Dzitkowski on 30.03.15.
//

#include <assert.h>
#include "sparse_vector.h"
#include <math.h>
#include <stdexcept>

sparse_vector sparse_vector::operator+(const sparse_vector &m) const
{
	auto result = sparse_vector(m);
	for (auto it = data.cbegin(); it != data.cend(); it++)
		result.add(*it);
	result.clean();
	return result;
}

sparse_vector sparse_vector::operator-(const sparse_vector &m) const
{
	auto result = sparse_vector(*this);
	for (auto it = m.data.cbegin(); it != m.data.cend(); it++)
		result.sub(*it);
	result.clean();
	return result;
}

std::vector<sparse_matrix_elem> sparse_vector::operator*(const sparse_vector &m) const
{
	auto result = std::vector<sparse_matrix_elem>{};
	auto value = double{ 0.0f };
	if (dir == row_wise && m.dir == column_wise) // result will be one double
	{
		for (int i = 0; i < length; i++)
			value += get(i) * m.get(i);
		if (value != 0)
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
		for (auto it = data.cbegin(); it != data.cend(); it++)
		{
			for (auto it_m = m.cbegin(); it_m != m.cend(); it_m++)
			{
				value = it->second * it_m->second;
				if (value != 0)
					result.push_back(sparse_matrix_elem{it_m->first, it->first, value});
			}
		}
	}
	return result;
}

sparse_vector sparse_vector::operator*(const double &m) const
{
	auto result = sparse_vector(*this);
	if (m == 0) result.clear();
	for (int i = 0; i < result.size(); i++)
		result.mul(i, m);
	return result;
}

sparse_vector sparse_vector::operator/(const double &m) const
{
	auto result = sparse_vector(*this);
	if (m == 0)
		throw std::runtime_error("division by zero");
	for (int i = 0; i < result.size(); i++)
		result.div(i, m);
	return result;
}

sparse_vector sparse_vector::operator+(const double &m) const
{
	auto result = sparse_vector(*this);
	if (m == 0) return result;
	for (int i = 0; i < result.size(); i++)
		result.add(i, m);
	result.clean();
	return result;
}

sparse_vector sparse_vector::operator-(const double &m) const
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
	auto value = v;
	if (v == 0) data.clear();
	else
		for (int i = 0; i < size(); i++)
			mul(i, value);
	return *this;
}

sparse_vector &sparse_vector::operator+=(const double &v)
{
	auto value = v;
	if (v == 0) return *this;
	for (int i = 0; i < size(); i++)
		add(i, value);
	clean();
	return *this;
}

sparse_vector &sparse_vector::operator-=(const double &v)
{
	auto value = v;
	if (v == 0) return *this;
	for (int i = 0; i < size(); i++)
		sub(i, value);
	clean();
	return *this;
}

sparse_vector &sparse_vector::operator/=(const double &v)
{
	auto value = v;
	if (v == 0)
		throw std::runtime_error("division by zero");
	for (int i = 0; i < size(); i++)
		div(i, value);
	return *this;
}

sparse_vector& sparse_vector::operator+=(const sparse_vector &v)
{
	assert(length == v.size());
	for(int i=0; i < length; i++)
		set(i, get(i)+v[i]);
	this->clean();
	return *this;
}

sparse_vector& sparse_vector::operator-=(const sparse_vector &v)
{
	assert(length == v.size());
	for(int i=0; i < length; i++)
		set(i, get(i)-v[i]);
	this->clean();
	return *this;
}

bool sparse_vector::operator==(const sparse_vector &v) const
{
	try
	{
		for (int i = 0; i < length; i++)
		{
			double diff = (*this)[i] - v[i];
			if (abs(diff) > 0.01) return false;
		}
		return true;
	}
	catch(...)
	{
		return false;
	}
}

bool sparse_vector::operator!=(const sparse_vector &v) const
{
	return !(*this == v);
}


