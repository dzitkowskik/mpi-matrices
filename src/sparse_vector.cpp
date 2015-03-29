//
// Created by ghash on 29.03.15.
//

#include "sparse_vector.h"

double sparse_vector::get(int nIndex) const
{
	if (nIndex >= length || nIndex < 0)
		throw std::runtime_error("index out of bounds");
	if (data.count(nIndex) > 0)
		return data.at(nIndex);
	return 0;
}

sparse_vector sparse_vector::operator+(const sparse_vector &m)
{
	auto result = sparse_vector(m);
	for (auto it = data.begin(); it != data.end(); it++)
		result.add(*it);
}

sparse_vector sparse_vector::operator-(const sparse_vector &m)
{
	auto result = sparse_vector(m);
	for (auto it = data.begin(); it != data.end(); it++)
		result.add(*it);
}

std::vector<sparse_matrix_elem> sparse_vector::operator*(const sparse_vector &m)
{
	auto result = std::vector<sparse_matrix_elem>{};
	auto value = 0.0f;
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
	if (m == 0) data.clear();
	auto result = sparse_vector(*this);
	for (int i = 0; i < result.size(); i++)
		result.mul(i, m);
	return result;
}

sparse_vector sparse_vector::operator+(const double &m)
{
	auto result = sparse_vector(*this);
	if (m == 0) return result;
	for (int i = 0; i < result.size(); i++)
	{
		result.add(i, m);
		if (get(i) == 0) data.erase(i);
	}
	return result;
}

sparse_vector sparse_vector::operator-(const double &m)
{
	auto result = sparse_vector(*this);
	if (m == 0) return result;
	for (int i = 0; i < result.size(); i++)
	{
		result.sub(i, m);
		if (get(i) == 0) data.erase(i);
	}
	return result;
}

void sparse_vector::set(std::pair<int, double> item)
{
	if (item.first >= length) length = item.first + 1;
	if (item.second == 0) return;
	if (data.count(item.first) > 0)
		data[item.first] = item.second;
	data.insert(item);
}

void sparse_vector::reset(int len)
{
	length = len;
	data.clear();
}

void sparse_vector::reset(int len, direction d)
{
	dir = d;
	reset(len);
}

void sparse_vector::set(int index, double value)
{
	this->set(std::make_pair(index, value));
}

int sparse_vector::size() const
{
	return length;
}

bool sparse_vector::isColumn() const
{
	return dir == column_wise;
}

bool sparse_vector::isRow() const
{
	return dir == row_wise;
}

direction sparse_vector::getDir() const
{
	return dir;
}

void sparse_vector::add(std::pair<int, double> item)
{
	if (item.first >= length)
		throw std::runtime_error("add value to non existing item");
	if (item.second == 0) return;
	if (data.count(item.first) > 0)
		data[item.first] += item.second;
	else data.insert(item);
}

void sparse_vector::add(int index, double value)
{
	this->add(std::make_pair(index, value));
}

void sparse_vector::sub(std::pair<int, double> item)
{
	item.second = -item.second;
	this->add(item);
}

void sparse_vector::sub(int index, double value)
{
	this->sub(std::make_pair(index, value));
}

void sparse_vector::mul(std::pair<int, double> item)
{
	if (item.second == 0) data.erase(item.first);
	if (data.count(item.first) > 0)
		data[item.first] *= item.second;
}

void sparse_vector::mul(int index, double value)
{
	this->mul(std::make_pair(index, value));
}

std::vector<sparse_matrix_elem> sparse_vector::getElements(direction d, int x) const
{
	std::vector<sparse_matrix_elem> result;
	if (d == column_wise)
		for (auto it = data.cbegin(); it != data.cend(); it++)
			result.push_back(sparse_matrix_elem{x, it->first, it->second});
	else
		for (auto it = data.cbegin(); it != data.cend(); it++)
			result.push_back(sparse_matrix_elem{it->first, x, it->second});
	return result;
}

void sparse_vector::setDir(direction d)
{
	dir = d;
}
