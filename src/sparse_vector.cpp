//
// Created by ghash on 29.03.15.
//

#include "sparse_vector.h"

// CONSTRUCTORS

sparse_vector::sparse_vector() : length(0)
{ }

sparse_vector::sparse_vector(int len) : length(len)
{ }

sparse_vector::sparse_vector(direction dir) : length(0), dir(dir)
{ }

sparse_vector::sparse_vector(int len, direction dir) : length(len), dir(dir)
{ }

sparse_vector::sparse_vector(int len, direction dir, std::vector<sparse_matrix_elem> elements)
		: length(len), dir(dir)
{
	for (auto it = elements.begin(); it != elements.end(); it++)
		set(dir == column_wise ? it->row : it->col, it->value);
}

sparse_vector::~sparse_vector()
{ }

sparse_vector::sparse_vector(const sparse_vector &other)
{
	data = other.data;
	length = other.length;
	dir = other.dir;
}

// GETTERS AND SETTERS

double sparse_vector::get(int nIndex) const
{
	if (nIndex >= length || nIndex < 0)
		throw std::runtime_error("index out of bounds");
	if (data.count(nIndex) > 0)
		return data.at(nIndex);
	return 0;
}

void sparse_vector::set(std::pair<int, double> item)
{
	if (item.first >= length) length = item.first + 1;
	if (item.second == 0) return;
	if (data.count(item.first) > 0)
		data[item.first] = item.second;
	data.insert(item);
}

void sparse_vector::set(int index, double value)
{
	this->set(std::make_pair(index, value));
}

void sparse_vector::reset(int len)
{
	length = len;
	//data.clear();
}

void sparse_vector::reset(int len, direction d)
{
	dir = d;
	reset(len);
}

int sparse_vector::size() const
{
	return length;
}

direction sparse_vector::getDir() const
{
	return dir;
}

void sparse_vector::setDir(direction d)
{
	dir = d;
}

void sparse_vector::clean()
{
	for (auto it = data.begin(); it != data.end();)
	{
		if (it->second == 0) data.erase(it++);
		else ++it;
	}
}

void sparse_vector::clear()
{
	data.clear();
}

// OPERATIONS

void sparse_vector::add(std::pair<int, double> item)
{
	if (item.first >= length)
		throw std::runtime_error("add/sub value to/from non existing item");
	if (item.second == 0) return;
	if (data.count(item.first) > 0)
		data[item.first] += item.second;
	else data.insert(item);
	clean();
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
	if (item.first >= length)
		throw std::runtime_error("mul a non existing item");
	if (data.count(item.first) > 0)
	{
		if (item.second == 0) data.erase(item.first);
		data[item.first] *= item.second;
	}
}

void sparse_vector::mul(int index, double value)
{
	this->mul(std::make_pair(index, value));
}

void sparse_vector::div(std::pair<int, double> item)
{
	if (item.second == 0)
		throw std::runtime_error("division by zero");
	if (item.first >= length)
		throw std::runtime_error("mul a non existing item");
	if (data.count(item.first) > 0)
	{
		if (item.second == 0) data.erase(item.first);
		data[item.first] /= item.second;
	}
}

void sparse_vector::div(int index, double value)
{
	this->div(std::make_pair(index, value));
}

// OTHER

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

void sparse_vector::print() const
{
	for (auto it = data.cbegin(); it != data.cend(); it++)
		printf("(%d)=>%f", it->first, it->second);
}

bool sparse_vector::operator==(const sparse_vector &v)
{
	try
	{
		for (int i = 0; i < length; i++)
			if (data[i] != v[i]) return false;
		return true;
	}
	catch(...)
	{
		return false;
	}
}

bool sparse_vector::operator!=(const sparse_vector &v)
{
	return !(*this == v);
}
