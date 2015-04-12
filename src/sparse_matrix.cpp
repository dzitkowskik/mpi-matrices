//
// Created by Karol Dzitkowski on 08.03.2015.
// Copyright (c) 2015 Karol Dzitkowski. All rights reserved.
//

#include <fstream>
#include <sstream>
#include <map>
#include "sparse_matrix.h"

using namespace std;

sparse_matrix::sparse_matrix(vector<sparse_matrix_elem> elements, int width, int height, direction d)
		: dir(d), width(width), height(height)
{
	init();
	if (dir == column_wise) createMatrixByCols(elements);
	else if (dir == row_wise) createMatrixByRows(elements);
}

sparse_matrix::sparse_matrix(const sparse_matrix &m)
{
	dir = m.dir;
	data = m.data;
	width = m.width;
	height = m.height;
}

sparse_matrix::sparse_matrix(int width, int height, direction d) : dir(d), width(width), height(height)
{ init(); }

sparse_matrix::sparse_matrix(int width, int height) : width(width), height(height)
{ }

sparse_matrix::sparse_matrix() : width(0), height(0)
{ }

sparse_matrix::~sparse_matrix()
{ }

void sparse_matrix::init()
{
	int w_size = dir == column_wise ? width : height;
	int h_size = dir == column_wise ? height : width;
	data.clear();
	data.resize(w_size);
	for (int i = 0; i < w_size; i++) data[i].reset(h_size, dir);
}

void sparse_matrix::resize(int w, int h)
{
	width = w;
	height = h;
	int w_size = dir == column_wise ? width : height;
	int h_size = dir == column_wise ? height : width;
	data.resize(w_size);
	for (int i = 0; i < w_size; i++)
		data[i].reset(h_size, dir);
}

void sparse_matrix::transpose()
{
	auto raw_data = getRawData();
	std::swap(width, height);
	if (dir == column_wise)
	{
		dir = row_wise;
		init();
		createMatrixByCols(raw_data);
	}
	else
	{
		dir = column_wise;
		init();
		createMatrixByRows(raw_data);
	}
}

void sparse_matrix::createMatrixByRows(vector<sparse_matrix_elem> elements)
{
	for (auto it = elements.begin(); it != elements.end(); it++)
	{
		auto elem = *it;
		if (height < elem.row + 1) resize(width, elem.row + 1);
		if (width < elem.col + 1) resize(elem.col + 1, height);
		data[elem.row].set(elem.col, elem.value);
	}
}

void sparse_matrix::createMatrixByCols(vector<sparse_matrix_elem> elements)
{
	for (auto it = elements.begin(); it != elements.end(); it++)
	{
		auto elem = *it;
		if (width < elem.col + 1) resize(elem.col + 1, height);
		if (height < elem.row + 1) resize(width, elem.row + 1);
		data[elem.col].set(elem.row, elem.value);
	}
}

vector<sparse_matrix_elem> readElements(const char *name)
{
	vector<sparse_matrix_elem> elements;
	string line;
	ifstream infile(name);

	while (getline(infile, line) && !line.empty())
	{
		istringstream iss(line);
		int col, row;
		double val;
		if (!(iss >> col >> row >> val))
		{
			throw std::runtime_error(string("Error while reading file: ") + name);
		}
		elements.push_back(sparse_matrix_elem{col, row, val});
	}
	return elements;
}

sparse_matrix sparse_matrix::fromFile(const char *name, direction d)
{
	return sparse_matrix(readElements(name), 0, 0, d);
}

void sparse_matrix::printSparse()
{
	auto raw_data = getRawData();
	printf("\ncol\trow\tvalue\n");
	for (auto it = raw_data.cbegin(); it != raw_data.cend(); it++)
	{
		printf("%d\t%d\t%f\n", it->col, it->row, it->value);
	}
}

vector<sparse_matrix> sparse_matrix::splitToN(int N) const
{
	int size = dir == column_wise ? width : height;

	vector<sparse_matrix> result;
	vector<sparse_matrix_elem> elements;

	double step = (double) size / N;
	double edge = step;

	int new_width = 0, new_height = 0;

	for (int i = 1; i <= size; i++)
	{
		if (i <= edge)
		{
			if (dir == column_wise)
			{
				new_width++;
				new_height = height;
			}
			else
			{
				new_height++;
				new_width = width;
			}
			auto tmp = data[i - 1].getElements(dir, i - 1);
			elements.insert(elements.end(), tmp.begin(), tmp.end());
		}
		else
		{
			result.push_back(sparse_matrix(elements, new_width, new_height, dir));
			elements.clear();
			edge += step;
			i--;
		}
	}

	result.push_back(sparse_matrix(elements, new_width, new_height, dir));
	return result;
}

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
	catch(exception e)
	{
		printf("sparse_matrix operator+ : %s", e.what());
	}
	return sparse_matrix(elements, width, height, dir);
}

sparse_matrix sparse_matrix::operator*(const sparse_matrix &m)
{
	sparse_matrix result(m.width, height, dir);

	for (int i = 0; i < width; i++)
	{
		sparse_matrix tmp = sparse_matrix(data[i] * m.data[i], m.width, height, dir);
		result = result + tmp;
	}

	return result;
}

vector<sparse_matrix_elem> sparse_matrix::getRawData() const
{
	std::vector<sparse_matrix_elem> elements;
	try
	{
		for(int i=0; i<data.size(); i++)
		{
			auto tmp = data[i].getElements(dir, i);
			if (tmp.size() > 0)
				elements.insert(elements.end(), tmp.begin(), tmp.end());
		}
	}
	catch(exception e)
	{
		printf("getRawData: %s", e.what());
	}
	return elements;
}

int sparse_matrix::getWidth() const
{ return width; }

int sparse_matrix::getHeight() const
{ return height; }

sparse_vector &sparse_matrix::operator[](size_t el)
{ return data[el]; }

const sparse_vector &sparse_matrix::operator[](size_t el) const
{ return data[el]; }

sparse_matrix sparse_matrix::getL()
{
	sparse_matrix result(width, height, dir);
	for(int i=0; i<width; i++)
		for(int j=i; j<height; j++)
		{
			if (i==j) result[i][i] = 1.0;
			else result[i][j] = (*this)[i][j];
		}
	result.clean();
	return result;
}

sparse_matrix sparse_matrix::getU()
{
	sparse_matrix result(width, height, dir);
	for(int i=0; i<width; i++)
		for(int j=0; j<=i; j++)
			result[i][j] = (*this)[i][j];
	result.clean();
	return result;
}

void sparse_matrix::clean()
{
	for(int i=0; i<width; i++)
		data[i].clean();
}
