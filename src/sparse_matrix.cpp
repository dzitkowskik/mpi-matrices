//
// Created by Karol Dzitkowski on 08.03.2015.
// Copyright (c) 2015 Karol Dzitkowski. All rights reserved.
//

#include <fstream>
#include <sstream>
#include <map>
#include "sparse_matrix.h"
#include <stdexcept>

using namespace std;

sparse_matrix::sparse_matrix(vector<sparse_matrix_elem> elements, int width, int height, direction d)
		: dir(d), width(width), height(height)
{ fill(elements); }

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

vector<sparse_matrix_elem> readSparseElements(const char *name)
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

vector<sparse_matrix_elem> readDenseElements(const char *name, int &w, int &h)
{
	vector<sparse_matrix_elem> elements;
	string line;
	ifstream file(name);
	int i = 0, j = 0;
	while( getline(file, line) )
	{
		std::istringstream stream(line);
		double value;
		while(stream >> value)
		{
			elements.push_back(sparse_matrix_elem{j, i, value});
			j++;
		}
		i++;
		j = 0;
	}
	w = j;
	h = i;
	return elements;
}

sparse_matrix sparse_matrix::fromSparseFile(const char *name, direction d)
{
	return sparse_matrix(readSparseElements(name), 0, 0, d);
}

sparse_matrix sparse_matrix::fromDenseFile(const char *name, direction d)
{
	int width, height;
	auto elements = readDenseElements(name, width, height);
	return sparse_matrix(elements, width, height, d);
}

void sparse_matrix::printSparse() const
{
	auto raw_data = getRawData();
	printf("col\trow\tvalue\n");
	for (auto it = raw_data.cbegin(); it != raw_data.cend(); it++)
	{
		printf("%d\t%d\t%f\n", it->col, it->row, it->value);
	}
	printf("\n");
}

vector<pair<sparse_matrix, int>> sparse_matrix::splitToN(int N) const
{
	int size = dir == column_wise ? width : height;

	vector<pair<sparse_matrix, int>> result;
	vector<sparse_matrix_elem> elements;

	int len = size / N;
	int split_size = 0;
	int new_width = width, new_height = height;

	for (int i = 1, j = 0, n = 1; i <= size; i++)
	{
		if (j < len || n == N)
		{
			auto tmp = data[i - 1].getElements(dir, i - 1);
			elements.insert(elements.end(), tmp.begin(), tmp.end());
			j++;
			split_size++;
		}
		else
		{
			result.push_back(make_pair(sparse_matrix(elements, new_width, new_height, dir), split_size));
			split_size = 0;
			elements.clear();
			j=0;
			n++;
			i--;
		}
	}

	result.push_back(make_pair(sparse_matrix(elements, new_width, new_height, dir), split_size));
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

direction sparse_matrix::getDir() const
{
	return dir;
}

sparse_vector sparse_matrix::getRow(int n)
{
	if (dir == row_wise) return data[n];
	else
	{
		sparse_vector result(width, row_wise);
		for(int i=0; i<width; i++)
			result.set(i, data[i][n]);
		return result;
	}
}

sparse_vector sparse_matrix::getCol(int n)
{
	if (dir == column_wise) return data[n];
	else
	{
		sparse_vector result(height, column_wise);
		for(int i=0; i<height; i++)
			result.set(i, data[n][i]);
		return result;
	}
}

void sparse_matrix::fill(vector<sparse_matrix_elem> elements)
{
	init();
	if (dir == column_wise) createMatrixByCols(elements);
	else if (dir == row_wise) createMatrixByRows(elements);
}
