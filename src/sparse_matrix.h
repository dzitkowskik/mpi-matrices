//
// Created by Karol Dzitkowski on 08.03.2015.
// Copyright (c) 2015 Karol Dzitkowski. All rights reserved.
//


#ifndef __sparse_matrix_H_
#define __sparse_matrix_H_

#include <cstdlib>
#include <tuple>
#include <vector>
#include <array>
#include "sparse_matrix_elem.h"
#include "direction.h"
#include "sparse_vector.h"

using namespace std;

class sparse_matrix
{
private:
	vector<sparse_matrix_elem> raw_data;
	vector<sparse_vector> data;
	direction dir;
	int width;
	int height;

public:
	sparse_matrix() : width(0), height(0)
	{ }

	sparse_matrix(int width, int height) : width(width), height(height)
	{ }

	sparse_matrix(int width, int height, direction d) : dir(d), width(width), height(height)
	{ init(); }

	sparse_matrix(vector<sparse_matrix_elem> elements, int width, int height, direction d)
			: raw_data(elements), dir(d), width(width), height(height)
	{
		init();
		if (dir == column_wise) createMatrixByCols(elements);
		else if (dir == row_wise) createMatrixByRows(elements);
	}

	~sparse_matrix()
	{ }

	sparse_matrix(const sparse_matrix &m)
	{
		dir = m.dir;
		data = m.data;
		raw_data = m.raw_data;
		width = m.width;
		height = m.height;
	}

	void resize(int w, int h);

	void transpose();

	void init();

	void printSparse();

	void createMatrixByRows(vector<sparse_matrix_elem> elements);

	void createMatrixByCols(vector<sparse_matrix_elem> elements);

	vector<sparse_matrix> splitToN(int N) const;

	static sparse_matrix fromFile(const char *name, direction d);

	sparse_matrix operator+(const sparse_matrix &m);

	sparse_matrix operator*(const sparse_matrix &m);

	const vector<sparse_matrix_elem> &getRawData() const
	{ return raw_data; }

	int numberOfElements()
	{ return raw_data.size(); }

	int getWidth() const
	{ return width; }

	int getHeight() const
	{ return height; }
};


#endif //__sparse_matrix_H_
