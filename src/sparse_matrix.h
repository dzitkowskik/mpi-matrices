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
// FIELDS
private:
	vector<sparse_vector> data;
	direction dir;
	int width;
	int height;

// CONSTRUCTORS
public:
	sparse_matrix();
	sparse_matrix(int width, int height);
	sparse_matrix(int width, int height, direction d);
	sparse_matrix(vector<sparse_matrix_elem> elements, int width, int height, direction d);
	~sparse_matrix();
	sparse_matrix(const sparse_matrix &m);

// OPERATORS
public:
	sparse_matrix operator+(const sparse_matrix &m);
	sparse_matrix& operator+=(const sparse_matrix &m);
	sparse_matrix operator-(const sparse_matrix &m);
	sparse_matrix& operator-=(const sparse_matrix &m);
	sparse_matrix operator*(const sparse_matrix &m);
	sparse_vector & operator[](size_t el);
	const sparse_vector & operator[](size_t el) const;
	bool operator==(const sparse_matrix &m);
	bool operator!=(const sparse_matrix &m);

// METHODS
public:
	void resize(int w, int h);
	void transpose();
	void init();
	void clean();
	void printSparse();
	void createMatrixByRows(vector<sparse_matrix_elem> elements);
	void createMatrixByCols(vector<sparse_matrix_elem> elements);
	vector<sparse_matrix> splitToN(int N) const;
	static sparse_matrix fromSparseFile(const char *name, direction d);
	static sparse_matrix fromDenseFile(const char *name, direction d);
	vector<sparse_matrix_elem> getRawData() const;
	int getWidth() const;
	int getHeight() const;
	direction getDir() const;
	sparse_matrix getL();
	sparse_matrix getU();
	sparse_vector getRow(int n);
	sparse_vector getCol(int n);
};

#endif //__sparse_matrix_H_