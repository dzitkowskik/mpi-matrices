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

	static sparse_matrix identity(int size, direction dir = column_wise);

// OPERATORS
public:
	sparse_matrix operator+(const sparse_matrix &m) const;
	sparse_matrix& operator+=(const sparse_matrix &m);
	sparse_matrix operator-(const sparse_matrix &m) const;
	sparse_matrix& operator-=(const sparse_matrix &m);
	sparse_matrix operator*(const sparse_matrix &m) const;
	sparse_vector operator*(const sparse_vector &v) const;

	sparse_vector & operator[](size_t el);
	const sparse_vector & operator[](size_t el) const;
	bool operator==(const sparse_matrix &m);
	bool operator!=(const sparse_matrix &m);

// METHODS
public:
	void fill(vector<sparse_matrix_elem> elements);
	void resize(int w, int h);
	void toggleDir();
	void transpose();
	void init();
	void clean();
	void printSparse() const;
	void printDense() const;
	void createMatrixByRows(vector<sparse_matrix_elem> elements);
	void createMatrixByCols(vector<sparse_matrix_elem> elements);
	vector<std::pair<sparse_matrix, int>> splitToN(int N) const;
	static sparse_matrix fromSparseFile(const char *name, direction d, int offset = 0);
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