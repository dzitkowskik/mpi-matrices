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

using namespace std;

enum direction { column_wise, row_wise };

struct sparse_elem
{
	int col;
	int row;
	double value;
};

class sparse_matrix
{
public:
	vector<sparse_elem> raw_data;
	vector<vector<tuple<int,double>>> data;
	direction dir;
  int width;
  int height;

  sparse_matrix() : width(0), height(0) {}
	sparse_matrix(int width, int height) : width(width), height(height) {}
	sparse_matrix(vector<sparse_elem> elements, int width, int height, direction d)
   : raw_data(elements), dir(d), width(width), height(height)
	{
		if(dir==column_wise) createMatrixByCols(elements);
		else if(dir==row_wise) createMatrixByRows(elements);
	}
	~sparse_matrix(){}
	sparse_matrix(const sparse_matrix& m)
	{
		dir = m.dir;
		data = m.data;
    raw_data = m.raw_data;
    width = m.width;
    height = m.height;
	}
  void printSparse();
	void createMatrixByRows(vector<sparse_elem> elements);
	void createMatrixByCols(vector<sparse_elem> elements);
  vector<sparse_matrix> splitToN(int N) const;

	static sparse_matrix fromFile(const char* name, direction d);
	static vector<sparse_elem> multiplyVectors(vector<tuple<int,double>> col, vector<tuple<int,double>> row);
	static vector<sparse_elem> addMatrices(vector<sparse_elem> a, vector<sparse_elem> b);
	static vector<sparse_elem> multiplyMatrices(sparse_matrix m_column, sparse_matrix m_row);
};


#endif //__sparse_matrix_H_
