//
// Created by Karol Dzitkowski on 08.03.2015.
// Copyright (c) 2015 Karol Dzitkowski. All rights reserved.
//

#ifndef __dense_matrix_H_
#define __dense_matrix_H_

#include <cstdlib>
#include <tuple>
#include <vector>
#include <array>
#include "sparse_matrix.h"

using namespace std;

class dense_matrix
{
public:
  double** data;
  int height;
  int width;

  dense_matrix(){}
  dense_matrix(int width, int height) : width(width), height(height)
  {
    initEmptyData();
  }
  dense_matrix(vector<sparse_elem> elements)
  {
    initSize(elements);
    initEmptyData();
    initData(elements);
  }
  dense_matrix(const sparse_matrix& matrix)
  {
    width = matrix.width;
    height = matrix.height;
    initEmptyData();
    initData(matrix.raw_data);
  }
  ~dense_matrix(){}
  dense_matrix(const dense_matrix& m)
  {
    data = m.data;
    height = m.height;
    width = m.width;
  }

  void initSize(vector<sparse_elem> elements)
  {
    height=0;
    width=0;
    for (auto it=elements.begin(); it != elements.end(); it++)
    {
      if(it->row > height) height = it->row;
      if(it->col > width) width = it->col;
    }
  }

  void initData(vector<sparse_elem> elements)
  {
    for (auto it=elements.begin(); it != elements.end(); it++)
      data[(it->col)-1][(it->row)-1] = it->value;
  }

  void initEmptyData()
  {
    data = new double*[width];
    for(int i = 0; i < width; i++)
      data[i] = new double[height];

    for(int i = 0; i < width; i++)
      for(int j = 0; j < height; j++)
        data[i][j] = 0;
  }

  void printDense();

  dense_matrix operator+(const dense_matrix &m);
  dense_matrix operator*(const dense_matrix &m);

  friend bool operator== (dense_matrix &m1, dense_matrix &m2);
  friend bool operator!= (dense_matrix &m1, dense_matrix &m2);
};


#endif //__dense_matrix_H_
