#include "generator.h"
#include <vector>
#include <random>
#include <iostream>
#include <chrono>
#include <algorithm>

using namespace std;



bool ExistElemKey(vector<sparse_matrix_elem> v, sparse_matrix_elem item)
{
    for(auto it = v.begin(); it != v.end(); it++)
      if(it->col == item.col && it->row == item.row)
        return true;
    return false;
}

MpiMatrix Generator::GenerateRandomMatrix(int width, int height, int num, direction dir)
{
  vector<sparse_matrix_elem> elements;

  for(int i = 0; i < num; i++)
  {
    int c = rand() % width;
    int r = rand() % height;
    double v = ((double) rand() / (RAND_MAX)) * 10;
    v = rand() % 10;

    if (ExistElemKey(elements, sparse_matrix_elem{ c, r, 0 }))
    {
      --i;
      continue;
    }
    elements.push_back(sparse_matrix_elem{ c, r, v });
  }

  sparse_matrix m(elements, width, height, dir);
  return MpiMatrix(rank, proc_cnt, m);
}
