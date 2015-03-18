#include "generator.h"
#include <vector>
#include <random>
#include <iostream>
#include <chrono>
#include <algorithm>

using namespace std;



bool ExistElemKey(vector<sparse_elem> v, sparse_elem item)
{
    for(auto it = v.begin(); it != v.end(); it++)
      if(it->col == item.col && it->row == item.row)
        return true;
    return false;
}

MpiMatrix Generator::GenerateRandomMatrix(int width, int height, int num, direction dir)
{
  vector<sparse_elem> elements;

  for(int i = 0; i < num; i++)
  {
    int c = rand() % width + 1;
    int r = rand() % height + 1;
    double v = ((double) rand() / (RAND_MAX));
    v = rand() % 9 + 1;
    if (ExistElemKey(elements, sparse_elem{ c, r, 0 }))
    {
      --i;
      continue;
    }
    elements.push_back(sparse_elem{ c, r, v });
  }

  sparse_matrix m(elements, width, height, dir);
  return MpiMatrix(rank, proc_cnt, m);
}
