//
// Created by Karol Dzitkowski on 08.03.2015.
// Copyright (c) 2015 Karol Dzitkowski. All rights reserved.
//

#include <fstream>
#include <sstream>
#include <map>
#include <stdexcept>
#include "sparse_matrix.h"

using namespace std;

void sparse_matrix::createMatrixByRows(vector<sparse_elem> elements)
{
	data.resize(height);

	for(auto it= elements.begin(); it != elements.end(); it++)
	{
		auto elem = *it;
		if(height < elem.row)
    {
      data.resize(elem.row);
      height = elem.row;
    }
		data[elem.row-1].push_back(make_tuple(elem.col, elem.value));
	}
}

void sparse_matrix::createMatrixByCols(vector<sparse_elem> elements)
{
	data.resize(width);

	for(auto it = elements.begin(); it != elements.end(); it++)
	{
		auto elem = *it;
		if(width < elem.col)
    {
      data.resize(elem.col);
      width = elem.col;
    }
		data[elem.col-1].push_back(make_tuple(elem.row, elem.value));
	}
}

vector<sparse_elem> readElements(const char *name)
{
	vector<sparse_elem> elements;
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
		elements.push_back(sparse_elem{col, row, val});
	}
	return elements;
}

sparse_matrix sparse_matrix::fromFile(const char *name, direction d)
{
	return sparse_matrix(readElements(name), 0, 0, d);
}

void addToMap(map<tuple<int,int>,double>* value_map, vector<sparse_elem>* v)
{
	for (auto it = v->begin(); it != v->end(); it++)
	{
		auto key = make_tuple(it->col, it->row);
		if(value_map->count(key)>0)
			(*value_map)[key]+=it->value;
		else value_map->insert(pair<tuple<int,int>,double>(key, it->value));
	}
}

vector<sparse_elem> sparse_matrix::addMatrices(vector<sparse_elem> a, vector<sparse_elem> b)
{
	vector<sparse_elem> v;
	map<tuple<int,int>,double> value_map;
	addToMap(&value_map, &a);
	addToMap(&value_map, &b);
	for(map<tuple<int,int>,double>::iterator it = value_map.begin(); it != value_map.end(); it++)
	{
		v.push_back(sparse_elem{get<0>(it->first), get<1>(it->first), it->second});
	}
	return v;
}

vector<sparse_elem> sparse_matrix::multiplyVectors(vector<tuple<int,double>> col, vector<tuple<int,double>> row)
{
	vector<sparse_elem> v;

	for (auto it_col = col.begin(); it_col != col.end(); it_col++)
	{
		for (auto it_row = row.begin(); it_row != row.end(); it_row++)
		{
			double value = get<1>(*it_row) * get<1>(*it_col);
			v.push_back(sparse_elem{get<0>(*it_row), get<0>(*it_col), value});
		}
	}
	return v;
}

vector<sparse_elem> sparse_matrix::multiplyMatrices(sparse_matrix m_column, sparse_matrix m_row)
{
	vector<sparse_elem> v;

	for(int i = 0; i < m_column.data.size(); i++)
	{
		auto x = sparse_matrix::multiplyVectors(m_column.data[i], m_row.data[i]);
		v = sparse_matrix::addMatrices(v, x);
	}

	return v;
}

void sparse_matrix::printSparse()
{
	printf("\ncol\trow\tvalue\n");
	for(auto it = raw_data.begin(); it != raw_data.end(); it++)
	{
		printf("%d\t%d\t%f\n", it->col, it->row, it->value);
	}
}

vector<sparse_matrix> sparse_matrix::splitToN(int N) const
{
  int size = dir == column_wise ? width : height;

  vector<sparse_matrix> result;
  vector<sparse_elem> elements;

  double step = (double)size/N;
  double edge = step;

  int new_width=0, new_height=0;

	for(int i = 1; i <= size; i++)
  {
    if(i <= edge)
    {
      if(dir == column_wise)
      {
        new_width++;
        new_height = height;
      }
      else
      {
        new_height++;
        new_width = width;
      }

			for(auto it=data[i-1].begin(); it != data[i-1].end(); it++)
			{
				if(dir == column_wise)
					elements.push_back(sparse_elem{i, get<0>(*it), get<1>(*it)});
				else
					elements.push_back(sparse_elem{get<0>(*it), i, get<1>(*it)});
			}
    }
    else
    {
      //printf("width %d height %d\n", width, height);
      //printf("new_width %d new_height %d\n", new_width, new_height);

      result.push_back(sparse_matrix(elements, new_width, new_height, dir));
      elements.resize(0);
      edge += step;
      i--;
    }
  }

  //printf("width %d height %d\n", width, height);
  //printf("new_width %d new_height %d\n", new_width, new_height);

	result.push_back(sparse_matrix(elements, new_width, new_height, dir));
  return result;
}
