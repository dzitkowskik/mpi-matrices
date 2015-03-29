//
// Created by Karol Dzitkowski on 29.03.15.
//

#ifndef MPI_MATRICES_SPARSE_VECTOR_H
#define MPI_MATRICES_SPARSE_VECTOR_H

#include <map>
#include <vector>
#include "direction.h"
#include "sparse_matrix_elem.h"

class sparse_vector
{
private:
	std::map<int, double> data;
	int length;
	direction dir;

public:
	sparse_vector() : length(0)
	{ }

	sparse_vector(int len) : length(len)
	{ }

	sparse_vector(direction dir) : length(0), dir(dir)
	{ }

	sparse_vector(int len, direction dir) : length(len), dir(dir)
	{ }

	sparse_vector(const sparse_vector &other)
	{
		data = other.data;
		length = other.length;
		dir = other.dir;
	}

	~sparse_vector()
	{ }

	void setDir(direction d);

	sparse_vector operator+(const sparse_vector &m);

	sparse_vector operator+(const double &m);

	sparse_vector operator-(const sparse_vector &m);

	sparse_vector operator-(const double &m);

	std::vector<sparse_matrix_elem> operator*(const sparse_vector &m);

	sparse_vector operator*(const double &m);

	double get(int nIndex) const;

	void set(std::pair<int, double> item);

	void set(int index, double value);

	void add(std::pair<int, double> item);

	void add(int index, double value);

	void sub(std::pair<int, double> item);

	void sub(int index, double value);

	void mul(std::pair<int, double> item);

	void mul(int index, double value);

	void reset(int len);

	void reset(int len, direction d);

	int size() const;

	direction getDir() const;

	bool isColumn() const;

	bool isRow() const;

	std::vector<sparse_matrix_elem> getElements(direction d, int x = 0) const;
};


#endif //MPI_MATRICES_SPARSE_VECTOR_H
