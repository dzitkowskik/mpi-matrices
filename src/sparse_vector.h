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
// FIELDS
private:
	std::map<int, double> data;
	int length;
	direction dir;

// CONSTRUCTORS
public:
	sparse_vector();
	sparse_vector(int len);
	sparse_vector(direction dir);
	sparse_vector(int len, direction dir);
	sparse_vector(int len, direction dir, std::vector<sparse_matrix_elem> elements);
	sparse_vector(const sparse_vector &other);
	~sparse_vector();

private:
	sparse_vector(double);
	sparse_vector(float);

// OPERATORS
public:
	double &operator[](int el);
	const double operator[](int el) const;

	sparse_vector operator+(const sparse_vector &v) const;
	sparse_vector operator-(const sparse_vector &v) const;
	std::vector<sparse_matrix_elem> operator*(const sparse_vector &v) const;

	sparse_vector operator+(const double &d) const;
	sparse_vector operator-(const double &d) const;
	sparse_vector operator*(const double &d) const;
	sparse_vector operator/(const double &d) const;

	sparse_vector& operator+=(const double &d);
	sparse_vector& operator-=(const double &d);
	sparse_vector& operator*=(const double &d);
	sparse_vector& operator/=(const double &d);

	sparse_vector& operator+=(const sparse_vector &v);
	sparse_vector& operator-=(const sparse_vector &v);

	bool operator==(const sparse_vector &v) const;
	bool operator!=(const sparse_vector &v) const;

// GETTERS AND SETTERS
public:
	double get(int nIndex) const;
	direction getDir() const;

	void set(std::pair<int, double> item);
	void set(int index, double value);
	void setDir(direction d);

// METHODS
public:
	// OPERATIONS
	void add(std::pair<int, double> item);
	void add(int index, double value);
	void sub(std::pair<int, double> item);
	void sub(int index, double value);
	void mul(std::pair<int, double> item);
	void mul(int index, double value);
	void div(std::pair<int, double> item);
	void div(int index, double value);

	// UTILITY
	void reset(int len);
	void reset(int len, direction d);
	int size() const;
	void clean();
	void clear();
	void print() const;

	// OTHER
	std::vector<sparse_matrix_elem> getElements(direction d, int x = 0) const;
	double l2_norm() const;
	double sum() const;
	double dot(const sparse_vector &other) const;
};

#endif //MPI_MATRICES_SPARSE_VECTOR_H
