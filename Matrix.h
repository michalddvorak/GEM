#pragma once

#include <cstdlib> //For size_t
#include <iostream>
#define M_EPSILON 10e-15

using namespace std;

class Matrix
{
public:
	Matrix(size_t rows, size_t cols);

	Matrix(const Matrix &matrix);

	Matrix(Matrix &&matrix) noexcept;

	virtual ~Matrix();

	double &operator ()(size_t row, size_t col);

	const double &operator ()(size_t row, size_t col) const;

	Matrix &operator =(const Matrix &m2);

	Matrix operator -() const;

	Matrix operator +(const Matrix &m2) const;

	Matrix &operator +=(const Matrix &m2);

	Matrix operator *(double k) const;

	Matrix &operator *=(double k);

	Matrix operator -(const Matrix &m2) const;

	Matrix &operator -=(const Matrix &m2);

	Matrix operator *(const Matrix &m2) const;

	bool operator ==(const Matrix &m2) const;

	Matrix Transpose();

	bool Inverse(Matrix & outres) const;

	bool Gem(Matrix & outres) const;

	bool FillFromStdin();

	void FillFromArray(const double *);

	void FillWithRandom(double lo, double hi);

	void FillMainDiagonal();

	void FillAntiDiagonal();

	void FillWithSingleValue(double val);

	size_t GetRows() const;

	size_t GetCols() const;

	double * GetMainDiagonal() const;

	double * GetAntiDiagonal() const;

	double * GetRow(int m) const;

	double * GetCol(int n) const;

	bool IsUpperTriangular() const;

	bool IsLowerTriangular() const;

	double Det() const;

	friend ostream& operator << (ostream& stream, const Matrix& matrix);

	Matrix();

private:
	double * internal_array;
	size_t rows, cols;
	size_t totalsize;
	bool issq() const;
};

//scalar multiplication k * m
Matrix operator *(double k, const Matrix &m);