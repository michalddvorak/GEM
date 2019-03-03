#include "Matrix.h"
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <limits>
#include <cstring>
#include <random>
#include <cfloat>

using namespace std;

Matrix::Matrix()
{
	this->rows           = this->cols = this->totalsize = 0;
	this->internal_array = nullptr;
}

Matrix::Matrix(size_t rows, size_t cols)
{
	this->rows           = rows;
	this->cols           = cols;
	this->totalsize      = rows * cols;
	this->internal_array = new double[rows * cols];
}


Matrix::Matrix(const Matrix &matrix)
{
	//cout << "Called copy" << endl;
	this->rows           = matrix.rows;
	this->cols           = matrix.cols;
	this->totalsize      = matrix.totalsize;
	this->internal_array = new double[rows * cols];
	memcpy(this->internal_array, matrix.internal_array, totalsize * sizeof(double));
}

Matrix::Matrix(Matrix &&matrix) noexcept
{
	//cout << "Called move" << endl;
	this->rows           = matrix.rows;
	this->cols           = matrix.cols;
	this->totalsize      = matrix.totalsize;
	this->internal_array = matrix.internal_array;
}

Matrix::~Matrix()
{
	delete[](internal_array);
}

void Matrix::FillWithSingleValue(double val)
{
	fill_n(internal_array, rows * cols, val);
}

void Matrix::FillWithRandom(double lo, double hi)
{
	uniform_real_distribution<double> unif(lo, hi);
	default_random_engine             re;
	for (size_t                       i = 0; i < totalsize; i++)
		internal_array[i] = unif(re);
}

bool Matrix::FillFromStdin()
{
	cout << "Enter matrix: (" << rows << " × " << cols << "):" << endl;

	for (size_t i = 0; i < rows * cols; i++)
	{
		if (!(cin >> internal_array[i]))
		{
			cout << "Invalid input.\n";
			return false;
		}
	}
	return true;

}

double &Matrix::operator ()(size_t row, size_t col)
{
	if (row < 0 || row >= rows)
		throw out_of_range("Rows out of range.");
	if (col < 0 || col >= cols)
		throw out_of_range("Cols out of range.");
	return internal_array[cols * row + col];
}

const double &Matrix::operator ()(size_t row, size_t col) const
{
	if (row < 0 || row >= rows)
		throw out_of_range("Rows out of range.");
	if (col < 0 || col >= cols)
		throw out_of_range("Cols out of range.");
	return internal_array[cols * row + col];
}

Matrix Matrix::operator +(const Matrix &m2) const
{
	Matrix res = *this; //udela deep copy, zavola se konstruktor Matrix(const& Matrix)
	//to same jako Matrix res = Matrix(*this)
	res += m2;
	return res;
}

Matrix &Matrix::operator +=(const Matrix &m2)
{
	if (m2.rows != this->rows || m2.cols != this->cols)
		throw length_error("Matrix dimensions mismatch");

	for (size_t i = 0; i < totalsize; i++)
		this->internal_array[i] += m2.internal_array[i];
	return *this;
}

Matrix Matrix::operator -() const
{
	Matrix res = *this;
	res *= -1;
	return res;
}

Matrix &Matrix::operator *=(double k)
{
	for (size_t i = 0; i < totalsize; i++)
		this->internal_array[i] *= k;
	return *this;
}

Matrix operator *(double k, const Matrix &m)
{
	Matrix res = m;
	res *= k;
	return res;
}

Matrix Matrix::operator *(double k) const
{
	return k * *this;
}

Matrix &Matrix::operator -=(const Matrix &m)
{
	return *this += -m;
}

Matrix Matrix::operator -(const Matrix &m) const
{
	return (*this) + (-m);
}

Matrix &Matrix::operator =(const Matrix &m2)
{
	if (&m2 == this)
		return *this;
	if (m2.totalsize == this->totalsize)
		memcpy(this->internal_array, m2.internal_array, totalsize * sizeof(double));
	else
	{
		delete[](this->internal_array);
		this->internal_array = new double[m2.totalsize];
		memcpy(this->internal_array, m2.internal_array, m2.totalsize * sizeof(double));
	}
	this->rows      = m2.rows;
	this->cols      = m2.cols;
	this->totalsize = m2.totalsize;
	return *this;

}

Matrix Matrix::operator *(const Matrix &m2) const
{
	if (this->cols != m2.rows)
	{
		throw length_error("Matrix dimensions mismatch");
	}
	Matrix   res = Matrix(this->rows, m2.cols);
	for (int i   = 0; i < this->rows; i++)
	{
		for (int j = 0; j < m2.cols; j++)
		{
			double   sum = 0;
			for (int r   = 0; r < this->cols; r++)
			{
				sum += (*this)(i, r) * m2(r, j);
			}
			res(i, j) = sum;
		}
	}
	return res;
}

bool Matrix::operator ==(const Matrix &m2) const
{
	if (m2.rows != this->rows || m2.cols != this->cols)
		return false;
	for (size_t i = 0; i < totalsize; i++)
		if (internal_array[i] != m2.internal_array[i])
			return false;
	return true;
}

/**
 * Transposes a matrix
 * @return Transposed matrix.
 */
Matrix Matrix::Transpose()
{
	Matrix mat = *this;
	mat.cols = this->rows;
	mat.rows = this->cols;
	for (size_t i = 0; i < mat.rows; i++)
	{
		for (size_t j = 0; j < mat.cols; j++)
		{
			mat(i, j) = (*this)(j, i);
		}
	}
	return mat;
}

/**
 * Fill a square matrix with 1's on the main diagonal.
 */
void Matrix::FillMainDiagonal()
{
	if (!this->issq())
		throw length_error("Matrix must be n × n");
	this->FillWithSingleValue(0);
	for (size_t i     = 0; i < rows; i++)
		(*this)(i, i) = 1;
}

/**
 * Determines if this matrix is a square matrix
 * @return true if matrix is square, false otherwise
 */
bool Matrix::issq() const
{
	return this->rows == this->cols;
}

/**
 * Calculates the determinant of a square matrix.
 * @return the determinant of this matrix.
 */
double Matrix::Det() const
{
	if (!this->issq())
		throw length_error("Matrix must be n × n");
	int      n   = (int)this->rows;
	Matrix   mat = *this;
	for (int k   = 0; k < n - 1; k++)
	{
		int r = k;
		while (abs(mat(r, k)) < M_EPSILON)
		{
			if (r == n - 1)
				return 0;
			else
				r++;
			cout << mat;
		}
		if (r > k)
		{
			for (int q = 0; q < n; q++)
			{
				double tmp = mat(r, q);
				mat(r, q) = mat(k, q);
				mat(k, q) = tmp;
			}
		}
		cout << mat;
		for (int i = k + 1; i < n; i++)
		{
			for (int j = k + 1; j < n; j++)
			{
				mat(i, j) -= (mat(i, k) / mat(k, k)) * mat(k, j);
			}
			mat(i, k)  = 0;
			cout << mat;
		}
	}
	double   det = 1;
	for (int i   = 0; i < n; i++)
		det *= mat(i, i);
	return det;
}

/**
 * Gaussian elimination. If matrix is in the form n × n + 1 (represents system of n equations of n variables),
 * solve the system and prints steps meanwhile.
 * @param outres reveives a copy of this matrix and all operations of GEM are
 * done on matrix outres
 * @return true if GEM was succesful and found one and only one solution, returns false otherwise
 */
bool Matrix::Gem(Matrix &outres) const
{
	const string l = string(50, '_') + "\n";
	cout << (*this) << l;
	if (this->rows != this->cols - 1)
	{
		cout << "Matrix must be n × n + 1";
		return false;;
	}
	const int &n = (int)this->rows; //amount of equations
	outres = *this;
	//Find nonzero pivot for k-th column
	for (int    k = 0; k < n - 1; k++)
	{
		int r = k;
		while (abs(outres(r, k)) < M_EPSILON)
		{
			if (r == n - 1) //If no nonzero pivot is found for k-th row, matrix is singular
			{
				cout << "Matrix is singular." << endl << outres;
				return false;
			} else
				r++;
			cout << outres << l;
		}
		if (r > k) // If we moved with r, swap lines
		{
			cout << "Swapping row 1 and " << r << endl;
			for (int q = 0; q < n + 1; q++)
			{
				double tmp = outres(r, q);
				outres(r, q) = outres(k, q);
				outres(k, q) = tmp;
			}
			cout << outres << l;
		}
		cout << outres << l;
		for (int i = k + 1; i < n; i++)
		{
			for (int j   = k + 1; j < n + 1; j++)
			{
				outres(i, j) -= (outres(i, k) / outres(k, k)) * outres(k, j);
			}
			outres(i, k) = 0;
			cout << outres << l;
		}
	}
	for (int    k = n - 1; k >= 0; k--)
	{
		if (abs(outres(k, k)) < M_EPSILON)
		{//If we still have zero pivot (might occur)
			outres(k, k) = 0;
			cout << "Matrix is singular." << endl << outres;
			return false;

		}
		for (int i = k - 1; i >= 0; i--)
		{
			outres(i, n) -= (outres(i, k) / outres(k, k)) * outres(k, n);
			outres(i, k) = 0;
			cout << outres << l;
		}
	}
	for (size_t k = 0; k < n; k++)
	{
		outres(k, n) = outres(k, n) / outres(k, k);
		outres(k, k) = 1;
		cout << outres << l;
	}
	cout << "Matrix is solvable:" << endl;
	cout << outres << l;
	return true;
}

/**
 * Gets the amount of rows of this matrix.
 * @return Rows of this matrix.
 */
size_t Matrix::GetRows() const { return this->rows; }

/**
 * Gets the amount of columns of this matrix
 * @return Columns of this matrix.
 */
size_t Matrix::GetCols() const { return this->cols; }

/**
 * Function to send matrix to outputstream.
 * @param os ostream to be printed to
 * @param matrix  matrix to be printed
 * @return reference to os
 */
ostream &operator <<(ostream &os, const Matrix &matrix)
{
	//os << "Matrix[" << matrix.rows << " × " << matrix.cols << "]:" << endl;
	for (size_t i = 0, q = 0; i < matrix.totalsize; q++)
	{
		q &= -(q != matrix.cols); //If q reached cols, reset it to 0
		if (matrix.internal_array[i] >= 0)
			os << " " << matrix.internal_array[i++];
		else
			os << matrix.internal_array[i++];
		if (q == matrix.cols - 1)
			os << endl;
		else
			os << "\t";
	}
	return os;
}
