#include "Matrix.h"
#include <iostream>
#include <iomanip>
#include <climits>


int main()
{
	size_t n;
	cout << "Enter amount of equations (n equations of n variables)." << endl;
	while (!(cin >> n))
	{
		cin.clear();
		cin.ignore(INT_MAX, '\n');
		cout << "Invalid input, try again." << endl;
		cout << "Enter amount of equations (n equations of n variables)." << endl;
	}
	Matrix m = Matrix(n, n + 1);
	if(!m.FillFromStdin())
	{
		cout << "Invalid input." << endl;
		return 1;
	}
	Matrix res = Matrix(n, n + 1);
	m.Gem(res);

	return 0;
}