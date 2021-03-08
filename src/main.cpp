#include "matrix.h"
#include <iostream>


int main(int argc, char *argv[]){
	Vec<3> V({2, 4, 2});
	Vec<3> U({3, 3, 3});

	Vec<3> CP = V * U;

	V[2] = 5;

	Mat<3, 3> A({
		{1, 3, 3},
		{1, 4, 3},
		{1, 3, 4}
	});

	Mat<4, 4> B({
		{2,  4,   7, 2},
		{4, 67,   7, 6},
		{6, -5,  12, 8},
		{9, -2, -25, 10}
	});

	std::cout << "Determinent: " << A.Determinant() << std::endl;

	std::cout << "Identity matrix: " << std::endl;
	
	Mat<3, 3, double>(A.Inverse() * A).print();

	return 0;
}