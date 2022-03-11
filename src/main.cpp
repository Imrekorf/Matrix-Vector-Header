#include "matrix.h"
#include <iostream>

int main(int argc, char *argv[]){
	/*
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
	*/

	
	// hypothesis if a point is inside a rectangle defined by constraints: A, B, C, D with A connected to B and C, B connected to A and D, C connected to A and D, and D connected to C and B.
	// C - D
	// |   |
	// A - B
	// then any vector X from point A will be within the rectangle if the distance is between 0 and the vector (A->B + A->C)

	
	// hypothesis if a point is inside a triangle T defined by constraints A, B, C
	// then any point P is inside the triangle T if both lengths of vectors B->P : V and C->P : S to V are smaller than A->B : D and A->C : E respectivly
	// ||V|| < ||D|| ^ ||S|| < ||E|| => P ∈ T 

	/*
	Vec<2> A(0, 0), B(0, 6), C(6, 0);
	Vec<2> a = B-A; Vec<2> b = C-A; Vec<2> c = C-B;

	Vec<2> P(-1, -1);
	Vec<2> d = P-A; Vec<2> e = P-B; Vec<2> f = P-C;

	if(
		(d.length() + e.length()) < (b.length() + c.length()) &&
		(d.length() + f.length()) < (a.length() + c.length()) &&
		(e.length() + f.length()) < (a.length() + b.length())
	)
		std::cout << "point in triangle" << std::endl;
	else
		std::cout << "point outside triangle" << std::endl;
	*/

	// Question, is the acceleration vector, integrated twice, the same length as the displacement vector?
	const double T = 0.01;
	Vec<3> Acc(1, 2, 3);
	Vec<3> Spe;
	Vec<3> Dis;

	Spe.x = Acc.x * T + 0;
	Spe.y = Acc.y * T + 0;
	Spe.z = Acc.z * T + 0;

	Dis.x = 0 + 0 + 0.5 * ( Acc.x * pow(T, 2) );
	Dis.y = 0 + 0 + 0.5 * ( Acc.y * pow(T, 2) );
	Dis.z = 0 + 0 + 0.5 * ( Acc.z * pow(T, 2) );

	std::cout << "Displacement Vector length: " << Dis.length() << std::endl;
	std::cout << "Acceleration Vector integrated twice length: " << (0 + 0 + 0.5 * ( Acc.length() * pow(T, 2) )) << std::endl;

	return 0;
}