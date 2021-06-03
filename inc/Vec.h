/* ---------------------------------------------------------------------------
MIT License

Copyright (c) 2021 Imre Korf

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
--------------------------------------------------------------------------- */

#pragma once

#include <memory>
#include <type_traits>
#include <cmath>
#include <iostream>

#define PI 3.14159265

enum AngleType{ 
	Radians,
	Degrees
};

template<int S, typename T = double>
class Vec{
private:
	union
	{
		struct {
			T x;
			T y;
			T z;
			T w;
		};
		T Varr[S];
	};

public:
	Vec() {
		static_assert(std::is_arithmetic<T>::value, "Values must be numeric");
		for (unsigned int i = 0; i < S; i++)
		{
			Varr[i] = 0;
		}
	}

	Vec(const T X, const T Y = 0, const T Z = 0, const T W = 0) : x(X), y(Y), z(Z), w(W){
		static_assert(std::is_arithmetic<T>::value, "Values must be numeric");
		for (unsigned int i = 3; i < S; i++)
		{
			Varr[i] = 0;
		}
	}

	Vec(std::initializer_list<T> v){
		static_assert(std::is_arithmetic<T>::value, "Values must be numeric");
		for(auto vi = v.begin(); vi != v.end(); vi++){
			Varr[(int)(vi - v.begin())] = *vi;
		}
	}

	~Vec(){}

	// declarations

	template<typename N = double> const Vec<S, T>  operator+ (const Vec<S, N> &U) const;
	template<typename N = double> 		Vec<S, T>& operator+=(const Vec<S, N> &U);					
	template<typename N = double> const Vec<S, T>  operator- (const Vec<S, N> &U) const;	
	template<typename N = double> 		Vec<S, T>& operator-=(const Vec<S, N> &U);
	template<typename N = double> 		double	   operator* (const Vec<S, N> &U) const;	// dot product
	template<typename N = double> const Vec<S, T>  operator* (const N R)		  const;	// scalar
	template<typename N = double> const Vec<S, T>  operator/ (const N R)		  const;	// scalar
	template<typename N = double> 		Vec<S, T>& operator*=(const N R);					// scalar
	template<typename N = double> 		Vec<S, T>& operator/=(const N R);					// scalar
	template<int S2, 
			 typename N = double> 	    Vec<S, T>& operator= (const Vec<S2, N> &U);	

	template<typename N = double> 		bool 	   operator==(const Vec<S, N> &U) const;
	template<typename N = double> 		bool 	   operator!=(const Vec<S, N> &U) const;
	template<typename Cast> 		   			   operator Vec<S, Cast> () const;
	
	T& operator[](int i){return Varr[i];}
	T operator[](int i) const {return Varr[i];}

	double length() const;
	double lengthsquared() const;
	template<typename N = double> double distance(const Vec<S, N> &U) const;
	template<typename N = double> Vec<S, T> crossproduct(const Vec<S, N> &U) const;
	template<AngleType degrees = Degrees, typename N = double> double angle(const Vec<S, N> &U) const;

	void print();

	T X() const{return x;}
	void X(T val){x = val;}

	T Y() const{return y;}
	void Y(T val){y = val;}

    T Z() const{return z;}
	void Z(T val){z = val;}

	T W() const{return w;}
	void W(T val){w = val;}
	
};

//? definitions

//! addition
//* Adds two vectors
template<int S, typename T>
template<typename N>
Vec<S, T>& Vec<S, T>::operator+=(const Vec<S, N> &U){
	for(unsigned int i = 0; i < S; Varr[i] += U.Varr[i], i++){}
	return *this;
}

//* Adds two vectors
template<int S, typename T>
template<typename N> 
const Vec<S, T> Vec<S, T>::operator+(const Vec<S, N> &U) const{
	return Vec<S, T>(*this) += U;
}

//! subtraction

//* Subtracts two vectors
template<int S, typename T>
template<typename N>
Vec<S, T>& Vec<S, T>::operator-=(const Vec<S, N> &U){
	for(unsigned int i = 0; i < S; Varr[i] -= U.Varr[i], i++){}
	return *this;
}

//* Subtracts two vectors
template<int S, typename T>
template<typename N>
const Vec<S, T> Vec<S, T>::operator-(const Vec<S, N> &U) const{
	return Vec<S, T>(*this) -= U;
}

//! multiplication

//* Calculates the Scalar
//* Ar = U
template<int S, typename T>
template<typename N>
Vec<S, T>& Vec<S, T>::operator*=(const N r){
	static_assert(std::is_arithmetic<N>::value, "Type must be numeric or Vector<T>");
	for(unsigned int i = 0; i < S; Varr[i] *= r, i++){}
	return *this;
}

//* Calculates the Scalar
//* Ar = U
template<int S, typename T>
template<typename N>
const Vec<S, T> Vec<S, T>::operator*(const N r) const{
	static_assert(std::is_arithmetic<N>::value, "Type must be numeric or Vector<T>");
	return Vec<S, T>(*this) *= r;
}

//* Calculates the Scalar
//* rA = U
template<int S, typename T, typename N>
const Vec<S, T> operator*(const N &r, const Vec<S, T> &V) {
	static_assert(std::is_arithmetic<N>::value, "Type must be numeric or Vector<T>");
	return Vec<S, T>(V) *= r;
}

//* Calculates the dot product
//* A • B = V[i] * U[i] + V[i+1] * U[i+1] + ...
template<int S, typename T>
template<typename N>
double Vec<S, T>::operator*(const Vec<S, N> &U) const{
	T j = 0; 
	for(unsigned int i = 0; i < S; j += Varr[i] * U.Varr[i], i++){}
	return j;
}

//* Calculates the Scalar
//* Ar = U
template<int S, typename T>
template<typename N>
Vec<S, T>& Vec<S, T>::operator/=(const N r){
	static_assert(std::is_arithmetic<N>::value, "Type must be numeric or Vector<T>");
	for(unsigned int i = 0; i < S; Varr[i] /= r, i++){}
	return *this;
}

//* Calculates the Scalar
//* A / r = U
template<int S, typename T>
template<typename N>
const Vec<S, T> Vec<S, T>::operator/(const N r) const{
	static_assert(std::is_arithmetic<N>::value, "Type must be numeric or Vector<T>");
	return Vec<S, T>(*this) /= r;
}

//! logical

template<int S, typename T>
template<typename N>
bool Vec<S, T>::operator==(const Vec<S, N> &U) const{
	bool t = 1;
	for(unsigned int i = 0; i < S && t; t = Varr[i] == U.Varr[i], i++){}
	return t;
}

template<int S, typename T>
template<typename N>
bool Vec<S, T>::operator!=(const Vec<S, N> &U) const{
	return !(*this == U);
}

//! Code based operators

//* assigns a vector to another
template<int S, typename T>
template<int S2, typename N> 
Vec<S, T>& Vec<S, T>::operator=(const Vec<S2, N> &U){
	for(unsigned int i = 0; i < S; i++){
		Varr[i] = U.Varr[i];
	}
	return *this;
}

template<int S, typename T>
template<typename Cast> 
Vec<S, T>::operator Vec<S, Cast> () const{
	Vec<S, Cast> j; 
	for(unsigned int i = 0; i < S; j[i] = static_cast<Cast>(Varr[i]), i++){}
	return j;
}

//! mathematical functions

//* returns the length of the vector squared
//* V[i]^2 + V[i+1]^2 + ...
template<int S, typename T>
double Vec<S, T>::lengthsquared() const{
	double j = 0;
	for(unsigned int i = 0; i < S; j += std::pow(Varr[i], 2), i++){}
	return j;
}

//* returns the length of the vector
//* √(V[i]^2 + V[i+1]^2 + ...)
template<int S, typename T>
double Vec<S, T>::length() const{
	return std::sqrt(lengthsquared());
}

//* returns the angle of the vector
//* acos((V • U) / ||V * U||) 
template<int S, typename T>
template<AngleType degrees, typename N> 
double Vec<S, T>::angle(const Vec<S, N> &U) const{
	return std::acos(
		this->dotproduct(U) / std::sqrt(this->lengthsquared() * U.lengthsquared())
	) * ((bool)degrees ? (180.0 / PI) : 1);
}

//* returns the cross product of the vector
//* A x B
template<int S, typename T>
template<typename N>
Vec<S, T> Vec<S, T>::crossproduct(const Vec<S, N> &U) const{
	static_assert(S == 3 || S == 2, "Cross Product Vectors must be of size 2, 3");
	Vec<3, T> W({Varr[2]*U[3] - Varr[3]*U[2], 
				 Varr[3]*U[1] - Varr[1]*U[3], 
				 Varr[1]*U[2] - Varr[2]*U[1]}); 
	return W;
}

//* returns the distance between two vectors
//* ||V - U||
template<int S, typename T>
template<typename N>
double Vec<S, T>::distance(const Vec<S, N> &U) const{
	return (*this - U).length();
}

//* prints vector values
template<int S, typename T>
void Vec<S, T>::print(){
	char arr[4] = {'x', 'y', 'z', 'w'};
	for(unsigned int i = 0; i < S;){
		if(i < 4){
			std::cout << arr[i];
		}
		else{
			std::cout << i;
		}
		std::cout << ": " << Varr[i];
		if(++i < S){std::cout << ", ";};
	}
	std::cout << std::endl;
}

#undef PI 