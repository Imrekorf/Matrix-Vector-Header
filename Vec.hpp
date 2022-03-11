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

namespace MVH {

/**
 * @brief Specifies if the angle is radian or degrees
 */
enum AngleType{ 
	Radians,
	Degrees
};

class Quaternion;

/**
 * @brief Class for vector arithmatic and logic
 * 
 * @tparam S The size of the vector
 * @tparam T The datatype of the vector
 */
template<int S, typename T = double>
class Vec{
public:
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

	Vec(){
		for(unsigned int i = 0; i < S; i++){
			Varr[i] = 0;
		}
	}

	/**
	 * @brief Construct a new Vec object
	 * 
	 * @param X The value for X, default 0
	 * @param Y The value for Y, default 0
	 * @param Z The value for Z, default 0
	 * @param W The value for W, default 0
	 */
	Vec(const T X, const T Y = 0, const T Z = 0, const T W = 0) : x(X), y(Y), z(Z), w(W){
		static_assert(std::is_arithmetic<T>::value, "Values must be numeric");
		for (unsigned int i = 4; i < S; i++)
		{
			Varr[i] = 0;
		}
	}

	/**
	 * @brief Construct a new Vec object based on an initializer list
	 * 
	 * @param v the list of values
	 */
	Vec(std::initializer_list<T> v){
		static_assert(std::is_arithmetic<T>::value, "Values must be numeric");
		for(auto vi = v.begin(); vi != v.end(); vi++){
			Varr[(int)(vi - v.begin())] = *vi;
		}
	}

	/**
	 * @brief Construct a new Vec object by copying v into the vector object
	 * 
	 * @param v the vector object
	 */
	Vec(const Vec<S, T> &v){
		for(int i = 0; i < S; i++){
			Varr[i] = v.Varr[i];
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
	
	/** @brief Sets the value of element i
	 *  @param i the index of the element
	 *  @return T reference to the value of the element */
	T& operator[](int i){return Varr[i];}
	/** @brief Gets the value of element i
	 *  @param i the index of the element
	 *  @return T the value of the element */
	T operator[](int i) const {return Varr[i];}

	double length() const;
	double lengthsquared() const;
	template<typename N = double> double distance(const Vec<S, N> &U) const;
	template<typename N = double> Vec<S, T> crossproduct(const Vec<S, N> &U) const;
	template<AngleType degrees = Degrees, typename N = double> double angle(const Vec<S, N> &U) const;

	
		  Vec<S, T>& normalize();
	const Vec<S, T>  getNormalized() const;
		  Vec<3, T>& rotate(const Quaternion& q);
		  Vec<3, T>  getRotated(const Quaternion& q) const;

	void print();	
};

//? definitions

//! addition
/**
 * @brief Adds vector U to the vector object
 * V + U
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @tparam N the datatype of vector U
 * @param U the 2nd vector to add from
 * @return Vec<S, T>& reference to the added vector object
 */
template<int S, typename T>
template<typename N>
Vec<S, T>& Vec<S, T>::operator+=(const Vec<S, N> &U){
	for(unsigned int i = 0; i < S; Varr[i] += U.Varr[i], i++){}
	return *this;
}

/**
 * @brief Adds vector U to the vector object
 * V + U
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @tparam N the datatype of vector U
 * @param U the 2nd vector to add from
 * @return const Vec<S, T> The added vector object
 */
template<int S, typename T>
template<typename N> 
const Vec<S, T> Vec<S, T>::operator+(const Vec<S, N> &U) const{
	return Vec<S, T>(*this) += U;
}

//! subtraction

//* Subtracts two vectors
/**
 * @brief Subtracts vector U from the vector object
 * V - U
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @tparam N the datatype of vector U
 * @param U the 2nd vector to substract from
 * @return Vec<S, T>& reference to the substracted vector object
 */
template<int S, typename T>
template<typename N>
Vec<S, T>& Vec<S, T>::operator-=(const Vec<S, N> &U){
	for(unsigned int i = 0; i < S; Varr[i] -= U.Varr[i], i++){}
	return *this;
}

//* Subtracts two vectors
/**
 * @brief Subtracts vector U from the vector object
 * V - U
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @tparam N the datatype of vector U
 * @param U the 2nd vector to substract from
 * @return const Vec<S, T> The substracted vector object
 */
template<int S, typename T>
template<typename N>
const Vec<S, T> Vec<S, T>::operator-(const Vec<S, N> &U) const{
	return Vec<S, T>(*this) -= U;
}

//! multiplication

/**
 * @brief Calculates the scalar of the vector object inplace
 * A * r = U
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @tparam N the datatype of vector U
 * @param r the scalar value
 * @return Vec<S, T>& reference to the scaled vector object
 */
template<int S, typename T>
template<typename N>
Vec<S, T>& Vec<S, T>::operator*=(const N r){
	static_assert(std::is_arithmetic<N>::value, "Type must be numeric or Vector<T>");
	for(unsigned int i = 0; i < S; Varr[i] *= r, i++){}
	return *this;
}

/**
 * @brief Calculates the scalar of the vector object
 * A * r = U
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @tparam N the datatype of vector U
 * @param r the scalar value
 * @return const Vec<S, T> scaled vector object
 */
template<int S, typename T>
template<typename N>
const Vec<S, T> Vec<S, T>::operator*(const N r) const{
	static_assert(std::is_arithmetic<N>::value, "Type must be numeric or Vector<T>");
	return Vec<S, T>(*this) *= r;
}


/**
 * @brief Calculates the scalar of the vector object
 * r * A = U
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @tparam N the datatype of vector U
 * @param r the scalar value
 * @param V the vector to scale
 * @return const Vec<S, T> scaled vector object
 */
template<int S, typename T, typename N>
const Vec<S, T> operator*(const N &r, const Vec<S, T> &V) {
	static_assert(std::is_arithmetic<N>::value, "Type must be numeric or Vector<T>");
	return Vec<S, T>(V) *= r;
}

/**
 * @brief Calculates the dot product between the vector object and vector U
 * A * U = V 
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @tparam N the datatype of vector U
 * @param U the 2nd vector to calculate the dot product between
 * @return double 
 */
template<int S, typename T>
template<typename N>
double Vec<S, T>::operator*(const Vec<S, N> &U) const{
	T j = 0; 
	for(unsigned int i = 0; i < S; j += Varr[i] * U.Varr[i], i++){}
	return j;
}

/**
 * @brief Calculates the scalar of the vector object inplace.
 * A / r = U
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @tparam N the datatype of vector U
 * @param r the scalar value
 * @return Vec<S, T>& reference to the scaled vector.
 */
template<int S, typename T>
template<typename N>
Vec<S, T>& Vec<S, T>::operator/=(const N r){
	static_assert(std::is_arithmetic<N>::value, "Type must be numeric or Vector<T>");
	for(unsigned int i = 0; i < S; Varr[i] /= r, i++){}
	return *this;
}

/**
 * @brief Calculates the scalar of the vector object vector
 * A / r = U
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @tparam N the datatype of vector U
 * @param r the scalar value
 * @return const Vec<S, T> the scaled vector.
 */
template<int S, typename T>
template<typename N>
const Vec<S, T> Vec<S, T>::operator/(const N r) const{
	static_assert(std::is_arithmetic<N>::value, "Type must be numeric or Vector<T>");
	return Vec<S, T>(*this) /= r;
}

//! logical
/**
 * @brief Checks if this vecgtor and vector U are equal
 * 
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @tparam N the datatype of vector U
 * @param U the 2nd vector to compare to
 * @return true if the vectors are equal
 * @return false if the vectors are not equal
 */
template<int S, typename T>
template<typename N>
bool Vec<S, T>::operator==(const Vec<S, N> &U) const{
	bool t = 1;
	for(unsigned int i = 0; i < S && t; t = Varr[i] == U.Varr[i], i++){}
	return t;
}

/**
 * @brief Checks if the vector object and vector U are not equal
 * 
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @tparam N the datatype of vector U
 * @param U the 2nd vector to compare to
 * @return true if the vectors are not equal
 * @return false if the vectors are equal
 */
template<int S, typename T>
template<typename N>
bool Vec<S, T>::operator!=(const Vec<S, N> &U) const{
	return !(*this == U);
}

//! Code based operators

//* assigns a vector to another
/**
 * @brief Copies vector U of length S2 to the vector object of size S
 * The smallest vector defines how many variables get copied. 
 * Starts from variable X.
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @tparam S2 the size of vector U
 * @tparam N the datatype of vector U
 * @param U the vector to copy from
 * @return Vec<S, T>& reference to the assigned vector object.
 */
template<int S, typename T>
template<int S2, typename N> 
Vec<S, T>& Vec<S, T>::operator=(const Vec<S2, N> &U){
	for(unsigned int i = 0; i < S && i < S2; i++){
		Varr[i] = U.Varr[i];
	}
	return *this;
}

/**
 * @brief Casts the vector object from type T to type Cast
 * 
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @tparam Cast the new datatype of the vector
 * @return Vec<S, Cast> The new vector containing the casted values.
 */
template<int S, typename T>
template<typename Cast> 
Vec<S, T>::operator Vec<S, Cast> () const{
	Vec<S, Cast> j; 
	for(unsigned int i = 0; i < S; j[i] = static_cast<Cast>(Varr[i]), i++){}
	return j;
}

//! mathematical functions

/**
 * @brief Calculates the length of the vector object but does not calculate the root of it.
 * V[i]^2 + V[i+1]^2 + ...
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @return double the length of the vector squared
 */
template<int S, typename T>
double Vec<S, T>::lengthsquared() const{
	double j = 0;
	for(unsigned int i = 0; i < S; j += std::pow(Varr[i], 2), i++){}
	return j;
}

/**
 * @brief Calculates the length of the vector object
 * √(V[i]^2 + V[i+1]^2 + ...)
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @return double the length of the vector
 */
template<int S, typename T>
double Vec<S, T>::length() const{
	return std::sqrt(lengthsquared());
}

/**
 * @brief Calculates the angle between two vector objects
 * acos((V • U) / ||V * U||) 
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @tparam N the datatype vector U
 * @tparam degrees specifies if the returned value should be in radians or degrees.
 * @param U the 2nd vector to calculate the angle between
 * @return double 
 */
template<int S, typename T>
template<AngleType degrees, typename N> 
double Vec<S, T>::angle(const Vec<S, N> &U) const{
	return std::acos(
		this->dotproduct(U) / std::sqrt(this->lengthsquared() * U.lengthsquared())
	) * ((bool)degrees ? (180.0 / PI) : 1);
}

/**
 * @brief Calculates the cross product between the vector object and vector U. 
 * Supported vector sizes are between 2 and 3.
 * A x B
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @tparam N the datatype vector U
 * @param U the 2nd vector to calculate the cross product with
 * @return Vec<S, T> the new vector containing the cross product
 */
template<int S, typename T>
template<typename N>
Vec<S, T> Vec<S, T>::crossproduct(const Vec<S, N> &U) const{
	static_assert(S == 3 || S == 2, "Cross Product Vectors must be of size 2, 3");
	Vec<3, T> W({Varr[2]*U[3] - Varr[3]*U[2], 
				 Varr[3]*U[1] - Varr[1]*U[3], 
				 Varr[1]*U[2] - Varr[2]*U[1]}); 
	return W;
}

/**
 * @brief Calculates the distance between the vector object and vector U. 
 * ||V - U||
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @tparam N the datatype vector U
 * @param U 2nd vector to calculate distance to
 * @return double the distance between the two vectors.
 */
template<int S, typename T>
template<typename N>
double Vec<S, T>::distance(const Vec<S, N> &U) const{
	return (*this - U).length();
}

// ? output functions

/**
 * @brief Writers a vector to stream
 * 
 * @tparam S The size of the vector
 * @tparam T The datatype of the vector
 * @param ostream the stream to write to
 * @param V the vector object to write to stream
 * @return std::ostream& reference to the stream containing the vector values.
 */
template<int S, typename T>
std::ostream& operator<<(std::ostream &ostream, const Vec<S, T> &V){
	//char arr[4] = {'x', 'y', 'z', 'w'};
	for(unsigned int i = 0; i < S;){
		ostream << V.Varr[i];
		if(++i < S){ostream << ", ";};
	}
	return ostream;
}

/** @brief prints vector values
 */
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


// ? Rotation functions

/**
 * @brief Specialized vector known as a quaternion, used for rotation calculations.
 * 
 */
class Quaternion : public Vec<4, double> {
public:
	Quaternion() : Vec<4, double>(1, 0, 0, 0) {}
	Quaternion(double w, double x, double y, double z) : Vec<4, double>(x, y, z, w) {}

	/**
	 * @brief Calculates the Hamilton product between two quaternions
	 * @param q the other quaternion to use for the product
	 * @return Quaternion the product of the two quaternions
	 */
	const Quaternion getHProduct(const Quaternion q) const{
		// Quaternion multiplication is defined by:
		//     (Q1 * Q2).w = (w1w2 - x1x2 - y1y2 - z1z2)
		//     (Q1 * Q2).x = (w1x2 + x1w2 + y1z2 - z1y2)
		//     (Q1 * Q2).y = (w1y2 - x1z2 + y1w2 + z1x2)
		//     (Q1 * Q2).z = (w1z2 + x1y2 - y1x2 + z1w2
		return Quaternion(
			w*q.w - x*q.x - y*q.y - z*q.z,  // new w
			w*q.x + x*q.w + y*q.z - z*q.y,  // new x
			w*q.y - x*q.z + y*q.w + z*q.x,  // new y
			w*q.z + x*q.y - y*q.x + z*q.w); // new z
	}
	/**
	 * @brief creates a conjugated version of the quaternion
	 * @return Quaternion the conjugated quaternion
	 */
	const Quaternion getConjugate() const {
		return Quaternion(w, -x, -y, -z);
	}
};

/**
 * @brief Rotates the vector object inplace by the quaternion vector
 * 
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @param q the quaternion rotation vector
 * @return Vec<3, T>& reference to the rotated vector object.
 */
template<int S, typename T>
Vec<3, T>& Vec<S, T>::rotate(const Quaternion& q){
	static_assert(S == 3, "Quaternion rotation needs a vector of size 3");
	// http://www.cprogramming.com/tutorial/3d/quaternions.html
	// http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/transforms/index.htm
	// http://content.gpwiki.org/index.php/OpenGL:Tutorials:Using_Quaternions_to_represent_rotation
	// ^ or: http://webcache.googleusercontent.com/search?q=cache:xgJAp3bDNhQJ:content.gpwiki.org/index.php/OpenGL:Tutorials:Using_Quaternions_to_represent_rotation&hl=en&gl=us&strip=1

	// P_out = q * P_in * conj(q)
	// - P_out is the output vector
	// - q is the orientation quaternion
	// - P_in is the input vector (a*aReal)
	// - conj(q) is the conjugate of the orientation quaternion (q=[w,x,y,z], q*=[w,-x,-y,-z])
	Quaternion p(0, x, y, z);

	// quaternion multiplication: q * p, stored back in p
	p = q.getHProduct(p);

	// quaternion multiplication: p * conj(q), stored back in p
	p = p.getHProduct(q.getConjugate());

	// p quaternion is now [0, x', y', z']
	x = p.x;
	y = p.y;
	z = p.z;
	return *this;
}

/**
 * @brief Rotates a copy of the vector object by quaternion and returns it.
 * 
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @param q the quaternion rotation vector
 * @return Vec<3, T> the rotated vector
 */
template<int S, typename T>
Vec<3, T> Vec<S, T>::getRotated(const Quaternion& q) const{
	static_assert(S == 3, "Quaternion rotation needs a vector of size 3");
	Vec<3, T> r(x, y, z);
	r.rotate(q);
	return r; 
}

/**
 * @brief Normalizes the vector object inplace
 * 
 * @tparam S the size of the vector
 * @tparam T the datatype of the vector
 * @return Vec<S, T>& reference to the normalized vector object.
 */
template<int S, typename T>
Vec<S, T>& Vec<S, T>::normalize(){
	*this /= lengthsquared();
	return *this;
}

/**
 * @brief Creates a new vector object and normalizes it.
 * 
 * @tparam S The size of the vector
 * @tparam T The datatype of the vector
 * @return const Vec<S, T> The unit vector
 */
template<int S, typename T>
const Vec<S, T> Vec<S, T>::getNormalized() const{
	Vec<S, T> r(*this);
	r.normalize();
	return r;
}

}

#undef PI 