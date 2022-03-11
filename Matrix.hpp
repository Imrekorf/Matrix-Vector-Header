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

#include <cmath>
#include "Vec.hpp"
#include <iostream>

namespace MVH {

/**
 * @brief Base matrix class
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 */
template<unsigned int H, unsigned int W = H, typename T = double>
class MatBase {
private:
	Vec<W, T> M[H];
protected:

public:
	MatBase();
	MatBase(const MatBase<H, W, T> &m1);
	MatBase(std::initializer_list<std::initializer_list<T>> v);
	~MatBase();

	Vec<W, T>& operator[] (int i) {return M[i];}
	Vec<W, T> operator[] (int i) const {return M[i];}
	Vec<W, T> GetRow(const unsigned int i) const;
	Vec<H, T> GetColumn(const unsigned int i) const;
	void SetRow(const unsigned int i, const Vec<W, T> V);
	void SetColumn(const unsigned int i, const Vec<H, T> V);
	MatBase<H, W, T> Transpose(void) const;

	template<typename N>		 		MatBase<H, W , T>& operator*=(const MatBase<H, W , N> &B);		// Matrix multiplication 
	template<int W2, typename N> const	MatBase<H, W2, T>  operator* (const MatBase<W, W2, N> &B) const;// Matrix multiplication
	template<typename N> 		const	MatBase<H, W , T>  operator* (const MatBase<H, W , N> &B) const;// Matrix multiplication
	template<typename N>		 		MatBase<H, W , T>& operator*=(const N &r);						// Scalar multiplication
	template<typename N>		 const	MatBase<H, W , T>  operator* (const N &r) 			  	  const;// Scalar multiplication
	template<typename N>		 		MatBase<H, W , T>& operator+=(const MatBase<H, W, N> &B);		// Matrix addition
	template<typename N>		 const  MatBase<H, W , T>  operator+ (const MatBase<H, W, N> &B)  const;// Matrix addition
	template<typename N>		 		MatBase<H, W , T>& operator-=(const MatBase<H, W, N> &B);		// Matrix subtraction
	template<typename N>		 const  MatBase<H, W , T>  operator- (const MatBase<H, W, N> &B)  const;// Matrix subtraction

	template<typename N> 		 const 	Vec<H, T> 	   	   operator* (const Vec<W, N> &V) 	 	  const;// Vector matrix multiplication
	template<typename Cast> 		   					   operator MatBase<H, W, Cast> () 		  const;
	/**
	 * @brief prints out the values of the matrix
	 */
	void print(){
		for(unsigned int i = 0; i < H; i++){
			for(unsigned int j = 0; j < W; j++){
				std::cout << (*this)[i][j] << " ";
			}
			std::cout << std::endl;
		}
	}
};

//! Data manipulation

/**
 * @brief Construct a new empty MatBase object
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 */
template<unsigned int H, unsigned int W, typename T>
MatBase<H, W, T>::MatBase() {
	static_assert(std::is_arithmetic<T>::value, "Values must be numeric");
}

/**
 * @brief Copy constructor
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @param m1 the MatBase object to copy from
 */
template<unsigned int H, unsigned int W, typename T>
MatBase<H, W, T>::MatBase(const MatBase<H, W, T> &m1) {
	for(unsigned int i = 0; i < H; i++){
		(*this)[i] = m1[i];
	}
}

/**
 * @brief Construct a new MatBase object based on an initializer list
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @param v the initializer list
 */
template<unsigned int H, unsigned int W, typename T>
MatBase<H, W, T>::MatBase(std::initializer_list<std::initializer_list<T>> v){
	static_assert(std::is_arithmetic<T>::value, "Values must be numeric");
	bool err = 0;
	unsigned int size = 0;
	try {
		if(v.size() != H){
			size = v.size();
			throw "Matrix input height does not match templated height\nMatrix height: ";
		}
		for(auto vi = v.begin(); vi != v.end(); vi++){
			for(auto vj = vi->begin(); vj != vi->end(); vj++){
				if((unsigned int)(vj - vi->begin()) >= W){
					err = 1;
					size = vi->size();
					throw "Matrix input width does not match templated width\nMatrix Width: ";
					break;
				}
				(*this)[(int)(vi - v.begin())][(int)(vj - vi->begin())] = *vj;
			}
		}
	} catch (const char* msg) {
		std::cerr << "Error occurred during Matrix constructor: " << std::endl;
		std::cerr << msg << (err ? W : H) << " Input Size: " << size << std::endl;
	}
}

/**
 * @brief Destroy the MatBase object
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 */
template<unsigned int H, unsigned int W, typename T>
MatBase<H, W, T>::~MatBase() {}

/**
 * @brief Returns the row at given position
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @param i the row position
 * @return Vec<W, T> vector containing the row
 */
template<unsigned int H, unsigned int W, typename T>
Vec<W, T> MatBase<H, W, T>::GetRow(const unsigned int i) const{
	return (*this)[i];
}

/**
 * @brief Returns the column at given position
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @param i the column position
 * @return Vec<H, T> vector containing the column
 */
template<unsigned int H, unsigned int W, typename T>
Vec<H, T> MatBase<H, W, T>::GetColumn(const unsigned int i) const{
	Vec<H, T> C;
	for(unsigned int j = 0; j < H; C[j] = (*this)[j][i], j++){}
	return C;
}

/**
 * @brief Set the row of the matrix
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @param i the row position
 * @param V the vector of values
 */
template<unsigned int H, unsigned int W, typename T>
void MatBase<H, W, T>::SetRow(const unsigned int i, const Vec<W, T> V){
	(*this)[i] = V;
}

/**
 * @brief Sets the column of the matrix
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @param i the column position
 * @param V the vector of values
 */
template<unsigned int H, unsigned int W, typename T>
void MatBase<H, W, T>::SetColumn(const unsigned int i, const Vec<H, T> V){
	for(int j = 0; j < H; (*this)[j][i] = V[j], j++){}
}

//Math
/**
 * @brief Transposes the matrix
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @return MatBase<H, W, T> The transposed matrix
 */
template<unsigned H, unsigned int W, typename T>
MatBase<H, W, T> MatBase<H, W, T>::Transpose(void) const{
	MatBase<W, H, T> Mt;
	for(unsigned int i = 0; i < H; i++){
		Mt.SetRow(i, this->GetColumn(i));
	}
	return Mt;
}

//! Operator overloads

/**
 * @brief Multiplies the matrix with another inplace
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @tparam N datatype of the 2nd matrix
 * @param B the 2nd matrix
 * @return MatBase<H, W, T>& reference to the multiplied matrix object
 */
template<unsigned int H, unsigned int W, typename T>	
template<typename N>
MatBase<H, W, T>& MatBase<H, W, T>::operator*=(const MatBase<H, W, N> &B){
	MatBase<H, W, T> C = *this * B;
	*this = C;
	return *this;
}


//* Matrix multiplication
/**
 * @brief Multiplies the matrix with another, allows for different matrix sizes
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @tparam W2 width of the 2nd matrix
 * @tparam N datatype of the 2nd matrix
 * @param B the 2nd matrix
 * @return const MatBase<H, W2, T> the multiplied matrix object
 */
template<unsigned int H, unsigned int W, typename T>
template<int W2, typename N> 		
const MatBase<H, W2, T> MatBase<H, W, T>::operator*(const MatBase<W, W2, N> &B) const{
	MatBase<H, W2, T> C;
	// loops over B's Columns
	for(unsigned int i = 0; i < W2; i++){
		Vec<W, T> ColB = B.GetColumn(i);
		// copy value to new column
		for(unsigned int j = 0; j < H; j++){
			C[j][i] = (*this * ColB)[j];
		}
	}
	return C;
}

//* Matrix multiplication
/**
 * @brief Multiplies the matrix with another
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @tparam N datatype of the 2nd matrix
 * @param B the 2nd matrix
 * @return const MatBase<H, W, T> the multiplied matrix object
 */
template<unsigned int H, unsigned int W, typename T>	
template<typename N>	
const MatBase<H, W, T> MatBase<H, W, T>::operator*(const MatBase<H, W, N> &B) const{
	MatBase<H, W, T> C;
	// loops over B's Columns
	for(unsigned int i = 0; i < W; i++){
		Vec<W, T> ColB = B.GetColumn(i);
		// copy value to new column
		for(unsigned int j = 0; j < H; j++){
			C[j][i] = (*this * ColB)[j];
		}
	}
	return C;
}

/**
 * @brief Multiplies a vector to a matrix
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @tparam N datatype of the vector
 * @param V the vector to multiply with
 * @return const Vec<H, T> 
 */
template<unsigned int H, unsigned int W, typename T>
template<typename N>
const Vec<H, T> MatBase<H, W, T>::operator*(const Vec<W, N> &V) const{
	Vec<H, T> ColA;
	// loops over A's columns and B's column-values
	for(unsigned int j = 0; j < W; j++){
		ColA += V[j] * GetColumn(j);
	}
	return ColA;
}

/**
 * @brief Scales matrix inplace
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @tparam N datatype of the scalar
 * @param r the scalar
 * @return MatBase<H, W , T>& reference to the scaled matrix
 */
template<unsigned int H, unsigned int W, typename T>
template<typename N>		 		
MatBase<H, W , T>& MatBase<H, W, T>::operator*=(const N &r){
	for(unsigned int i = 0; i < H; i++){
		(*this)[i] *= r;
	}
	return *this;
}

/**
 * @brief Scales matrix
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @tparam N datatype of the scalar
 * @param r the scalar
 * @return const MatBase<H, W , T> the scaled matrix
 */
template<unsigned int H, unsigned int W, typename T>
template<typename N>		 		
const MatBase<H, W , T> MatBase<H, W, T>::operator*(const N &r) const{
	return MatBase<H, W, T>(*this) *= r;
}

/**
 * @brief Adds two matrices together inplace
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @tparam N datatype of the 2nd matrix
 * @param B the 2nd matrix
 * @return MatBase<H, W , T>& reference to the matrix
 */
template<unsigned int H, unsigned int W, typename T>
template<typename N>		 		
MatBase<H, W , T>& MatBase<H, W, T>::operator+=(const MatBase<H, W, N> &B){
	for(unsigned int i = 0; i < H; i++){
		(*this)[i] += B[i];
	}
	return *this;
}

/**
 * @brief Adds two matrices together
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @tparam N datatype of the 2nd matrix
 * @param B the 2nd matrix
 * @return const MatBase<H, W , T> the resulting matrix
 */
template<unsigned int H, unsigned int W, typename T>
template<typename N>		 		
const MatBase<H, W , T> MatBase<H, W, T>::operator+(const MatBase<H, W, N> &B) const{
	return MatBase<H, W, T>(*this) += B;
}

/**
 * @brief Subtracts matrix B from the matrix inplace
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @tparam N datatype of the 2nd matrix
 * @param B the 2nd matrix
 * @return MatBase<H, W , T>& reference to the matrix
 */
template<unsigned int H, unsigned int W, typename T>
template<typename N>		 		
MatBase<H, W , T>& MatBase<H, W, T>::operator-=(const MatBase<H, W, N> &B){
	for(unsigned int i = 0; i < H; i++){
		(*this)[i] -= B[i];
	}
	return *this;
}

/**
 * @brief Subtracts matrix B from the matrix
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @tparam N datatype of the 2nd matrix
 * @param B the 2nd matrix
 * @return const MatBase<H, W , T> the resulting matrix
 */
template<unsigned int H, unsigned int W, typename T>
template<typename N>		 		
const MatBase<H, W , T> MatBase<H, W, T>::operator-(const MatBase<H, W, N> &B) const{
	return MatBase<H, W, T>(*this) -= B;
}

/**
 * @brief Casts matrix to different matrix type
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 * @tparam Cast the new datatype of the matrix
 * @return MatBase<H, W, Cast> the casted matrix
 */
template<unsigned int H, unsigned int W, typename T>
template<typename Cast>
MatBase<H, W, T>::operator MatBase<H, W, Cast> () const {
	MatBase<H, W, Cast> CM;
	for(unsigned int i = 0; i < H; CM[i] = static_cast<Vec<W, Cast>>((*this)[i]), i++){}
	return CM;
}



//?								 ?//
//?			Matrix Class		 ?//
//? 							 ?//

/**
 * @brief Matrix class for a non-square matrix
 * 
 * @tparam H height of the matrix
 * @tparam W width of the matrix
 * @tparam T datatype of the matrix
 */
template<unsigned int H, unsigned int W = H, typename T= double>
class Mat : public MatBase<H, W, T> {
public:
	Mat() : MatBase<H, W, T>() {}
	Mat(const Mat<H, W, T> &m1) : MatBase<H, W, T>(m1) {}
	Mat(const MatBase<H, W, T> &m1) : MatBase<H, W, T>(m1) {}
	Mat(std::initializer_list<std::initializer_list<T>> v) : MatBase<H, W, T>(v) {}
	~Mat() {}

	template<int W2, typename N> 		Mat<H, W2, T>& operator*=(const Mat<W, W2, N> &B) 		{return (MatBase<H, W, T>)(*this) *= (MatBase<W, W2, N>)B;} // Matrix multiplication
	template<typename N>		 		Mat<H, W , T>& operator*=(const Mat<H, W , N> &B) 		{return (MatBase<H, W, T>)(*this) *= (MatBase<H, W , N>)B;}	// Matrix multiplication
	template<int W2, typename N> const	Mat<H, W2, T>  operator* (const Mat<W, W2, N> &B) const {return (MatBase<H, W, T>)(*this) *  (MatBase<W, W2, N>)B;}	// Matrix multiplication
	template<typename N> 		 const	Mat<H, W , T>  operator* (const Mat<H, W , N> &B) const {return (MatBase<H, W, T>)(*this) *  (MatBase<H, W , N>)B;}	// Matrix multiplication
	template<typename N>		 		Mat<H, W , T>& operator*=(const N &r)					{return (MatBase<H, W, T>)(*this) *= r;}					// Scalar multiplication
	template<typename N>		 const	Mat<H, W , T>  operator* (const N &r) 			  const {return (MatBase<H, W, T>)(*this) *  r;}					// Scalar multiplication
	template<typename N>		 		Mat<H, W , T>& operator+=(const Mat<H, W , N> &B)		{return (MatBase<H, W, T>)(*this) += (MatBase<H, W , N>)B;} // Matrix addition
	template<typename N>		 const  Mat<H, W , T>  operator+ (const Mat<H, W , N> &B) const	{return (MatBase<H, W, T>)(*this) +  (MatBase<H, W , N>)B;}	// Matrix addition
	template<typename N>		 		Mat<H, W , T>& operator-=(const Mat<H, W , N> &B)		{return (MatBase<H, W, T>)(*this) -= (MatBase<H, W , N>)B;}	// Matrix subtraction
	template<typename N>		 const  Mat<H, W , T>  operator- (const Mat<H, W , N> &B) const	{return (MatBase<H, W, T>)(*this) -  (MatBase<H, W , N>)B;}	// Matrix subtraction

										Mat<H, W , T>  Transpose(void) 					  const {return (MatBase<H, W, T>)(*this).Transpose();}

	template<typename N> 		 const 	Vec<H, T> 	   operator* (const Vec<W, N> &V) 	  const	{return (MatBase<H, W, T>)(*this) * V;}	// Vector matrix multiplication
	template<typename Cast> 		   				   operator Mat<H, W, Cast> () 	  	  const {return (Mat<H, W, Cast>)((MatBase<H, W, T>)(*this));}
};

/**
 * @brief Specialized Square matrix
 * 
 * @tparam H height and width of the matrix
 * @tparam T datatype of the matrix
 */
template<unsigned H, typename T>
class Mat<H, H, T> : public MatBase<H, H, T> {
public:
	Mat() : MatBase<H, H, T>() {}
	Mat(const Mat<H, H, T> &m1) : MatBase<H, H, T>(m1) {}
	Mat(const MatBase<H, H, T> &m1) : MatBase<H, H, T>(m1) {}
	Mat(std::initializer_list<std::initializer_list<T>> v) : MatBase<H, H, T>(v) {}
	~Mat() {}

	template<typename N> 		Mat<H, H, T>& operator*=(const Mat<H, H, N> &B) 		{return (MatBase<H, H, T>)(*this) *= (MatBase<H, H, N>)B;} 	// Matrix multiplication
	template<typename N> const	Mat<H, H, T>  operator* (const Mat<H, H, N> &B) const 	{return (MatBase<H, H, T>)(*this) *  (MatBase<H, H, N>)B;} 	// Matrix multiplication
	template<typename N>		Mat<H, H, T>& operator*=(const N &r)					{return (MatBase<H, H, T>)(*this) *= r;}					// Scalar multiplication
	template<typename N> const	Mat<H, H, T>  operator* (const N &r) 			const 	{return (MatBase<H, H, T>)(*this) *  r;}					// Scalar multiplication
	template<typename N> 		Mat<H, H, T>& operator+=(const Mat<H, H, N> &B)			{return (MatBase<H, H, T>)(*this) += (MatBase<H, H, N>)B;}	// Matrix addition
	template<typename N> const  Mat<H, H, T>  operator+ (const Mat<H, H, N> &B) const	{return (MatBase<H, H, T>)(*this) +  (MatBase<H, H, N>)B;}	// Matrix addition
	template<typename N> 		Mat<H, H, T>& operator-=(const Mat<H, H, N> &B)			{return (MatBase<H, H, T>)(*this) -= (MatBase<H, H, N>)B;}	// Matrix subtraction
	template<typename N> const  Mat<H, H, T>  operator- (const Mat<H, H, N> &B) const	{return (MatBase<H, H, T>)(*this) -  (MatBase<H, H, N>)B;}	// Matrix subtraction
	template<typename Cast> 		   		  operator Mat<H, H, Cast> () 	  	const 	{return (Mat<H, H, Cast>)((MatBase<H, H, T>)(*this));}
	template<typename N> const 	Vec<H, T> 	  operator* (const Vec<H, N> &V) 	const	{return (MatBase<H, H, T>)(*this) * V;}						// Vector matrix multiplication

								Mat<H, H, T>  Transpose(void) 					const 	{return Mat<1, 1, T>(*this);}
	Mat<H-1, H-1, T> SubMatrix(const unsigned int i, const unsigned int j) const;
	Mat<H, H, T> CoFactorMatrix(void) const;
	Mat<H, H, T> Adjoint(void) const;
	Mat<H, H, T> Inverse(void) const;
	double CoFactor(const unsigned int i, const unsigned int j) const;
	double Minor(const unsigned int i, const unsigned int j) const;
	double Determinant(void) const;
};

//! Math

/**
 * @brief Gets a submatrix from the square matrix
 * copies the matrix to a new matrix object with [H-1] size and excludes passed row and column
 * @tparam H height and width of the matrix
 * @tparam T datatype of the matrix
 * @param i the row to skip
 * @param j the column to skip
 * @return Mat<H-1, H-1, T> submatrix of the original matrix
 */
template<unsigned H, typename T>
Mat<H-1, H-1, T> Mat<H, H, T>::SubMatrix(const unsigned int i, const unsigned int j) const {
	Mat<H-1, H-1, T> Ms;
	for(unsigned int i2 = 0, i3 = 0; i2 < H; i2++){
		if(i2 == i) // skip row
			continue;
		for(unsigned int j2 = 0, j3 = 0; j2 < H; j2++){
			if(j2 == j) // skip column
				continue;
			Ms[i3][j3] = (*this)[i2][j2];
			j3++;
		}
		i3++;
	}
	return Ms;
}

/**
 * @brief Calculates the Minor of the matrix
 * 
 * @tparam H height and width of the matrix
 * @tparam T datatype of the matrix
 * @param i the row to skip on the submatrix
 * @param j the column to skip on the submatrix
 * @return double the Minor of the matrix
 */
template<unsigned H, typename T>
double Mat<H, H, T>::Minor(const unsigned int i, const unsigned int j) const {
	return this->SubMatrix(i, j).Determinant();
}

/**
 * @brief Calculates the CoFactor of the matrix
 * 
 * @tparam H height and width of the matrix
 * @tparam T datatype of the matrix
 * @param i the row to skip on the submatrix
 * @param j the column to skip on the submatrix
 * @return double the CoFactor of the matrix
 */
template<unsigned H, typename T>
double Mat<H, H, T>::CoFactor(const unsigned int i, const unsigned int j) const {
	return ((i+j) % 2 ? -1 : 1) * Minor(i, j);
}

/**
 * @brief Calculates the Determinant of the matrix
 * 
 * @tparam H height and width of the matrix
 * @tparam T datatype of the matrix
 * @return double the Determinant of the matrix
 */
template<unsigned H, typename T>
double Mat<H, H, T>::Determinant(void) const {
	// gaussian elimination
	double det = 1;
	double detDiv = 1;
	// use a long double as big matrices can get large with this calculation methode
	Mat<H, H, long double> GJR = (*this);	// Matrix to execute the reduction on

	unsigned int i = 0;
	for(; i < H-1; i++){
		// PN = PN * P0[0] - P0 * PN[0] => [0, 0, ...,] 0, a, b, c, ...
		// divide determinant by GJR[i][i] for every loop
		// det /= (H - (i + 1)) * GJR[i][i];
		for(unsigned int j = i+1; j < H; j++){
			detDiv *= GJR[i][i];
			GJR[j] = (GJR[j] * GJR[i][i]) - (GJR[i] * GJR[j][i]);
		}
	}
	
	for(unsigned int i = 0; i < H; i++){
		det *= GJR[i][i];
	}

	det /= detDiv;

	return det;
}

/**
 * @brief Calculates the cofactor matrix of the matrix
 * 
 * @tparam H height and width of the matrix
 * @tparam T datatype of the matrix
 * @return Mat<H, H, T> The cofactor matrix
 */
template<unsigned H, typename T>
Mat<H, H, T> Mat<H, H, T>::CoFactorMatrix(void) const {
	Mat<H, H, T> cM;
	for(unsigned int i = 0; i < H; i++){
		for(unsigned int j = 0; j < H; j++){
			cM[i][j] = this->CoFactor(i, j);
		}
	}
	return cM;
}

/**
 * @brief Calculates the Adjoint of the matrix
 * 
 * @tparam H height and width of the matrix
 * @tparam T datatype of the matrix 
 * @return Mat<H, H, T> the Adjointed matrix
 */
template<unsigned H, typename T>
Mat<H, H, T> Mat<H, H, T>::Adjoint(void) const {
	return CoFactorMatrix().Transpose();
}

/**
 * @brief Calculates the inverse of the matrix
 * 
 * @throw runtime_error when matrix is not invertible
 * @tparam H height and width of the matrix
 * @tparam T datatype of the matrix
 * @return Mat<H, H, T> The inversed matrix
 */
template<unsigned H, typename T>
Mat<H, H, T> Mat<H, H, T>::Inverse(void) const {
	double det = Determinant();
	if(det == 0){
		throw std::runtime_error("Matrix is not invertible, determinant of 0");
	}
	return Adjoint() * (1 / det);
}

/**
 * @brief Base case for recursion
 * 
 * @tparam T datatype of the matrix
 */
template<typename T>
class Mat<1, 1, T> : public MatBase<1, 1, T> {

public:
	Mat() : MatBase<1, 1, T>() {}
	Mat(const Mat<1, 1, T> &m1) : MatBase<1, 1, T>(m1) {}
	Mat(const MatBase<1, 1, T> &m1) : MatBase<1, 1, T>(m1) {}
	Mat(std::initializer_list<std::initializer_list<T>> v) : MatBase<1, 1, T>(v) {}
	~Mat() {}

	template<typename N> 		Mat<1, 1, T>& operator*=(const Mat<1, 1, N> &B) 		{return (MatBase<1, 1, T>)(*this) *= (MatBase<1, 1, N>)B;} // Matrix multiplication
	template<typename N> const	Mat<1, 1, T>  operator* (const Mat<1, 1, N> &B) const 	{return (MatBase<1, 1, T>)(*this) *  (MatBase<1, 1, N>)B;} // Matrix multiplication
	template<typename N>		Mat<1, 1, T>& operator*=(const N &r)					{return (MatBase<1, 1, T>)(*this) *= r;}					// Scalar multiplication
	template<typename N> const	Mat<1, 1, T>  operator* (const N &r) 			const 	{return (MatBase<1, 1, T>)(*this) *  r;}					// Scalar multiplication
	template<typename N> 		Mat<1, 1, T>& operator+=(const Mat<1, 1, N> &B)			{return (MatBase<1, 1, T>)(*this) += (MatBase<1, 1, N>)B;}	// Matrix addition
	template<typename N> const  Mat<1, 1, T>  operator+ (const Mat<1, 1, N> &B) const	{return (MatBase<1, 1, T>)(*this) +  (MatBase<1, 1, N>)B;}	// Matrix addition
	template<typename N> 		Mat<1, 1, T>& operator-=(const Mat<1, 1, N> &B)			{return (MatBase<1, 1, T>)(*this) -= (MatBase<1, 1, N>)B;}	// Matrix subtraction
	template<typename N> const  Mat<1, 1, T>  operator- (const Mat<1, 1, N> &B) const	{return (MatBase<1, 1, T>)(*this) -  (MatBase<1, 1, N>)B;}	// Matrix subtraction
								Mat<1, 1, T>  Transpose(void) 					const 	{return (MatBase<1, 1, T>)(*this).Transpose();}

	Mat<1, 1, T> SubMatrix(const unsigned int i, const unsigned int j) const;
	double Determinant(void) const;
};

/**
 * @brief base case for submatrix recursion
 * 
 * @tparam T the submatrix
 * @param i the row to skip ( ignored )
 * @param j the column to skip ( ignored ) 
 * @return Mat<1, 1, T> the original matrix
 */
template<typename T>
Mat<1, 1, T> Mat<1, 1, T>::SubMatrix(const unsigned int i, const unsigned int j) const {
	return *this;
}

/**
 * @brief Returns the value of [0,0] in the matrix
 * This specialisation is needed for calculating the determinant
 * @tparam T the datatype of the matrix
 * @return double the determinant value
 */
template<typename T>
double Mat<1, 1, T>::Determinant(void) const {
	return (*this)[0][0];
}

}