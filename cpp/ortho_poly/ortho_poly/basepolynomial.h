#ifndef BASE_POLYNOMIAL_H
#define BASE_POLYNOMIAL_H

#include "common.h"

//! This is the base class providing the general polynomial interfaces
class ORTHOPOLY_API BasePolynomial
{
public:
	//! constructors
	BasePolynomial() = default;
	BasePolynomial(
		const MatrixXXd& X,/*!< [in] X coordinates*/
		const MatrixXXd& Y /*!< [in] Y coordinates*/
	);
	BasePolynomial(
		const MatrixXXd& X, /*!< X coordinates*/
		const MatrixXXd& Y, /*!< Y coordinates*/
		const VectorXd& j , /*!< order vector of the polynomials*/
		const VectorXd& c   /*!< coefficients for each order*/
	);
	BasePolynomial(const BasePolynomial& bp);
	BasePolynomial& operator = (const BasePolynomial& bp);

	virtual ~BasePolynomial();

public:
	//! interface methods
	virtual BasePolynomial& normalize(const MatrixXXd& X, const MatrixXXd& Y) = 0;

protected:

protected:
	MatrixXXd m_X; /*!< Cartesian X */
	MatrixXXd m_Y; /*!< Cartesian Y*/
	VectorXd m_j;  /*!< Orders*/
	VectorXd m_c;  /*!< Coefficients*/
};

#endif // !BASE_POLYNOMIAL_H


