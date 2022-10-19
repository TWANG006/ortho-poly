#ifndef BASE_POLYNOMIAL_H
#define BASE_POLYNOMIAL_H

#include "common.h"

//! This is the base class providing the general polynomial interfaces
class ORTHOPOLY_API BasePolynomial
{
public:
	//! default constructor
	BasePolynomial() = default;

	//! construct with the Cartesian X, Y coordinates
	BasePolynomial(
		const MatrixXXd& X,/*!< [in] X coordinates*/
		const MatrixXXd& Y /*!< [in] Y coordinates*/
	);

	//! construct with the X, Y, order and coefficients
	BasePolynomial(
		const MatrixXXd& X, /*!< X coordinates*/
		const MatrixXXd& Y, /*!< Y coordinates*/
		const VectorXd& j , /*!< order vector of the polynomials*/
		const VectorXd& c   /*!< coefficients for each order*/
	);

	//! copy constructor
	BasePolynomial(const BasePolynomial& bp);

	//! copy assignment
	BasePolynomial& operator = (const BasePolynomial& bp);

	virtual ~BasePolynomial();

protected:
	MatrixXXd m_X; /*!< Cartesian X */
	MatrixXXd m_Y; /*!< Cartesian Y*/
	VectorXd m_j;  /*!< Orders*/
	VectorXd m_c;  /*!< Coefficients*/
};

#endif // !BASE_POLYNOMIAL_H



