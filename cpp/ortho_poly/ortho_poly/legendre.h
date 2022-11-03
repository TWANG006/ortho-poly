#ifndef LEGENDRE_H
#define LEGENDRE_H

#include "common.h"
#include "basepolynomial.h"

//! Class for generating and fitting the unassociated Legendre polynomials
/*! 
* Details:
*           1    n   ( n )       n-k        k
*	Ln(x) = --- Sigma (   )(x - 1)    (x + 1)
*           2^n  k=0  ( k )
* 
* Prom 1D to 2D:
*	Pj(x,y) = Ln(x)Lm(y)
* 
* Recursive equation:
*	(n + 1)Ln+1(x) = (2n + 1)xLn(x) - nLn-1(x), n >= 1
*/
class ORTHOPOLY_API Legendre : public BasePolynomial
{
public:
	/*! Only enable the default construction behaviors*/
	Legendre() = default;
	Legendre(const Legendre& lg) = default;
	Legendre& operator=(const Legendre& lg) = default;
	virtual ~Legendre();

public:
	//! Generate 1d Legendre polynomials with the orders in the `orders`
	/*!
	* NOTEs:
	*   1. Duplicates are accepted in the `orders`.
	*   2. Must be proceded with either the () operator or the `normalize` to assign
	*      values to `x`; e.g., lg(x).gen_1d_p(orders), otherwise, either no polynomials
	*      are generated or the existing `x` is used.
	* \return: vector of polynomials corresponding to the `orders`
	*/
	virtual vec_v gen_1d_p(
		vec_i& orders/*!< [in] vector of orders*/
	) override;

	//! Generate 1D Legendre polynomials and their derivatives
	/*!
	* NOTEs:
	*   1. Duplicates are accepted in the `orders`.
	*   2. Must be proceded with either the () operator or the `normalize` to assign
	*      values to `x`; e.g., lg(x).gen_1d_p(orders), otherwise, either no polynomials
	*      are generated or the existing `x` is used.
	* \return: tuple of <polynomials, polynomial derivatives> corresponding to
	*          the `orders`
	*/
	virtual std::tuple<vec_v, vec_v> gen_1d_p_dp(
		vec_i& orders/*!< [in] vector of orders*/
	) override;

	//! Generate 2D Legendre polynomials and their derivatives
	/*!
	* NOTEs:
	*   1. Duplicates are eliminated by default, since they are passed in as
	*      `set<int>`, which does not allow duplicate entries. This makes sense
	*      since 1 order only corresponds to 1 coefficient.
	*   2. Must be proceded with either the () operator or the `normalize` to assign 
	*      values to `X`, `Y`; e.g., lg(X, Y).gen_2d_p(orders), otherwise, either no 
	*      polynomials are generated or the existing `x` is used.
	* \return: vector of polynomials of the `orders`.
	*/
	virtual vec_m gen_2d_p(
		const set_i& j_orders/*!< [in] set of orders*/
	) override;

	//! Generate 2D surface from the Legendre polynomials with their correspnding coefficients
	/*! 
	* Using the `jorder_coeff` map, generate the corresponding Legendre polynomials, which are then
	* summed up to generate the 2D surface.
	* \return: tuple of the 2D surface and the underlying Legendre polynomials.
	*/
	virtual std::tuple<MatrixXXd, vec_m> gen_2d_p(
		const map_id& jorder_coeff/*!< [in] map of <j, c> pairs*/
	) override;

};

#endif // !LEGENDRE_H
