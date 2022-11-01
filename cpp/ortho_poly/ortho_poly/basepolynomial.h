#ifndef BASE_POLYNOMIAL_H
#define BASE_POLYNOMIAL_H

#include "common.h"

//! This is the base class providing the general polynomial interfaces
class ORTHOPOLY_API BasePolynomial
{
public:
	//! constructors
	BasePolynomial() = default;
	BasePolynomial(const BasePolynomial& bp) = default;
	BasePolynomial& operator = (const BasePolynomial& bp) = default;
	virtual ~BasePolynomial();

	virtual BasePolynomial& operator() (const MatrixXXd& X, const MatrixXXd& Y);
	virtual BasePolynomial& operator() (const VectorXd& x);

public:
	//! fit the surface `Z` with the `j_orders` polynomials
	BasePolynomial& fit(const MatrixXXd& Z, const set_i& j_orders);

	//! predict the fitted surface using the existing (X, Y) orders and coefficients
	MatrixXXd predict();

	//! predict the fitted surface with the existing (X, Y) and the input {orders, coeff}
	MatrixXXd predict(const map_id& order_coeff_map);

	//! fit and predict based on the input surface `Z` and `j_orders`.
	MatrixXXd fit_predict(const MatrixXXd& Z, const set_i& j_orders);

	/*! Decompose the 1D order j into 2D order [n, m] in x and y directions */
	std::tuple<int_t, int_t> _mn_from_j(const int& j);

public:
	//! Getters, for debug only. Will be deleted later.
	const MatrixXXd& GetX() const { return m_X; }
	const MatrixXXd& GetY() const { return m_Y; }
	const map_id& coeffs() const { return m_order_coeff_map; }

public:
	//! interface methods
	virtual BasePolynomial& normalize(const MatrixXXd& X, const MatrixXXd& Y) = 0;
	virtual BasePolynomial& normalize(const VectorXd& x) = 0;
	virtual vec_v gen_1d_p(vec_i& orders) = 0;
	virtual std::tuple<vec_v, vec_v> gen_1d_p_dp(vec_i& orders) = 0;
	virtual vec_m gen_2d_p(const set_i& j_orders) = 0;
	virtual std::tuple<MatrixXXd, vec_m> gen_2d_p(const map_id& jorder_coeff) = 0;


private:
	/*! Build the Ax = b linear system for polynomial fitting */
	VectorXd _build_solve_Axb(const vec_m& Ps, const MatrixXXd& Z);

protected:
	VectorXd m_x;/*!< 1d Cartesian coordinates*/
	MatrixXXd m_X;/*!< 2d Cartesian X */
	MatrixXXd m_Y;/*!< 2d Cartesian Y*/
	map_id m_order_coeff_map;/*!< ordered (order, coeff) pair*/
};

#endif // !BASE_POLYNOMIAL_H
