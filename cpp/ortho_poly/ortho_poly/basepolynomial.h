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
	//! common methods
	BasePolynomial& fit(const MatrixXXd& Z, const set_i& j_orders);

public:
	//! Getters, for debug only. Will be deleted later.
	const MatrixXXd& GetX() const { return m_X; }
	const MatrixXXd& GetY() const { return m_Y; }

public:
	//! interface methods
	virtual BasePolynomial& normalize(const MatrixXXd& X, const MatrixXXd& Y) = 0;
	virtual BasePolynomial& normalize(const VectorXd& x) = 0;
	virtual vec_v gen_1d_p(vec_i& orders) = 0;
	virtual std::tuple<vec_v, vec_v> gen_1d_p_dp(vec_i& orders) = 0;
	virtual vec_m gen_2d_p(const set_i& j_orders) = 0;
	virtual std::tuple<MatrixXXd, vec_m> gen_2d_p(const map_id& jorder_coeff) = 0;

private:
	VectorXd _build_solve_Axb(const vec_m& Ps, const MatrixXXd& Z);

protected:
	VectorXd m_x;/*!< 1d Cartesian coordinates*/
	MatrixXXd m_X;/*!< 2d Cartesian X */
	MatrixXXd m_Y;/*!< 2d Cartesian Y*/
	map_id m_order_coeff_map;/*!< ordered (order, coeff) pair*/
};

#endif // !BASE_POLYNOMIAL_H
