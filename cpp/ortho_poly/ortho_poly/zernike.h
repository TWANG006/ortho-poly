#ifndef ZERNIKE_H
#define ZERNIKE_H

#include "common.h"
#include "basepolynomial.h"

//! Class for generating and fitting the Zernike polynomials
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
class ORTHOPOLY_API Zernike : public BasePolynomial
{
public:
	Zernike() = default;
	Zernike(const Zernike& bp) = default;
	Zernike& operator = (const Zernike& bp) = default;
	virtual ~Zernike();

public:
	virtual vec_m gen_2d_p(const set_i& j_orders) override;
	virtual std::tuple<MatrixXXd, vec_m> gen_2d_p(const map_id& jorder_coeff) override;

protected:
	virtual std::tuple<int_t, int_t> _mn_from_j(const int& j) override;

private:
	virtual vec_v gen_1d_p(vec_i& orders) override { return vec_v(); };
	virtual std::tuple<vec_v, vec_v> gen_1d_p_dp(vec_i& orders) override { return std::tuple<vec_v, vec_v>(); };
	MatrixXXd _calculate_Rmn(
		map_iim& mn_p_map, 
		const MatrixXXd& R,
		const MatrixXXd& R0, 
		const int_t& m_order, 
		const int_t& n_order);
};

#endif // !ZERNIKE_H
