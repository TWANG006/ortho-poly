#ifndef LEGENDRE_H
#define LEGENDRE_H

#include "common.h"
#include "basepolynomial.h"

class ORTHOPOLY_API Legendre : public BasePolynomial
{
public:
	Legendre() = default;
	Legendre(const Legendre& lg) = default;
	Legendre& operator=(const Legendre& lg) = default;

	virtual ~Legendre();

public:
	/*! Normalization to the [-1, 1] unit square*/
	virtual Legendre& normalize(const MatrixXXd& X, const MatrixXXd& Y) override;
	virtual Legendre& normalize(const VectorXd& x) override;

	/*! 1D methods*/
	virtual vec_v gen_1d_p(vec_i& orders) override;
	virtual std::tuple<vec_v, vec_v> gen_1d_p_dp(vec_i& orders) override;

	/*! 2D methods*/
	virtual vec_m gen_2d_p(set_i& j_orders) override;
};

#endif // !LEGENDRE_H
