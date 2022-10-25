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
	virtual std::vector<VectorXd> gen_1d_poly(vec_i& orders) override;
};

#endif // !LEGENDRE_H



