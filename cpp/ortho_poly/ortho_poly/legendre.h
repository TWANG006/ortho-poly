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
	/*! Normalize the X, Y into [-1, 1] unit square*/
	virtual Legendre& normalize(const MatrixXXd& X, const MatrixXXd& Y) override;

};

#endif // !LEGENDRE_H



