#ifndef LEGENDRE_H
#define LEGENDRE_H

#include "common.h"
#include "basepolynomial.h"

class ORTHOPOLY_API Legendre : public BasePolynomial
{
public:
	Legendre() = default;
	Legendre(
		const MatrixXXd& X,
		const MatrixXXd& Y
	);
	Legendre(
		const MatrixXXd& X,
		const MatrixXXd& Y,
		const VectorXd& j,
		const VectorXd& c
	);
	Legendre(const Legendre& lg);
	Legendre& operator=(const Legendre& lg);
};

#endif // !LEGENDRE_H



