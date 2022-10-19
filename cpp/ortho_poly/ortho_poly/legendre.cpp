#include "pch.h"
#include "legendre.h"

Legendre::Legendre(const MatrixXXd& X, const MatrixXXd& Y)
	: BasePolynomial(X, Y)
{
}

Legendre::Legendre(const MatrixXXd& X, const MatrixXXd& Y, const VectorXd& j, const VectorXd& c)
	: BasePolynomial(X, Y, j, c)
{
}

Legendre::Legendre(const Legendre& lg)
	: BasePolynomial(lg)
{
}

Legendre& Legendre::operator=(const Legendre& lg)
{
	BasePolynomial::operator=(lg);
	return *this;
}
