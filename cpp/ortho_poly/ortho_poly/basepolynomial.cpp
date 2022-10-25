#include "pch.h"
#include "basepolynomial.h"


BasePolynomial::~BasePolynomial()
{
}

BasePolynomial& BasePolynomial::operator()(const MatrixXXd& X, const MatrixXXd& Y)
{
	m_X = X;
	m_Y = Y;
	return *this;
}

BasePolynomial& BasePolynomial::operator()(const VectorXd& x)
{
	m_x = x;
	return *this;
}
