#include "pch.h"
#include "basepolynomial.h"

BasePolynomial::BasePolynomial(const MatrixXXd& X, const MatrixXXd& Y)
	: m_X(X)
	, m_Y(Y)
	, m_j(VectorXd())
	, m_c(VectorXd())
{
}

BasePolynomial::BasePolynomial(const MatrixXXd& X, const MatrixXXd& Y, const VectorXd& j, const VectorXd& c)
	: m_X(X)
	, m_Y(Y)
	, m_j(j)
	, m_c(c)
{
}

BasePolynomial::BasePolynomial(const BasePolynomial& bp)
	: m_X(bp.m_X)
	, m_Y(bp.m_Y)
	, m_j(bp.m_j)
	, m_c(bp.m_c)
{
}

BasePolynomial& BasePolynomial::operator=(const BasePolynomial& bp)
{
	m_X = bp.m_X;
	m_Y = bp.m_Y;
	m_j = bp.m_j;
	m_c = bp.m_c;

	return *this;
}

BasePolynomial::~BasePolynomial()
{
}
