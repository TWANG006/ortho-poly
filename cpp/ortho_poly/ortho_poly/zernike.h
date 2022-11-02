#ifndef ZERNIKE_H
#define ZERNIKE_H

#include "common.h"
#include "basepolynomial.h"

class ORTHOPOLY_API Zernike : public BasePolynomial
{
public:
	Zernike() = default;
	Zernike(const Zernike& bp) = default;
	Zernike& operator = (const Zernike& bp) = default;
	virtual ~Zernike();
};

#endif // !ZERNIKE_H
