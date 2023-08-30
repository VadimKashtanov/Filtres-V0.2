#pragma once

//	Optimisation & Gains

#include "mdl.h"

//	Analyze seconde
float dp2(Mdl_t * mdl, uint depart, uint p0, uint p1);
float estimer_alpha(Mdl_t * mdl, uint depart, uint N);

//=somme(gains); w -= f'(x) * alpha; alpha = 1/moy(ddf(x))
float score(Mdl_t * mdl);