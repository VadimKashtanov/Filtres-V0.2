#pragma once

#include "mdl.h"

//	Controle Directe
void verifier_derivee(Mdl_t * mdl);
void plume_mdl(Mdl_t * mdl);


//	Observation comportement
void comportement(Mdl_t * mdl);

//	Ecrire derivees et derivee secondes
void derivee_et_seconde(Mdl_t * mdl, uint depart);