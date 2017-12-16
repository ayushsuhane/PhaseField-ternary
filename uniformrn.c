#include"setup_fn.h"
#include"Global_params.h"
#include <gsl/gsl_rng.h>
void UniformRNGen(){
//extern int randomseed;
const gsl_rng_type * URNG;
extern gsl_rng *URN;
gsl_rng_env_setup();

	URNG = gsl_rng_default;
	URN = gsl_rng_alloc (URNG);
//	gsl_rng_set(URN,randomseed);
//	gsl_rng_set(URN);
}
