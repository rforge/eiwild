#include <stdio.h>
#include <Rmath.h>
#include "useful.h"


double llBeta(int Z_rc, int Z_rC, 
				double beta_rc, double beta_rC,
				int T_c, int T_C,
				double sumbeta_rc, double sumbeta_rC,
				double alpha_rc, double alpha_rC){
	double ret;
	ret= Z_rc*log(beta_rc) + Z_rC*log(beta_rC)
			+ T_c*log(sumbeta_rc) + T_C*log(sumbeta_rC)
			+ (alpha_rc-1)*log(beta_rc) + (alpha_rC-1)*log(beta_rC);
	return ret;
}

double llAlphaGamma(int P, double sumalpha,
				double alpharc,
				double sumlogbetarc,
				double la1, double la2){
	double ret;
	// Rprintf("%f \n",la1);
	ret= P*(lgammafn(sumalpha)) - P*lgammafn(alpharc) + alpharc*sumlogbetarc + (la1-1)*log(alpharc) - alpharc*la2;
	 // - P*lgammafn(alpharc)   ;
	// ll =   +  - lgammafn(alpha)) ;

	return ret;
}

double llAlphaExpo(int P, double sumalpha,
				double alpharc,
				double sumlogbetarc,
				double la){
	double ret;
	// Rprintf("%0.10f \n",la);
	ret= P*(lgammafn(sumalpha)) - P*lgammafn(alpharc) + alpharc*sumlogbetarc + log(la) - la*alpharc;
	
	return ret;
}
