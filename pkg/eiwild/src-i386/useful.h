#ifndef USEFUL_H_   /* Include guard */
#define USEFUL_H_

double llBeta(int Z_rc, int Z_rC, 
				double beta_rc, double beta_rC,
				int T_c, int T_C,
				double sumbeta_rc, double sumbeta_rC,
				double alpha_rc, double alpha_rC);
				
double llAlphaGamma(int P, double sumalpha,
				double alpharc,
				double sumlogbetarc,
				double la1, double la2);
				
double llAlphaExpo(int P, double sumalpha,
				double alpharc,
				double sumlogbetarc,
				double la);

#endif