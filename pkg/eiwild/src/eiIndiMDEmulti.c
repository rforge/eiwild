#include <stdio.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "useful.h"


// Einziger Unterschied zu eiIndiMDG1 ist dass überall
// statt la1 =  la1[rr+cc*r] und la2 = la2[rr+cc*r] steht

// R-werte Gross und c-variablen klein
SEXP eiIndiMDEmulti(SEXP C, SEXP R, SEXP AGGPREC,
				SEXP ALPHAS, SEXP BETAS,
				SEXP TDF, SEXP ZRC, SEXP XDF, SEXP NDF,
				SEXP LA, SEXP ALPHAVARS, SEXP BETAVARS,
				SEXP BURNIN, SEXP THINNING, SEXP SAMPLE,
				SEXP VERBOSE, SEXP RETBETA){

////////////////////////////////////////////////////////////
////////// Einlesen und Umwandeln aller Vektoren und Werte
	int c, r, aggp;
	c=INTEGER(C)[0];
	r=INTEGER(R)[0];
	aggp=INTEGER(AGGPREC)[0];
	
	// int n[aggp];
	// for(int i=0;i<aggp;i++)
		// n[i]=INTEGER(N)[i];
	
	double alphas[r*c], betas[r*c*aggp], counts[r*c];
	for(int i=0;i<r*c;i++)
		alphas[i]=REAL(ALPHAS)[i];
	for(int i=0;i<r*c*aggp;i++)
		betas[i]=REAL(BETAS)[i];
	
	int tdf[c*aggp];
	for(int i=0;i<c*aggp;i++)
		tdf[i]=INTEGER(TDF)[i];
		
	int zrc[r*c*aggp];
	for(int i=0;i<r*c*aggp;i++)
		zrc[i]=INTEGER(ZRC)[i];
		
	double xdf[r*aggp];
	for(int i=0;i<r*aggp;i++)
		xdf[i]=REAL(XDF)[i];
		
	int ndf[r*aggp];
	for(int i=0; i<r*aggp; i++)
		ndf[i]=INTEGER(NDF)[i];

		
	double la[r*c];
	for(int i=0;i<r*c;i++){
		la[i]=REAL(LA)[i];
		// Rprintf("%0.10f -- ",la[i]);
	}
	
	double alphaVars[r*c], betaVars[r*(c-1)*aggp];
	for(int i=0;i<r*c;i++)
		alphaVars[i]=REAL(ALPHAVARS)[i];
	for(int i=0;i<r*(c-1)*aggp;i++){
		betaVars[i]=REAL(BETAVARS)[i];
		// Rprintf("%f", betaVars[i]);
	}
		
	int iter, burnin, thin, sample;
	burnin=INTEGER(BURNIN)[0];
	thin  =INTEGER(THINNING)[0];
	sample=INTEGER(SAMPLE)[0];
	iter=burnin + thin*sample;
	
	int verbose;
	verbose = INTEGER(VERBOSE)[0];
	int retbeta;
	retbeta = INTEGER(RETBETA)[0];
	
	SEXP betaret, alpharet, countsret;
	PROTECT(betaret = allocMatrix(REALSXP, sample, r*c*aggp));
	PROTECT(alpharet= allocMatrix(REALSXP, sample, r*c));
	PROTECT(countsret = allocMatrix(REALSXP,sample,r*c));
	
	for(int ii=0; ii<sample;ii++){
		for(int qq=0; qq<r*c*aggp; qq++){
			REAL(betaret)[ii + qq*sample]=-1;
			// Rprintf("%f -- ", REAL(betaret)[itercurr+1 + qq*sample]);
		}
	}
	
	SEXP alpha_acc, beta_acc;
	PROTECT(alpha_acc = allocVector(REALSXP, r*c));
    PROTECT(beta_acc = allocVector(REALSXP, r*(c-1)*aggp));
	for(int qq=0; qq<r*c; qq++){
		REAL(alpha_acc)[qq] = 0;
	}
	for(int qq=0; qq<r*(c-1)*aggp; qq++){
		REAL(beta_acc)[qq] = 0;
	}
		
	GetRNGstate();
	
	//Beginn des Metropolis Algorithmus
	double betacurr, betarefcurr, betasumcurr, betarefsumcurr;
	double betanew, betarefnew, betasumnew, betarefsumnew;
	double llbetacurr, llbetanew;
	double u;
	double alphacurr,alphasumcurr;
	double alphanew,alphasumnew;
	double logbetasum[r*c];
	double llalphacurr, llalphanew;
	int itercurr =0;
	int ii,pp,rr,cc;
	double tmp; //sumcounter

for(ii=0; ii<iter;ii++){
	for(pp=0; pp<aggp; pp++){
		for(rr=0; rr<r; rr++){
			for(cc=0; cc<(c-1); cc++){
				// int pp,rr,cc;
				// pp=4;
				// rr=1;
				// cc=2;
				// Rprintf("i=%0d -- r=%0d -- c=%0d\n",pp,rr,cc);
				betacurr = betas[r*c*pp + rr + r*cc];
				betarefcurr =betas[r*c*pp + rr + r*(c-1)];
				// \sum_r^R beta_rC * X_r
				//evtl das vor die iterationschleife, da es nur aggp elemente sind

				// Rprintf("%f %f %f\n", betacurr, betarefcurr,betarefsumcurr);
				betanew = rnorm(betacurr, betaVars[r*(c-1)*pp + rr + r*cc]);
				// Rprintf("%f - ", betanew);
				betarefnew = betacurr+betarefcurr - betanew;
					// \sum_r^R beta_rC * X_r (old) - beta^rC(old)*Xr + beta^rC(new)*Xr
				// Rprintf("%f %f\n", betacurr, betarefsumcurr);
				// Rprintf("%f %f\n", betanew,betarefsumnew);
				// muss zwischen 0 und referenz+aktuell sein bspl 0.2 0.3, 0.1. 
				// 0.2 soll ersetzt werden und hat somit zwischen 0 und 0.6 "Platz"
				// 0.6 =referenz+aktuell=0.4+0.2
				if((betanew>0) & (betanew<(betarefcurr+betacurr))){
					betasumcurr=0;
					betasumnew=0;
					betarefsumcurr=0;
					for(int qq=0;qq<r;qq++){
						betasumcurr +=betas[r*c*pp + qq + r*cc] * xdf[qq + (aggp-1)*qq + pp];
						betarefsumcurr += betas[r*c*pp + qq + r*(c-1)] * xdf[qq + (aggp-1)*qq + pp];
					}
					betasumnew = betasumcurr - betacurr*xdf[rr+(aggp-1)*rr+pp] +
												betanew*xdf[rr+(aggp-1)*rr+pp];
					betarefsumnew = betarefsumcurr - betarefcurr*xdf[rr+(aggp-1)*rr+pp] + 
												betarefnew*xdf[rr+(aggp-1)*rr+pp];
					
					// Rprintf("%f %f\n", betacurr, betasumcurr);
					// Rprintf("%f %f\n", betanew,betasumnew);
					
					// Rprintf("Zrc %d -- ZrC %d\n", zrc[r*c*pp + rr + r*cc], zrc[r*c*pp + rr + r*(c-1)]); 
					// Rprintf("logbetarc %f -- logbetarC %f\n", log(betacurr), log(betarefcurr));
					// Rprintf("Tc %d -- TC %d\n", tdf[pp + (aggp-1)*cc + cc], tdf[pp+(aggp-1)*(c-1)+(c-1)]); //1695
					// Rprintf("logsumx %f -- logsumxref %f\n", log(betasumcurr),log(betarefsumcurr) );
					// Rprintf("alpha_rc %f -- alphaRC %f", alphas[rr+r*cc],alphas[rr+r*(c-1)]);
					// Rprintf("\n");
					llbetacurr = llBeta(zrc[r*c*pp + rr + r*cc], zrc[r*c*pp + rr + r*(c-1)],
										betacurr, betarefcurr,
										tdf[pp + (aggp-1)*cc + cc], tdf[pp+(aggp-1)*(c-1)+(c-1)],
										betasumcurr, betarefsumcurr,
										alphas[rr+r*cc],alphas[rr+r*(c-1)]);
					llbetanew = llBeta(zrc[r*c*pp + rr + r*cc], zrc[r*c*pp + rr + r*(c-1)],
										betanew, betarefnew,
										tdf[pp + (aggp-1)*cc + cc], tdf[pp+(aggp-1)*(c-1)+(c-1)],
										betasumnew, betarefsumnew,
										alphas[rr+r*cc],alphas[rr+r*(c-1)]);


					u=log(runif(0,1));
					if(u< llbetanew-llbetacurr){
						betas[r*c*pp + rr + r*cc] = betanew;
						betas[r*c*pp + rr + r*(c-1)] = betarefnew;
						// Rprintf("%f %f - ",betanew,betarefnew);
						REAL(beta_acc)[r*(c-1)*pp + rr + r*cc] += 1;
					}	
					// Rprintf("%f %f - ", betas[r*c*pp + rr + r*cc],  betas[r*c*pp + rr + r*(c-1)]);
				}
			}
		}
	}  
 	
	// sum_1:P log(\beta_i^rc)
	
	for(rr=0; rr<r; rr++){
		for(cc=0; cc<c; cc++){
			tmp=0;
			for(int qq=0; qq<aggp; qq++){
				tmp += log(betas[r*c*qq + rr + r*cc]);
			}
			logbetasum[rr+cc*r] =  tmp;
			// Rprintf("%f ", logbetasum[rr+cc*r]);
			// Rprintf("\n");
		}
	}
	
	for(rr=0; rr<r; rr++){
		for(cc=0; cc<c; cc++){
			alphacurr = alphas[rr+cc*r];
			alphanew = rnorm(alphacurr,alphaVars[rr+cc*r]);
			// Rprintf("%f - ", alphanew);
			
			alphasumcurr = 0;
			for(int qq=0; qq<c; qq++){ // \sum_c:C \alpha^rc)
				alphasumcurr += alphas[rr+qq*r];
			}
			alphasumnew = alphasumcurr - alphacurr + alphanew;
			
			if(alphanew>0){
				// Rprintf("%d %d %0.10f \n ",rr+1,cc+1,la[rr+cc*r]);
				llalphacurr = llAlphaExpo(aggp, alphasumcurr, alphacurr, logbetasum[rr+cc*r], la[rr+cc*r]);
				llalphanew  = llAlphaExpo(aggp, alphasumnew, alphanew, logbetasum[rr+cc*r], la[rr+cc*r]);
				// Rprintf("%f\n",llalphanew);
				u=log(runif(0,1));
				// Rprintf("%f -- ",alphacurr);
				if(u< llalphanew-llalphacurr){
					alphas[rr+cc*r] = alphanew;
					REAL(alpha_acc)[rr+cc*r] +=1;
				}		
				// Rprintf("%f - ",alphas[rr+cc*r]);
			}
		}
		// Rprintf("\n");
	}
	
	// Rprintf("%d ", ii>=burnin && ((ii+1) % thin)==0);
	// Rprintf("%d - ", (ii+1) % thin);
	if(ii==0 && ii==burnin){
		Rprintf("no Burnin\n");
	}
	if(ii==0 && ii!=burnin){
		Rprintf("Burnin start\n");
	}
	if((ii+1)==(burnin/2) && burnin>1000){
		Rprintf("Burnin half time\n");
	}
	if((ii+1)==(burnin)){
		Rprintf("Burnin finished\n");
	}
	
	if(ii>=burnin && (ii % thin)==0){
		
		// Berechnung der cellCounts
		for(int rrr=0; rrr<r; rrr++){
			for(int ccc=0; ccc<c; ccc++){
				tmp=0;
				for(int qq=0; qq<aggp; qq++){
					tmp += betas[rrr + r*ccc + r*c*qq]*ndf[qq+aggp*rrr];
					// Rprintf("%f x %d--", betas[rrr + r*ccc + r*c*qq],ndf[qq+aggp*rrr]);
				}
				counts[rrr + r*ccc] = tmp;
			}
		}
		
		
		if((itercurr+1)%verbose == 0)
			Rprintf("%2d von %d\n", itercurr+1,sample);
		//		zeilen sind iterationen
		for(int qq=0; qq<r*c*aggp; qq++){
			REAL(betaret)[itercurr + qq*sample]=betas[qq];
			// Rprintf("%f -- ", REAL(betaret)[itercurr+1 + qq*sample]);
		}
		for(int qq=0; qq<r*c; qq++){
			REAL(alpharet)[itercurr+qq*sample]=alphas[qq];
			REAL(countsret)[itercurr+qq*sample]=counts[qq];
			// Rprintf("%f ",alphas[qq]);
		}
		itercurr++;
	} 
	 
	R_CheckUserInterrupt();
	}
	
	for(int qq=0; qq < r*c; qq++){
	    REAL(alpha_acc)[qq] = REAL(alpha_acc)[qq]/iter;
	}
	for(int qq=0; qq < r*(c-1)*aggp; qq++){
	    REAL(beta_acc)[qq] = REAL(beta_acc)[qq]/iter;
	}
	
	// Liste erstellen
	SEXP ret;
	if(retbeta==0){
		int listlength=4;
		PROTECT(ret = allocVector(VECSXP, listlength));
		SET_VECTOR_ELT(ret,0,alpharet);
		SET_VECTOR_ELT(ret,1,countsret);
		SET_VECTOR_ELT(ret,2,alpha_acc);
		SET_VECTOR_ELT(ret,3,beta_acc);
	
		// Listennamen geben
		SEXP ret_names;
		PROTECT(ret_names=allocVector(STRSXP,listlength));
		char *retnames[4] ={"alphaDraws","cellCounts","alphaAcc","betaAcc"};
		for(int i=0;i<listlength;i++){
			SET_STRING_ELT(ret_names,i,mkChar(retnames[i]));
		}
		setAttrib(ret, R_NamesSymbol,ret_names);	
	} else {
		int listlength=5;
		PROTECT(ret = allocVector(VECSXP, listlength));
		SET_VECTOR_ELT(ret,0,betaret);
		SET_VECTOR_ELT(ret,1,alpharet);
		SET_VECTOR_ELT(ret,2,countsret);
		SET_VECTOR_ELT(ret,3,alpha_acc);
		SET_VECTOR_ELT(ret,4,beta_acc);
		
		// Listennamen geben
		SEXP ret_names;
		PROTECT(ret_names=allocVector(STRSXP,listlength));
		char *retnames[5] ={"betaDraws","alphaDraws","cellCounts","alphaAcc","betaAcc"};
		for(int i=0;i<listlength;i++){
			SET_STRING_ELT(ret_names,i,mkChar(retnames[i]));
		}
		setAttrib(ret, R_NamesSymbol,ret_names);
	}

	PutRNGstate();
	UNPROTECT(7);
	return ret;
	// return(R_NilValue);
}













