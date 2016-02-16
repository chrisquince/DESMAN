/* C functions for running SampleTau from Cython*/

/*System includes*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

/*GSL includes*/
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

/*User includes*/
//#include "c_sample_tau.h"

static gsl_rng *ptGSLRNG;

void c_initRNG()
{
    //const gsl_rng_type * T;
    
    gsl_rng_env_setup();

    //T = gsl_rng_default;
    ptGSLRNG = gsl_rng_alloc (gsl_rng_mt19937);
}

void c_setRNG(unsigned long int seed)
{
    printf("GSL RNG initialise %lu\n",seed);
    gsl_rng_set (ptGSLRNG, seed);
}

void c_freeRNG()
{
    gsl_rng_free (ptGSLRNG);
}


void normaliseLog4(double *adLogProb)
{
    double dMax = adLogProb[0], dSum = 0.0;
    int b = 0;
    
    for(b = 1; b < 4; b++){
        if(adLogProb[b] > dMax){
            dMax = adLogProb[b];
        }
    }
    
    for(b = 0; b < 4; b++){
        adLogProb[b] = adLogProb[b] - dMax;
        dSum += exp(adLogProb[b]);
    }
    
    for(b = 0; b < 4; b++){
        adLogProb[b] = exp(adLogProb[b])/dSum;
        //printf("b=%d,p=%f\n",b,adLogProb[b]);
    }
       
    return;
}

int sample4(double *adProb, double dU){
    double adCProb[4];
    
    adCProb[0] = adProb[0];
    adCProb[1] = adProb[1] + adCProb[0];
    adCProb[2] = adProb[2] + adCProb[1];
    adCProb[3] = 1.0;
    
    if(dU < adCProb[0]){
        return 0;
    }
    else if(dU < adCProb[1])
        return 1;
    else if(dU < adCProb[2]){
        return 2;
    }
    else{
        return 3;
    }
}



int c_sample_tau (long *anTau, double* adPi, double *adEta, long* anVariants, int nV, int nG, int nS)
{
    int a = 0, b = 0;
    int g = 0, h = 0, v = 0, s = 0, t = 0;
    int nchange = 0;
    double adPSB[nS][4];
    double adPSBStore[nS][4];
    double dLogProb[4];
    int** anTauIndex = NULL;
    double u = 0.0;

    
    anTauIndex = (int **) malloc(nV*sizeof(int*));
    if(!anTauIndex)
        goto memoryError;
        
    for(v = 0; v < nV; v++){
        anTauIndex[v] = (int *) malloc(nG*sizeof(int));
        if(!anTauIndex[v])
            goto memoryError;
        for(g = 0; g < nG; g++){
            for(b = 0; b < 4; b++){
                int vIndex = v*4*nG + 4*g + b;
                
                if(anTau[vIndex] == 1){
                    anTauIndex[v][g] = b;
                    break;
                }
            }
            //printf("%d,%d,%d\n",v,g,anTauIndex[v][g]);
        }
        
    }
    
    //loop V positions
    for(v = 0; v < nV; v++){
        
        //loop G strains
        for(g = 0; g < nG; g++){
    
            //calc contribution from all other strains
            for(s = 0; s < nS; s++){
                for(b = 0; b < 4; b++){
                    adPSBStore[s][b] = 0.0;

                    for(h = 0; h < nG; h++){
                        if(h != g){
                            int pIndex = s*nG + h;
                            int eIndex = anTauIndex[v][h]*4 + b;
                            
                            adPSBStore[s][b] += adEta[eIndex]*adPi[pIndex]; 
                        }
                    }                      
                    
                }
            }

            for(a = 0; a < 4; a++){
                for(s = 0; s < nS; s++){
                    for(b = 0; b < 4; b++){                    
                        adPSB[s][b] = adPSBStore[s][b];
                        
                        adPSB[s][b] += adEta[a*4 + b]*adPi[s*nG + g];
                    }
                }
                dLogProb[a] = 0.0;
                for(s = 0; s < nS; s++){
                    for(b = 0; b < 4; b++){
                        int vIndex = v*4*nS + 4*s + b;
                        double temp = ((float) anVariants[vIndex])*log(adPSB[s][b]);
                        //printf("s=%d,b=%d,cont=%f\n",s,b,temp);
                        dLogProb[a] += temp;
                    }
                }
            }            
  //          printf("u=%f,p1=%f,p2=%f,p3=%f,p4 =%f\n",u,dLogProb[0],dLogProb[1],dLogProb[2],dLogProb[3]);
            
            normaliseLog4(dLogProb);
            
            u = gsl_rng_uniform (ptGSLRNG);
       //     printf("u=%f,p1=%f,p2=%f,p3=%f,p4 =%f\n",u,dLogProb[0],dLogProb[1],dLogProb[2],dLogProb[3]);
            t = sample4(dLogProb, u);
         //   printf("v=%d,g=%d,tnew=%d,told=%d\n",v,g,t,anTauIndex[v][g]);
            if(t != anTauIndex[v][g]){
                int vIndex = v*4*nG + 4*g;
                anTau[vIndex + anTauIndex[v][g]] = 0; 
                anTau[vIndex + t] = 1;
                
                nchange++;
                anTauIndex[v][g] = t;
            }
        } //finish sampling strain g
    
    }//finish sampling position v
    
    
    //Free up tau indices
    for(v = 0; v < nV; v++){
        free(anTauIndex[v]);
    }
    free(anTauIndex);
 

    return nchange;

    memoryError:
    fprintf(stderr, "Failed allocating memory in c_sample_tau\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
}

