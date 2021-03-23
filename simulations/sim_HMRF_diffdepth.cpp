//standard libraries
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

//stantard containers
#include <string.h>
#include <deque>
#include <queue>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <iterator>
#include <stack>
//algorithm
#include <algorithm>
//stream
#include <iostream>
#include <fstream>
#include <sstream>
//limits
#include <limits>

//memory usage
#include <sys/resource.h>
#include <sys/time.h>

using namespace std;

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_cdf.h>
#include <gsl_sf_gamma.h>
#include <gsl_sf_exp.h>

#define max_length 10000
#define PI 3.1415926

int OutputDirName;	//Directory for output files
char StrOutDir[max_length];
char StrOutCur[max_length];
int NumPoint;		//size of Hi-C contact matrix
int NumTissue;    //number of sample
int NumGibbs;		//number of gibbs sample
int InputAlphaName;
int AlphaSize;
int AlphaHTSize;
vector<double> alphaHT;

double Theta;
double Phi;
double IsingGamma;
double IsingPsi;
vector< vector<int> > dataHT;	//peak status for the hiden tissue
vector< vector<int> > datamodeHT;
vector< vector< vector<int> > > data;  //peak status

//define gsl random number generator
gsl_rng *rng;
int seed;


double GetLogLikelihood_HT( vector < vector<int> > & zHT )	//zinput, peak status matrix
{
 int i, j;
 double sum = 0;
 
 for(i=0; i<NumPoint-1; i++)
 {
  for(j=i+1; j<NumPoint; j++)
  {
   int count = 0;
   if(i-1 >= 0)       count += zHT[i-1][j];
   if(i+1 <= NumPoint-1)    count += zHT[i+1][j];
   if(j-1 >= 0)       count += zHT[i][j-1];
   if(j+1 <= NumPoint-1)  count += zHT[i][j+1];  

   sum += -log( 1.0 + exp(-2 * IsingGamma * zHT[i][j] - 2 * IsingPsi * zHT[i][j] * count) );
  }
 }

 return sum;
}

void GibbsUpdate_HT( vector< vector<int> > & zHT, int i, int j )
{
 int count = 0;
 if(i-1 >= 0)       count += zHT[i-1][j];
 if(i+1 <= NumPoint-1)    count += zHT[i+1][j];
 if(j-1 >= 0)       count += zHT[i][j-1];
 if(j+1 <= NumPoint-1)  count += zHT[i][j+1];

 double valueP = 1.0/(1.0 + exp(-2 * IsingGamma - 2 * IsingPsi * count));
 
 double value = gsl_rng_uniform(rng); 
 if(value < valueP)
 {
  zHT[i][j] = 1;
  zHT[j][i] = 1;
 }
 else if(value >= valueP)
 {
  zHT[i][j] = -1;
  zHT[j][i] = -1;
 }
}

void PeakUpdate( vector< vector < vector<int> > > & zinput, vector< vector<int> > & zHT, int SampleID )
{
 int i, j;
 
 for(i=0; i<NumPoint-1; i++)
 {
  for(j=i+1; j<NumPoint; j++)
  {
   double prob = alphaHT[0] * 2;
   double value = gsl_rng_uniform(rng);
   if(value < prob)
   {
    zinput[i][j][SampleID] = zHT[i][j];
   }
   else if(value >= prob)
   {
    zinput[i][j][SampleID] = -zHT[i][j];
   }
  }
 }
}

double GetLogLL_para( vector< vector < vector<int> > > & zinput, double gamma, double psi, int SampleID )      //zinput, peak status matrix
{
 int i, j;
 double sum = 0;

 for(i=0; i<NumPoint-1; i++)
 {
  for(j=i+1; j<NumPoint; j++)
  {
   int count = 0;
   if(i-1 >= 0)       count += zinput[i-1][j][SampleID];
   if(i+1 <= NumPoint-1)    count += zinput[i+1][j][SampleID];
   if(j-1 >= 0)       count += zinput[i][j-1][SampleID];
   if(j+1 <= NumPoint-1)  count += zinput[i][j+1][SampleID];

   sum += -log( 1.0 + exp(-2 * gamma * zinput[i][j][SampleID] - 2 * psi * zinput[i][j][SampleID] * count) );
  }
 }

 return sum;
}
 

int main(int argc, char **argv)
{
 time_t time_start, time_end;
 struct tm * timeinfo_start, * timeinfo_end;

 time(&time_start);
 timeinfo_start=localtime(&time_start);

 int i, j, k;
 char * pch=NULL;
 char str[max_length];
 
 for(i=1; i<argc-1; i++)
 {
  if(strcmp(argv[i],"-O")==0){OutputDirName=i+1;i++;}
  else if(strcmp(argv[i],"-NP")==0){NumPoint=atoi(argv[i+1]);i++;}
  else if(strcmp(argv[i],"-NG")==0){NumGibbs=atoi(argv[i+1]);i++;}
  else if(strcmp(argv[i],"-NT")==0){NumTissue=atoi(argv[i+1]);i++;}
  else if(strcmp(argv[i],"-Theta")==0){Theta=atof(argv[i+1]);i++;}
  else if(strcmp(argv[i],"-Phi")==0){Phi=atof(argv[i+1]);i++;}
  else if(strcmp(argv[i],"-Gamma")==0){IsingGamma=atof(argv[i+1]);i++;}
  else if(strcmp(argv[i],"-Psi")==0){IsingPsi=atof(argv[i+1]);i++;}
  else if(strcmp(argv[i],"-Alpha")==0){InputAlphaName=i+1;i++;}
  else if(strcmp(argv[i],"-SEED")==0){seed=atoi(argv[i+1]);i++;}
 }

 strcpy( StrOutDir, argv[ OutputDirName ] );
 strcpy( StrOutCur, StrOutDir );

 FILE * s_time;
 s_time=fopen(strcat(StrOutCur, "/time.txt"), "w");
 fprintf(s_time, "start time:\t%s\n", asctime(timeinfo_start));

 rng = gsl_rng_alloc(gsl_rng_rand48);
 //seed = time(NULL)*getpid();
 gsl_rng_set(rng, seed);
 
 //size of alpha vector = 2^NumTissue
 AlphaSize = (int) pow(2, NumTissue);
 
 printf("NumPoint = %d, NumTissue = %d, AlphaSize = %d\n", NumPoint, NumTissue, AlphaSize);

 //initialize alpha vector
 int NumTissueHT = 2;
 AlphaHTSize = (int) pow(2, 2);
 vector< vector<double> > alphaInput;
 for(i=0; i<AlphaHTSize; i++)
 {
  alphaInput.push_back( vector<double>() );
  for(j=0; j<NumTissueHT+3; j++)  alphaInput[i].push_back(0);
 }

 //load alpha
 i=0;
 j=0;
 fstream file_alpha_input(argv[ InputAlphaName ],ios::in);
 while(!file_alpha_input.eof() && i<AlphaHTSize)
 {
  file_alpha_input.getline(str,max_length);
  pch=strtok(str,"\t");
  while(pch!=NULL && j<NumTissueHT+3)
  {
   alphaInput[i][j] = atof(pch);
   pch=strtok(NULL, "\t");
   j++;
  }
   alphaHT.push_back( alphaInput[i][NumTissueHT+2] );       //define alpha vector
   j=0;
   i++;
 }
 file_alpha_input.close();

 //define peak status matrix for hidden tissue
 for(i=0; i<NumPoint; i++)
 {
  dataHT.push_back( vector<int>() );
  for(j=0; j<NumPoint; j++)	dataHT[i].push_back(0);
 } 

 //initialize peak status for hidden tissue
 for(i=0; i<NumPoint-1; i++)
 {
  for(j=i+1; j<NumPoint; j++)
  {
   double value = gsl_rng_uniform(rng);
   if( value < 0.5 )
   {
    dataHT[i][j] = 1;
    dataHT[j][i] = dataHT[i][j];
   }
   else if( value >= 0.5 )
   {
    dataHT[i][j] = -1;
    dataHT[j][i] = dataHT[i][j];
   }
  }
 }

 //define peak status matrix
 for(i=0; i<NumPoint; i++)
 {
  data.push_back( vector< vector<int> >() );
  for(j=0; j<NumPoint; j++)
  {
   data[i].push_back( vector<int>() );
   for(k=0; k<NumTissue; k++) data[i][j].push_back(0);
  }
 }
 
 //initialize peak status
 for(k=0; k<NumTissue; k++)
 {
  for(i=0; i<NumPoint-1; i++)
  {
   for(j=i+1; j<NumPoint; j++)
   {
    double value = gsl_rng_uniform(rng);   
    if( value < 0.5 )
    {
     data[i][j][k] = 1;
     data[j][i][k] = data[i][j][k];  
    }
    else if( value >= 0.5 )
    {
     data[i][j][k] = -1;
     data[j][i][k] = data[i][j][k]; 
    }  
   }  
  }
 }
 
 printf("Initial Loglike = %.8lf\n", GetLogLikelihood_HT(dataHT));

 //Use Gibbs sample to update peak status
 datamodeHT = dataHT;
 double loglikeHT = GetLogLikelihood_HT(dataHT);
 double loglikemodeHT = loglikeHT;

 strcpy( StrOutCur, StrOutDir );
 FILE * sht = fopen(strcat(StrOutCur, "/Record_loglike_hiden_tissue.txt"), "w");
 for(int index=0; index<NumGibbs; index++)
 {
  for(i=0; i<NumPoint-1; i++)
  {
   for(j=i+1; j<NumPoint; j++)	GibbsUpdate_HT(dataHT, i, j);
  }

  loglikeHT = GetLogLikelihood_HT(dataHT);

  if(loglikeHT > loglikemodeHT)
  {
   datamodeHT = dataHT;
   loglikemodeHT = loglikeHT;
 //  printf("%.4lf\n", GetAlpha(data, alpha));
  }

  fprintf(sht, "%.8lf\t%.8lf\n", loglikeHT, loglikemodeHT);

  if(index%1000==0) printf("index = %d\tloglike = %.8lf\tloglikemode = %.8lf\n", index, loglikeHT, loglikemodeHT);
 }
 fclose(sht);

 //Update peak status based on hidden tissue
 for(k=0; k<NumTissue; k++)	PeakUpdate(data, dataHT, k);

 //simulate Hi-C count data observed and expected
 vector< vector< vector<int> > > x;
 vector< vector< vector<double> > > e;
 for(i=0; i<NumPoint; i++)
 {
  x.push_back(vector< vector<int> >());
  e.push_back(vector< vector<double> >());
  for(j=0; j<NumPoint; j++)
  {
   x[i].push_back(vector<int>());
   e[i].push_back(vector<double>());
   for(k=0; k<NumTissue; k++)
   {
    x[i][j].push_back(0);
    e[i][j].push_back(0);
   }
  }
 }
 
 vector<double> NBLL; //calculate complete log likelihood of the negative binomial model. use zmax as peak indicator
 for(k=0; k<NumTissue; k++)	NBLL.push_back(0); 

 //The first (NumTissue-1) tissue have equal sequencing depth
 for(i=0; i<NumPoint-1; i++)
 {
  for(j=i+1; j<NumPoint; j++)
  {
   for(k=0; k<NumTissue-1; k++)
   {
    e[i][j][k] = (double)40/(j-i); // background model depends on 1D genomic distance.
    e[j][i][k] = e[i][j][k];

    double nbmean = e[i][j][k] * exp( Theta * (data[i][j][k]+1)/2 );
    x[i][j][k] = gsl_ran_negative_binomial(rng, Phi/(nbmean + Phi), Phi);
    x[j][i][k] = x[i][j][k];

    NBLL[k] += gsl_sf_lngamma(x[i][j][k] + Phi) - gsl_sf_lngamma(Phi) + Phi * log(Phi) + x[i][j][k] * log(nbmean) - (x[i][j][k] + Phi) * log(nbmean + Phi);
   }
  }
 }

 //The last tissue has a lower sequencing depth
 k = NumTissue-1;
 for(i=0; i<NumPoint-1; i++)
 {
  for(j=i+1; j<NumPoint; j++)
  {
   e[i][j][k] = (double)20/(j-i); // background model depends on 1D genomic distance.
   e[j][i][k] = e[i][j][k];

   double nbmean = e[i][j][k] * exp( Theta * (data[i][j][k]+1)/2 );
   x[i][j][k] = gsl_ran_negative_binomial(rng, Phi/(nbmean + Phi), Phi);
   x[j][i][k] = x[i][j][k];

   NBLL[k] += gsl_sf_lngamma(x[i][j][k] + Phi) - gsl_sf_lngamma(Phi) + Phi * log(Phi) + x[i][j][k] * log(nbmean) - (x[i][j][k] + Phi) * log(nbmean + Phi);
  }
 }

 vector<double> gammamode;
 vector<double> psimode;
 vector<double> LLmode;
 for(k=0; k<NumTissue; k++)
 {
  gammamode.push_back(0);
  psimode.push_back(0);
  LLmode.push_back(0);
 }

 for(k=0; k<NumTissue; k++)
 {
  double gamma = -0.01;
  double psi = 0.01; 
  LLmode[k] = GetLogLL_para(data, gamma, psi, k);
  gammamode[k] = gamma;
  psimode[k] = psi;  

  for(i=0; i<99; i++)
  {
   gamma = gamma - 0.01;
   for(j=0; j<99; j++)
   {
    psi = psi + 0.01;
    double LL = GetLogLL_para(data, gamma, psi, k);
    
    if(LL >= LLmode[k])
    {
     gammamode[k] = gamma;
     psimode[k] =  psi;
     LLmode[k] = LL;
    }
   }
   psi = 0.01;
  }
 }

 vector<int> TmpValue;
 for(i=0; i<NumTissue; i++)     TmpValue.push_back(0);

 strcpy( StrOutCur, StrOutDir );
 FILE * snui = fopen(strcat(StrOutCur, "/Record_nui.txt"), "w");
 fprintf(snui, "NumPoint = %d\n", NumPoint);
 fprintf(snui, "NumTissue = %d\n", NumTissue);
 fprintf(snui, "NumGibbs = %d\n", NumGibbs);
 for(k=0; k<NumTissue; k++)
 {
  fprintf(snui, "Tissue = %d\n", k);
  fprintf(snui, "Max Loglikelihood = %.8lf\n", LLmode[k]);
  fprintf(snui, "Full Loglikelihood= %.8lf\n", NBLL[k]+LLmode[k]);
  fprintf(snui, "Theta = %.4lf\n", Theta);
  fprintf(snui, "Phi = %.4lf\n", Phi);
  fprintf(snui, "IsingGamma = %.4lf\n", gammamode[k]);
  fprintf(snui, "IsingPsi = %.4lf\n", psimode[k]);
 }
 fprintf(snui, "SEED = %d\n", seed);
 fclose(snui);

 strcpy( StrOutCur, StrOutDir );
 FILE * sdata = fopen(strcat(StrOutCur, "/Record_peak.txt"), "w");
 for(i=0; i<NumPoint-1; i++)
 {
  for(j=i+1; j<NumPoint; j++)
  {
   fprintf(sdata, "%d\t%d\t", i, j);
   for(k=0; k<NumTissue; k++) fprintf(sdata, "%d\t", data[i][j][k]);
   fprintf(sdata, "\n");
  }
 }
 fclose(sdata);
 
 strcpy( StrOutCur, StrOutDir );
 FILE * slong = fopen(strcat(StrOutCur, "/Record_long_format.txt"), "w");
 for(k=0; k<NumTissue; k++)
 {
  for(i=0; i<NumPoint-1; i++)
  {
   for(j=i+1; j<NumPoint; j++)  fprintf(slong, "%d\t%d\t%d\t%d\t%.8lf\t%d\n", k, i, j, x[i][j][k], e[i][j][k], data[i][j][k]);
  }
 }
 fclose(slong);

 //QC datamode
 vector<int> countmode;
 for(i=0; i<AlphaSize; i++) countmode.push_back(0); // count the frequency of each configuration

 for(i=0; i<NumPoint-1; i++)
 {
  for(j=i+1; j<NumPoint; j++)
  {
   int numpeak=0;
   for(k=0; k<NumTissue; k++)
   {
    if(data[i][j][k]==1) numpeak += (int) pow(2, k);
   }
   countmode[numpeak]++;
  }
 }

 strcpy( StrOutCur, StrOutDir );
 FILE *sprob = fopen(strcat(StrOutCur, "/Record_prob.txt"), "w");
 for(i=0; i<AlphaSize; i++)
 {
  //convert decimal to binary, from i (decimal) -> TmpValue (binary, length 4)
  int Value = i;
  for(j=NumTissue-1; j>=0; j--)
  {
   TmpValue[j] = (int) floor( Value / pow(2, j) );
   Value = Value - TmpValue[j]*(int)pow(2, j);
  }
  fprintf(sprob, "%d\t", i);
  for(j=0; j<NumTissue; j++) fprintf(sprob, "%d\t", TmpValue[j]);
  
  int CountS = 0;
  for(j=0; j<NumTissue; j++) CountS += TmpValue[j]; //count the number of tissues with peak, defined as S
  fprintf(sprob, "%d\t%.8lf\n", CountS, countmode[i]/( NumPoint*(NumPoint-1)/2.0 ));
 }
 fclose(sprob);
/*
 strcpy( StrOutCur, StrOutDir );
 FILE *sprob = fopen(strcat(StrOutCur, "/Record_prob.txt"), "w");
 for(i=0; i<(int)alpha.size(); i++) fprintf(sprob, "%d\t%d\t%.8lf\n", i, countmode[i], countmode[i]/( NumPoint*(NumPoint-1)/2.0 ));
 fclose(sprob);
*/
 time(&time_end);
 timeinfo_end=localtime(&time_end);

 fprintf(s_time, "end time:\t%s\n", asctime(timeinfo_end));
 fprintf(s_time, "run time:\t%d\tseconds\n", (int)(time_end-time_start));

 fclose(s_time);

 printf("END\n");
 return 0;
}

