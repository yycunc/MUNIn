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

int InputObsExpName; //name of the input observed and expected PLAC-Seq contact frequency, long format
int InputParaName;   //name of the input initial nuisance parmemters
vector< vector<double> >  DataInput;
int OutputDirName;
char StrOutDir[max_length];
char StrOutCur[max_length];
int NumInput;
int NumTissue;      //number of tissues;
int NumPoint;	    //size of Hi-C contact matrix
int NumGibbs;       //number of gibbs sample
int NumTune;        //size of Metropolis-Hasting tuning
int bininitial;	    //Initial position of all the bins
int binsize;        //Size of each bin

int InputAlphaName;
int ThetaInput;
int PhiInput;
int GammaInput;
int PsiInput;
int AlphaSize;
vector<double> theta;
vector<double> phi;
vector<double> IsingGamma;
vector<double> IsingPsi;
vector<double> alpha;

vector<double> theta_orig;
vector<double> phi_orig;
vector<double> IsingGamma_orig;
vector<double> IsingPsi_orig; 

vector< vector< vector<int> > > data; //peak indicator: z=1 peak, z=-1 background
vector< vector< vector<int> > > dataMode; //best mode of peak indicator
vector< vector< vector<int> > > x;		//observed PLAC-Seq contact frequency
vector< vector< vector<double> > > e;	//expected PLAC-Seq contact frequency
vector< vector< vector<double> > > pvalue;	
vector<vector<double> > InputPara; // NumSample * 4 matrix, keep initial values for input parameters.
vector< vector<double> > pvalueSP;
vector< vector<double> > pvalueGS;
vector< vector<double> > pvalueIS;
vector< vector<double> > loglikeSP;
vector< vector<double> > loglikeGS;
vector< vector<double> > loglikeIS;
vector< vector< vector<double> > > p;

vector< vector< vector<int> > > dataModetest;
vector< vector< vector<double> > > ptest;

//define gsl random number generator
gsl_rng *rng;
int seed;

#include "MUNIn_head.h"
#include "MUNIn_toolbox.cpp"

int main(int argc, char **argv)
{
 time_t time_start, time_end;
 struct tm * timeinfo_start, * timeinfo_end;

 time(&time_start);
 timeinfo_start=localtime(&time_start);

 int i,j,k;
 char * pch=NULL;
 char str[max_length];

 for(i=1;i<argc-1;i++)
 {
  if(strcmp(argv[i],"-I")==0){InputObsExpName=i+1;i++;}
  else if(strcmp(argv[i],"-O")==0){OutputDirName=i+1;i++;}
  else if(strcmp(argv[i],"-NP")==0){NumPoint=atoi(argv[i+1]);i++;}  
  else if(strcmp(argv[i],"-NT")==0){NumTissue=atoi(argv[i+1]);i++;} 
  else if(strcmp(argv[i],"-NG")==0){NumGibbs=atoi(argv[i+1]);i++;}
  else if(strcmp(argv[i],"-Bininitial")==0){bininitial=atoi(argv[i+1]);i++;} 
  else if(strcmp(argv[i],"-Binsize")==0){binsize=atoi(argv[i+1]);i++;}
  else if(strcmp(argv[i],"-Theta")==0){ThetaInput=i+1;i++;}
  else if(strcmp(argv[i],"-Phi")==0){PhiInput=i+1;i++;}
  else if(strcmp(argv[i],"-Gamma")==0){GammaInput=i+1;i++;}
  else if(strcmp(argv[i],"-Psi")==0){PsiInput=i+1;i++;}
  else if(strcmp(argv[i],"-Alpha")==0){InputAlphaName=i+1;i++;}  
  else if(strcmp(argv[i],"-SEED")==0){seed=atoi(argv[i+1]);i++;}
 }

 strcpy( StrOutDir, argv[ OutputDirName ] );
 strcpy( StrOutCur, StrOutDir );

 FILE * s_time;
 s_time=fopen(strcat(StrOutCur, "/Record_time.txt"), "w");
 fprintf(s_time, "start time:\t%s\n", asctime(timeinfo_start));
 
 rng = gsl_rng_alloc(gsl_rng_rand48);
 //seed = time(NULL)*getpid();
 gsl_rng_set(rng, seed);
 
 //begin read in observed and expected contact frequency matrix
 for(i=0; i<NumPoint; i++)
 {
  data.push_back( vector< vector<int> >() );
  dataMode.push_back( vector< vector<int> >() );
  x.push_back( vector< vector<int> >() );
  e.push_back( vector< vector<double> >() );
  pvalue.push_back( vector< vector<double> >() );
  p.push_back( vector< vector<double> >() );

  dataModetest.push_back( vector< vector<int> >() );
  ptest.push_back( vector< vector<double> >() );
  for(j=0; j<NumPoint; j++)
  {
   data[i].push_back( vector<int>() );
   dataMode[i].push_back( vector<int>() );
   x[i].push_back( vector<int>() );
   e[i].push_back( vector<double>() );
   pvalue[i].push_back( vector<double>() );
   p[i].push_back( vector<double>() );

   dataModetest[i].push_back( vector<int>() );
   ptest[i].push_back( vector<double>() );
   for(k=0; k<NumTissue; k++)
   {
    data[i][j].push_back(0);
    dataMode[i][j].push_back(0);
    x[i][j].push_back(0);
    e[i][j].push_back(0);
    pvalue[i][j].push_back(0);
    p[i][j].push_back(0);

    dataModetest[i][j].push_back(0);
    ptest[i][j].push_back(0);
   } 
  }
 }
 
 for(i=0; i<NumPoint; i++)
 {
  pvalueSP.push_back( vector<double>() );
  pvalueGS.push_back( vector<double>() );
  pvalueIS.push_back( vector<double>() );

  loglikeSP.push_back( vector<double>() );
  loglikeGS.push_back( vector<double>() );
  loglikeIS.push_back( vector<double>() );
  for(j=0; j<NumPoint; j++)
  {
   pvalueSP[i].push_back(0);
   pvalueGS[i].push_back(0);
   pvalueIS[i].push_back(0);

   loglikeSP[i].push_back(0);
   loglikeGS[i].push_back(0);
   loglikeIS[i].push_back(0);
  }
 }

 //count the size of the input file
 i = 0;
 fstream file_count_input(argv[ InputObsExpName ],ios::in);
 while(!file_count_input.eof() ) 
 {
  file_count_input.getline(str,max_length);
  i++;
 }         
 file_count_input.close();

 NumInput=i-1; 
 printf("NumInput = %d\n", NumInput);	
 
 //initialize  DataInput;
 for(i=0; i<NumInput; i++)
 {
  DataInput.push_back( vector<double> () );
  for(j=0; j<7; j++) DataInput[i].push_back(0);
 }
 
 //size of alpha vector = 2^NumTissue
 AlphaSize = (int) pow(2, NumTissue);

 printf("NumPoint = %d, NumTissue = %d, AlphaSize = %d\n", NumPoint, NumTissue, AlphaSize);
 
 //initialize alpha vector 
 vector< vector<double> > alphaInput;
 for(i=0; i<AlphaSize; i++)
 {
  alphaInput.push_back( vector<double>() );
  for(j=0; j<NumTissue+3; j++)	alphaInput[i].push_back(0);
 }

 //load alpha
 i=0;
 j=0;
 fstream file_alpha_input(argv[ InputAlphaName ],ios::in);
 while(!file_alpha_input.eof() && i<AlphaSize)
 {
  file_alpha_input.getline(str,max_length);
  pch=strtok(str,"\t");
  while(pch!=NULL && j<NumTissue+3)
  {
   alphaInput[i][j] = atof(pch);
   pch=strtok(NULL, "\t");
   j++;
  }
   //alpha.push_back( alphaInput[i][NumTissue+2] ); 
   alpha.push_back( (alphaInput[i][NumTissue+2] * (double)NumInput/NumTissue + 1.0)/( (double)NumInput/NumTissue + 4.0) );	//define alpha vector 
   //printf("alpha = %.8lf", alpha[i]);
   j=0;
   i++;
 }
 file_alpha_input.close();

 i=0;
 j=0;
 fstream file_op_input(argv[ InputObsExpName ],ios::in);
 while(!file_op_input.eof() && i<NumInput) 
 {
  file_op_input.getline(str,max_length); 
  pch=strtok(str,"\t");
  while(pch!=NULL && j<7)
  {
   DataInput[i][j]=atof(pch);
   pch=strtok(NULL, "\t");
   j++;
  } 
   j=0;
   i++;   
 }     
 file_op_input.close();

 // fill in x and e
 for(i=0; i<NumInput; i++)
 {
  int SampleID = (int) DataInput[i][0];
  int idi = ( (int) DataInput[i][1] - bininitial)/binsize;
  int idj = ( (int) DataInput[i][2] - bininitial)/binsize;
  int obscount = (int) DataInput[i][3];
  double expcount = DataInput[i][4];
  int peak = (int) DataInput[i][5];
  double p = DataInput[i][6];

  x[idi][idj][SampleID] = obscount;
  x[idj][idi][SampleID] = obscount;
//printf("%d\n", x[idi][idj][SampleID]);
  e[idi][idj][SampleID] = expcount;
  e[idj][idi][SampleID] = expcount; 
  
  data[idi][idj][SampleID] = peak;
  data[idj][idi][SampleID] = peak;

  pvalue[idi][idj][SampleID] = p;
  pvalue[idi][idj][SampleID] = p;
 }
 
 for(i=0; i<NumPoint-1; i++)
 {
  for(j=i+1; j<NumPoint; j++)
  {
   for(k=0; k<NumTissue; k++)
   {
    dataMode[i][j][k] = data[i][j][k];
    dataMode[j][i][k] = dataMode[i][j][k];
   }
  }
 }

 for(i=0; i<NumTissue; i++)
 {
  theta.push_back(0);
  phi.push_back(0);
  IsingGamma.push_back(0);
  IsingPsi.push_back(0);

  theta_orig.push_back(0);
  phi_orig.push_back(0);
  IsingGamma_orig.push_back(0);
  IsingPsi_orig.push_back(0);
 }

 i=0;
 fstream file_theta_input(argv[ ThetaInput ],ios::in);
 while(!file_theta_input.eof() && i<NumTissue)
 {
  file_theta_input.getline(str,max_length);
  theta[i] = atof(str);
  theta_orig[i] = atof(str);
  i++;
 }
 file_theta_input.close(); 

 i=0;
 fstream file_phi_input(argv[ PhiInput ],ios::in);
 while(!file_phi_input.eof() && i<NumTissue)
 {
  file_phi_input.getline(str,max_length);
  phi[i] = atof(str);
  phi_orig[i] = atof(str);
  i++;
 }
 file_phi_input.close();

 i=0;
 fstream file_gamma_input(argv[ GammaInput ],ios::in);
 while(!file_gamma_input.eof() && i<NumTissue)
 {
  file_gamma_input.getline(str,max_length);
  IsingGamma[i] = atof(str);
  IsingGamma_orig[i] = atof(str);
  i++;
 }
 file_phi_input.close();
 
 i=0;
 fstream file_psi_input(argv[ PsiInput ],ios::in);
 while(!file_psi_input.eof() && i<NumTissue)
 {
  file_psi_input.getline(str,max_length);
  IsingPsi[i] = atof(str);
  IsingPsi_orig[i] = atof(str);
  i++;
 }
 file_psi_input.close();

 component Example;
 
// initialization
 Example.InitialPara();
 //Example.InitialZ();
 
 component ExampleMode;
 ExampleMode = Example;

 // Gibbs sample  
 Example.GibbsRefine(ExampleMode);
 
 time(&time_end);
 timeinfo_end=localtime(&time_end);

 fprintf(s_time, "end time:\t%s\n", asctime(timeinfo_end));
 fprintf(s_time, "run time:\t%d\tseconds\n", (int)(time_end-time_start));

 fclose(s_time);

 printf("END!\n");
 return 0;
}
