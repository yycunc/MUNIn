void component::InitialPara()
{
 int i;

 TuningPara Example;

 for(i=0; i<NumTissue; i++)
 {  
  CountTheta.push_back(0);
  CountPhi.push_back(0);
  CountGamma.push_back(0);
  CountPsi.push_back(0);  
  
  TuneTheta.push_back(Example);
  TunePhi.push_back(Example); 
  TuneGamma.push_back(Example);
  TunePsi.push_back(Example);

 }
  
 loglikefull=0;
}


double component::CalcMedian(vector<double> & scores)
{
  size_t size = scores.size();

  if (size == 0)
  {
    return 0;  // Undefined, really.
  }
  else if (size == 1)
  {
    return scores[0];
  }
  else
  {
    sort(scores.begin(), scores.end());
    if (size % 2 == 0)
    {
      return (scores[size / 2 - 1] + scores[size / 2]) / 2;
    }
    else 
    {
      return scores[size / 2];
    }
  }
}

double component::GetNB( vector< vector < vector<int> > > & zinput, double ValueTheta, double ValuePhi, int SampleID )
{
 double sum = 0;
 int i, j; 
 
 for(i=0; i<NumPoint-1; i++)
 {
  for(j=i+1; j<NumPoint; j++)
  {
   double nbmean = e[i][j][SampleID] * exp( ValueTheta * (zinput[i][j][SampleID]+1)/2 );
   sum += gsl_sf_lngamma(x[i][j][SampleID]+ValuePhi) - gsl_sf_lngamma(ValuePhi) + ValuePhi*log(ValuePhi) + x[i][j][SampleID]*log(nbmean) - (x[i][j][SampleID]+ValuePhi)*log(nbmean+ValuePhi);
  }
 }
 
 return sum;
}

double component::GetLogLikelihood( vector< vector < vector<int> > > & zinput, vector<double> & a, double ValueGamma, double ValuePsi, int SampleID )
{
 double sum = 0;
 int i, j;

 for(i=0; i<NumPoint-1; i++)
 {
  for(j=i+1; j<NumPoint; j++)
  {
   int count = 0;
   if(i-1 >= 0)				count += zinput[i-1][j][SampleID];
   if(i+1 <= NumPoint-1)   	count += zinput[i+1][j][SampleID];
   if(j-1 >= 0)				count += zinput[i][j-1][SampleID];
   if(j+1 <= NumPoint-1)	count += zinput[i][j+1][SampleID];  

   vector<int> peak; // index vector, the k th element is peak
   vector<int> back; // index vector, the k th element is background

   peak = zinput[i][j];
   back = zinput[i][j];

   peak[SampleID] = 1;
   back[SampleID] = -1;

   int id, CountPeak=0, CountBack=0;
   for(id=0; id<NumTissue; id++)
   {
    if(peak[id] == 1)     CountPeak += (int) pow(2, id);
    if(back[id] == 1)     CountBack += (int) pow(2, id);
   }
   sum += -log( 1.0 + exp(-2 * ValueGamma * zinput[i][j][SampleID] - 2 * ValuePsi * zinput[i][j][SampleID] * count) * pow(a[CountBack]/a[CountPeak], zinput[i][j][SampleID]) );    
  }
 }
  
 return sum;
}
/*
void component::GibbsUpdate(vector< vector< vector<int> > > & zinput, vector<double> & a, int i, int j, int SampleID)
{
 int count = 0;
 if(i-1 >= 0)                         count += zinput[i-1][j][SampleID];
 if(i+1 <= NumPoint-1)        count += zinput[i+1][j][SampleID];
 if(j-1 >= 0)                         count += zinput[i][j-1][SampleID];
 if(j+1 <= NumPoint-1)        count += zinput[i][j+1][SampleID];
 
 vector<int> peak;	//index vector, the k th element is peak
 vector<int> back;	//index vector, the k th element is background

 peak = zinput[i][j];
 back = zinput[i][j];
 
 peak[SampleID] = 1;
 back[SampleID] = -1;

 int id, CountPeak = 0, CountBack = 0;
 for(id=0; id<NumTissue; id++)
 {
  if(peak[id] == 1)     CountPeak += (int) pow(2, id);
  if(back[id] == 1)     CountBack += (int) pow(2, id);
 }

 //peak log likelihood 
 double valueP = x[i][j][SampleID] * theta[SampleID] - (x[i][j][SampleID] + phi[SampleID]) * log(e[i][j][SampleID] * exp(theta[SampleID]) + phi[SampleID]) - log(1.0 + exp(-2 * IsingGamma[SampleID] - 2 * IsingPsi[SampleID] * count) * a[CountBack]/a[CountPeak]);
	           
 //background log likelihood
 double valueQ = -(x[i][j][SampleID] + phi[SampleID]) * log(e[i][j][SampleID] + phi[SampleID]) - log(1.0 + exp(2 * IsingGamma[SampleID] + 2 * IsingPsi[SampleID] * count) * alpha[CountPeak]/alpha[CountBack]);

//if(i*binsize+bininitial==24140000 && j*binsize+bininitial==24620000) printf("%d\t%d\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\t%.8lf\t%.8lf\t%.4lf\t%.8lf\t%.8lf\t%.4lf\n", SampleID, zinput[i][j][SampleID], count, x[i][j][SampleID], e[i][j][SampleID], a[CountBack], a[CountPeak], valueP, valueQ, p, x[i][j][SampleID] * theta[SampleID] - (x[i][j][SampleID] + phi[SampleID]) * log(e[i][j][SampleID] * exp(theta[SampleID]) + phi[SampleID]) - log(1.0 + exp(-2 * IsingGamma[SampleID] - 2 * IsingPsi[SampleID] * count)), -(x[i][j][SampleID] + phi[SampleID]) * log(e[i][j][SampleID] + phi[SampleID]) - log(1.0 + exp(2 * IsingGamma[SampleID] + 2 * IsingPsi[SampleID] * count)), );
 
 double p = 0.0;
 if(valueP >= valueQ) p = 1.0/(1.0 + exp(valueQ - valueP));
 else if(valueP < valueQ) p = exp(valueP - valueQ)/( 1.0 + exp(valueP - valueQ) );
 
//double valueP_no = x[i][j][SampleID] * theta[SampleID] - (x[i][j][SampleID] + phi[SampleID]) * log(e[i][j][SampleID] * exp(theta[SampleID]) + phi[SampleID]) - log(1.0 + exp(-2 * IsingGamma[SampleID] - 2 * IsingPsi[SampleID] * count));

//double valueQ_no = -(x[i][j][SampleID] + phi[SampleID]) * log(e[i][j][SampleID] + phi[SampleID]) - log(1.0 + exp(2 * IsingGamma[SampleID] + 2 * IsingPsi[SampleID] * count));

//double p_no = 0.0;
//if(valueP_no >= valueQ_no) p_no = 1.0/(1.0 + exp(valueQ_no - valueP_no));
//else if(valueP_no < valueQ_no) p_no = exp(valueP_no - valueQ_no)/( 1.0 + exp(valueP_no - valueQ_no) );

//if(i*binsize+bininitial==24140000 && j*binsize+bininitial==24620000) printf("%d\t%d\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\t%.8lf\t%.8lf\t%.4lf\t%.8lf\t%.8lf\t%.4lf\n", SampleID, zinput[i][j][SampleID], count, x[i][j][SampleID], e[i][j][SampleID], a[CountBack], a[CountPeak], valueP, valueQ, p, valueP_no, valueQ_no, p_no);

 double value = gsl_rng_uniform(rng);  
if(i*binsize+bininitial==198785000 && j*binsize+bininitial==199555000) printf("%d\t%d\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\n", SampleID, zinput[i][j][SampleID], count, x[i][j][SampleID], e[i][j][SampleID], p, value); 
 if(value < p)
 {
  zinput[i][j][SampleID] = 1;
  zinput[j][i][SampleID] = zinput[i][j][SampleID];	
 }
 else if(value >= p)
 {
  zinput[i][j][SampleID] = -1;
  zinput[j][i][SampleID] = zinput[i][j][SampleID];	
 }
 
}
*/

void component::GibbsUpdate(vector< vector< vector<int> > > & zinput, vector<double> & a, int i, int j)
{
 double B11 = 0.0;
 double B01 = 0.0;
 double B10 = 0.0;
 double B00 = 0.0;

 for(int SampleID=0; SampleID<NumTissue; SampleID++)
 {
  int count = 0;
  if(i-1 >= 0)                         count += zinput[i-1][j][SampleID];
  if(i+1 <= NumPoint-1)        count += zinput[i+1][j][SampleID];
  if(j-1 >= 0)                         count += zinput[i][j-1][SampleID];
  if(j+1 <= NumPoint-1)        count += zinput[i][j+1][SampleID];

  vector<int> peak;      //index vector, the k th element is peak
  vector<int> back;      //index vector, the k th element is background

  peak = zinput[i][j];
  back = zinput[i][j];

  peak[SampleID] = 1;
  back[SampleID] = -1;

  int id, CountPeak = 0, CountBack = 0;
  for(id=0; id<NumTissue; id++)
  {
   if(peak[id] == 1)     CountPeak += (int) pow(2, id);
   if(back[id] == 1)     CountBack += (int) pow(2, id);
  }

  //peak log likelihood
  double valueP = x[i][j][SampleID] * theta[SampleID] - (x[i][j][SampleID] + phi[SampleID]) * log(e[i][j][SampleID] * exp(theta[SampleID]) + phi[SampleID]) - log(1.0 + exp(-2 * IsingGamma[SampleID] - 2 * IsingPsi[SampleID] * count) * a[CountBack]/a[CountPeak]);

  //background log likelihood
  double valueQ = -(x[i][j][SampleID] + phi[SampleID]) * log(e[i][j][SampleID] + phi[SampleID]) - log(1.0 + exp(2 * IsingGamma[SampleID] + 2 * IsingPsi[SampleID] * count) * alpha[CountPeak]/alpha[CountBack]);

  if(SampleID == 0)
  {
   B11 += valueP;
   B10 += valueP;
   B01 += valueQ;
   B00 += valueQ;
  } else if(SampleID == 1)
  {
    B11 += valueP;
    B10 += valueQ;
    B01 += valueP;
    B00 += valueQ;
   }
  }

  // pvalueShareP
  double pvalue11 = 1.0/(1.0 + exp(B10-B11) + exp(B01-B11) + exp(B00-B11));
  // pvalueFstS
  double pvalue10 = 1.0/(1.0 + exp(B11-B10) + exp(B01-B10) + exp(B00-B10));
  // pvalueSecS
  double pvalue01 = 1.0/(1.0 + exp(B11-B01) + exp(B10-B01) + exp(B00-B01));
  double pvalue00 = 1.0 - pvalue11 - pvalue10 - pvalue01;
  
  double bin1 = pvalue11;
  double bin2 = pvalue11 + pvalue10;
  double bin3 = pvalue11 + pvalue10 + pvalue01;
  double value = gsl_rng_uniform(rng);
  //printf("%.4lf\n", value);
 
  if(value <= bin1){
     zinput[i][j][0] = 1;
     zinput[i][j][1] = 1;
  } else if(value > bin1 && value <= bin2){
     zinput[i][j][0] = 1;
     zinput[i][j][1] = -1;
  } else if(value > bin2 && value <= bin3){
     zinput[i][j][0] = -1;
     zinput[i][j][1] = 1;
  } else if(value > bin3){
     zinput[i][j][0] = -1;
     zinput[i][j][1] = -1;
  }
  zinput[j][i][0] = zinput[i][j][0];
  zinput[j][i][1] = zinput[i][j][1];
//  printf("%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%d\t%d\n", pvalue11, pvalue10, pvalue01, pvalue00, value, zinput[i][j][0], zinput[i][j][1]);

}

double component::UpdateTheta(int SampleID)
{
  double Theta = theta[SampleID];
  double ValueLargest = GetNB(data, Theta, phi[SampleID], SampleID);
  double range = theta_orig[SampleID]*0.2;
  double stepsize = range/20;
  int i = 1;

  while(i<=20)
  {
    double ThetaNew1 = theta_orig[SampleID] + i * stepsize;
    double ValueNew1 = GetNB(data, ThetaNew1, phi[SampleID], SampleID);
    if(ValueNew1 > ValueLargest)
    {
      Theta = ThetaNew1;
      ValueLargest = ValueNew1;
    }
    double ThetaNew2 = theta_orig[SampleID] - i * stepsize;
    double ValueNew2 = GetNB(data, ThetaNew2, phi[SampleID], SampleID);
    if(ValueNew2 > ValueLargest)
    {
      Theta = ThetaNew2;
      ValueLargest = ValueNew2;
    }
    i++;
  }

  return Theta;
}

double component::UpdatePhi(int SampleID)
{
 double Phi = phi[SampleID];
 double ValueLargest = GetNB(data, theta[SampleID], Phi, SampleID);
 double range = phi_orig[SampleID]*0.2;
 double stepsize = range/20;
 int i = 1;

 while(i<=20)
 {
  double PhiNew1 = phi_orig[SampleID] + i * stepsize;
  double ValueNew1 = GetNB(data, theta[SampleID], PhiNew1, SampleID);
  if(ValueNew1 > ValueLargest)
  {
   Phi = PhiNew1;
   ValueLargest = ValueNew1;
  }
  double PhiNew2 = phi_orig[SampleID] - i * stepsize;
  double ValueNew2 = GetNB(data, theta[SampleID], PhiNew2, SampleID);
  if(ValueNew2 > ValueLargest)
  {
   Phi = PhiNew2;
   ValueLargest = ValueNew2;
  }
  i++;
 }

 return Phi;
}

double component::UpdateGamma(int SampleID)
{
 double Gamma = IsingGamma[SampleID];
 double ValueLargest = GetLogLikelihood(data, alpha, Gamma, IsingPsi[SampleID], SampleID);
 double range = -IsingGamma_orig[SampleID]*0.2;
 double stepsize = range/20;
 int i = 1;

 while(i<=20)
 {
  double GammaNew1 = IsingGamma_orig[SampleID] + i * stepsize;
  if(GammaNew1 > 0)
  {
    GammaNew1 = -GammaNew1;
  }
  double ValueNew1 = GetLogLikelihood(data, alpha, GammaNew1, IsingPsi[SampleID], SampleID);
  if(ValueNew1 > ValueLargest)
  {
   Gamma = GammaNew1;
   ValueLargest = ValueNew1;
  }
  double GammaNew2 = IsingGamma_orig[SampleID] - i * stepsize;
  double ValueNew2 = GetLogLikelihood(data, alpha, GammaNew2, IsingPsi[SampleID], SampleID);
  if(ValueNew2 > ValueLargest)
  {
   Gamma = GammaNew2;
   ValueLargest = ValueNew2;
  }
  i++;
 }

 return Gamma;
}

double component::UpdatePsi(int SampleID)
{
 double Psi = IsingPsi[SampleID];
 double ValueLargest = GetLogLikelihood(data, alpha, IsingGamma[SampleID], Psi, SampleID);
 double range = IsingPsi_orig[SampleID]*0.2;
 double stepsize = range/20;
 int i = 1;

 while(i<=20)
 {
  double PsiNew1 = IsingPsi_orig[SampleID] + i * stepsize;
  double ValueNew1 = GetLogLikelihood(data, alpha, IsingGamma[SampleID], PsiNew1, SampleID);
  if(ValueNew1 > ValueLargest)
  {
   Psi = PsiNew1;
   ValueLargest = ValueNew1;
  }
  double PsiNew2 = IsingPsi_orig[SampleID] - i * stepsize;
  if(PsiNew2 < 0)
  {
    PsiNew2 = -PsiNew2;
  }
  double ValueNew2 = GetLogLikelihood(data, alpha, IsingGamma[SampleID], PsiNew2, SampleID);
  if(ValueNew2 > ValueLargest)
  {
   Psi = PsiNew2;
   ValueLargest = ValueNew2;
  }
  i++;
 }

 return Psi;
}

void SummaryStat(const vector< vector< vector<int> > > & xx, vector< vector<double> > & vv)
{
 vv.clear();
 int i, j, SampleID;

 for(i=0; i<9; i++)
 {
  vv.push_back( vector<double>() );
  for(SampleID=0; SampleID<NumTissue; SampleID++)	vv[i].push_back(0);
 }

 for(SampleID=0; SampleID<NumTissue; SampleID++)
 {
  vector<double> yy;
  for(i=0; i<NumPoint-1; i++)
  {
   for(j=i+1; j<NumPoint; j++) 
   {
    yy.push_back( xx[i][j][SampleID]+0.0 );
   }
  }
  
  vector<double> yycopy (yy);
  sort(yycopy.begin(), yycopy.end());

  int size = yy.size();
  vv[0][SampleID] = yycopy[0];                //min
  vv[1][SampleID] = yycopy[(int)(size*0.1)];  //10% percentile
  vv[2][SampleID] = yycopy[(int)(size*0.25)]; //25% percentile
  vv[3][SampleID] = yycopy[(int)(size*0.5)];  //50% percentile
  vv[4][SampleID] = yycopy[(int)(size*0.75)];   //75% percentile
  vv[5][SampleID] = yycopy[(int)(size*0.9)];    //90% percentile
  vv[6][SampleID] = yycopy[size-1];             //max

  double sum=0, sum2=0;
  for(i=0; i<size; i++)
  {
   sum += yy[i];
   sum2 += yy[i]*yy[i];
  }
  sum = sum/size; //mean

  vv[7][SampleID] = sum; //mean
  vv[8][SampleID] = sum2/size - sum*sum; //var
 }
}

void component::SimCount( vector< vector< vector<int> > > & xsim )
{
 xsim.clear();
 int i,j, SampleID;
 for(i=0; i<NumPoint; i++)
 {
  xsim.push_back( vector< vector<int> >());
  for(j=0; j<NumPoint; j++)
  {
   xsim[i].push_back( vector<int>() );
   for(SampleID=0; SampleID<NumTissue; SampleID++)	xsim[i][j].push_back(0);
  }
 }

 for(SampleID=0; SampleID<NumTissue; SampleID++)
 {
  for(i=0; i<NumPoint-1; i++)
  {
   for(j=i+1; j<NumPoint; j++)
   {
    double nbmean = e[i][j][SampleID] * exp( theta[SampleID] * (data[i][j][SampleID]+1)/2 );
    xsim[i][j][SampleID] = gsl_ran_negative_binomial(rng, phi[SampleID]/(phi[SampleID]+nbmean), phi[SampleID]);
    xsim[j][i][SampleID] = xsim[i][j][SampleID];
   }
  }
 }
}

double component::CalcPr(vector< vector< vector<int> > > & zinput, vector<double> & a, int i, int j, int SampleID)
{
 int count = 0;
 if(i-1 >= 0)                         count += zinput[i-1][j][SampleID];
 if(i+1 <= NumPoint-1)        count += zinput[i+1][j][SampleID];
 if(j-1 >= 0)                         count += zinput[i][j-1][SampleID];
 if(j+1 <= NumPoint-1)        count += zinput[i][j+1][SampleID];

 vector<int> peak;      //index vector, the k th element is peak
 vector<int> back;      //index vector, the k th element is background

 peak = zinput[i][j];
 back = zinput[i][j];

 peak[SampleID] = 1;
 back[SampleID] = -1;

 int id, CountPeak = 0, CountBack = 0;
 for(id=0; id<NumTissue; id++)
 {
  if(peak[id] == 1)     CountPeak += (int) pow(2, id);
  if(back[id] == 1)     CountBack += (int) pow(2, id);
 }

 //peak log likelihood
 double valueP = x[i][j][SampleID] * theta[SampleID] - (x[i][j][SampleID] + phi[SampleID]) * log(e[i][j][SampleID] * exp(theta[SampleID]) + phi[SampleID]) - log(1.0 + exp(-2 * IsingGamma[SampleID] - 2 * IsingPsi[SampleID] * count) * a[CountBack]/a[CountPeak]);
 
 //background log likelihood
 double valueQ = -(x[i][j][SampleID] + phi[SampleID]) * log(e[i][j][SampleID] + phi[SampleID]) - log(1.0 + exp(2 * IsingGamma[SampleID] + 2 * IsingPsi[SampleID] * count) * alpha[CountPeak]/alpha[CountBack]);

 double p = 0.0;
 if(valueP >= valueQ) p = 1.0/(1.0 + exp(valueQ - valueP));
 else if(valueP < valueQ) p = exp(valueP - valueQ)/( 1.0 + exp(valueP - valueQ) );
 return p;
}

double component::PvalueCalSP(vector< vector< vector<int> > > & zinput, vector<double> & a, int i, int j)
{
 double B11 = 0.0;
 double B01 = 0.0;
 double B10 = 0.0;
 double B00 = 0.0;

 for(int SampleID=0; SampleID<NumTissue; SampleID++)
 {
  int count = 0;
  if(i-1 >= 0)                         count += zinput[i-1][j][SampleID];
  if(i+1 <= NumPoint-1)        count += zinput[i+1][j][SampleID];
  if(j-1 >= 0)                         count += zinput[i][j-1][SampleID];
  if(j+1 <= NumPoint-1)        count += zinput[i][j+1][SampleID];

  vector<int> peak;      //index vector, the k th element is peak
  vector<int> back;      //index vector, the k th element is background

  peak = zinput[i][j];
  back = zinput[i][j];

  peak[SampleID] = 1;
  back[SampleID] = -1;

  int id, CountPeak = 0, CountBack = 0;
  for(id=0; id<NumTissue; id++)
  {
   if(peak[id] == 1)     CountPeak += (int) pow(2, id);
   if(back[id] == 1)     CountBack += (int) pow(2, id);
  }

  //peak log likelihood
  double valueP = x[i][j][SampleID] * theta[SampleID] - (x[i][j][SampleID] + phi[SampleID]) * log(e[i][j][SampleID] * exp(theta[SampleID]) + phi[SampleID]) - log(1.0 + exp(-2 * IsingGamma[SampleID] - 2 * IsingPsi[SampleID] * count) * a[CountBack]/a[CountPeak]);

  //background log likelihood
  double valueQ = -(x[i][j][SampleID] + phi[SampleID]) * log(e[i][j][SampleID] + phi[SampleID]) - log(1.0 + exp(2 * IsingGamma[SampleID] + 2 * IsingPsi[SampleID] * count) * alpha[CountPeak]/alpha[CountBack]);
  
  if(SampleID == 0)
  {
   B11 += valueP;
   B10 += valueP;
   B01 += valueQ;
   B00 += valueQ;
  } else if(SampleID == 1)
  {
   B11 += valueP;
   B10 += valueQ;
   B01 += valueP;
   B00 += valueQ;
  }
 }

 pvalueSP[i][j] = 1.0/(1.0 + exp(B10-B11) + exp(B01-B11) + exp(B00-B11));

 double p11 = 0.0;
 p11 = log(1.0/(exp(B10-B11) + exp(B01-B11) + exp(B00-B11)));
 return p11;
}

// Pvalue for first-tissue specific
double component::PvalueCalFstS(vector< vector< vector<int> > > & zinput, vector<double> & a, int i, int j)
{
 double B11 = 0.0;
 double B01 = 0.0;
 double B10 = 0.0;
 double B00 = 0.0;

 for(int SampleID=0; SampleID<NumTissue; SampleID++)
 {
  int count = 0;
  if(i-1 >= 0)                         count += zinput[i-1][j][SampleID];
  if(i+1 <= NumPoint-1)        count += zinput[i+1][j][SampleID];
  if(j-1 >= 0)                         count += zinput[i][j-1][SampleID];
  if(j+1 <= NumPoint-1)        count += zinput[i][j+1][SampleID];

  vector<int> peak;      //index vector, the k th element is peak
  vector<int> back;      //index vector, the k th element is background

  peak = zinput[i][j];
  back = zinput[i][j];

  peak[SampleID] = 1;
  back[SampleID] = -1;

  int id, CountPeak = 0, CountBack = 0;
  for(id=0; id<NumTissue; id++)
  {
   if(peak[id] == 1)     CountPeak += (int) pow(2, id);
   if(back[id] == 1)     CountBack += (int) pow(2, id);
  }

  //peak log likelihood
  double valueP = x[i][j][SampleID] * theta[SampleID] - (x[i][j][SampleID] + phi[SampleID]) * log(e[i][j][SampleID] * exp(theta[SampleID]) + phi[SampleID]) - log(1.0 + exp(-2 * IsingGamma[SampleID] - 2 * IsingPsi[SampleID] * count) * a[CountBack]/a[CountPeak]);
 
  //background log likelihood
  double valueQ = -(x[i][j][SampleID] + phi[SampleID]) * log(e[i][j][SampleID] + phi[SampleID]) - log(1.0 + exp(2 * IsingGamma[SampleID] + 2 * IsingPsi[SampleID] * count) * alpha[CountPeak]/alpha[CountBack]);

  if(SampleID == 0)
  {
   B11 += valueP;
   B10 += valueP;
   B01 += valueQ;
   B00 += valueQ;
  } else if(SampleID == 1)
  {
   B11 += valueP;
   B10 += valueQ;
   B01 += valueP;
   B00 += valueQ;
  }
 }
  
 pvalueGS[i][j] = 1.0/(1.0 + exp(B11-B10) + exp(B01-B10) + exp(B00-B10));

 double p10 = 0.0;
 p10 = log(1.0/(exp(B11-B10) + exp(B01-B10) + exp(B00-B10)));
 return p10;
}

// Pvalue for second-tissue specific
double component::PvalueCalSecS(vector< vector< vector<int> > > & zinput, vector<double> & a, int i, int j)
{
 double B11 = 0.0;
 double B01 = 0.0;
 double B10 = 0.0;
 double B00 = 0.0;

 for(int SampleID=0; SampleID<NumTissue; SampleID++)
 {
  int count = 0;
  if(i-1 >= 0)                         count += zinput[i-1][j][SampleID];
  if(i+1 <= NumPoint-1)        count += zinput[i+1][j][SampleID];
  if(j-1 >= 0)                         count += zinput[i][j-1][SampleID];
  if(j+1 <= NumPoint-1)        count += zinput[i][j+1][SampleID];
 vector<int> peak;      //index vector, the k th element is peak
  vector<int> back;      //index vector, the k th element is background

  peak = zinput[i][j];
  back = zinput[i][j];

  peak[SampleID] = 1;
  back[SampleID] = -1;

  int id, CountPeak = 0, CountBack = 0;
  for(id=0; id<NumTissue; id++)
  {
   if(peak[id] == 1)     CountPeak += (int) pow(2, id);
   if(back[id] == 1)     CountBack += (int) pow(2, id);
  }

  //peak log likelihood
  double valueP = x[i][j][SampleID] * theta[SampleID] - (x[i][j][SampleID] + phi[SampleID]) * log(e[i][j][SampleID] * exp(theta[SampleID]) + phi[SampleID]) - log(1.0 + exp(-2 * IsingGamma[SampleID] - 2 * IsingPsi[SampleID] * count) * a[CountBack]/a[CountPeak]);

  //background log likelihood
  double valueQ = -(x[i][j][SampleID] + phi[SampleID]) * log(e[i][j][SampleID] + phi[SampleID]) - log(1.0 + exp(2 * IsingGamma[SampleID] + 2 * IsingPsi[SampleID] * count) * alpha[CountPeak]/alpha[CountBack]);

  if(SampleID == 0)
  {
   B11 += valueP;
   B10 += valueP;
   B01 += valueQ;
   B00 += valueQ;
  } else if(SampleID == 1)
  {
   B11 += valueP;
   B10 += valueQ;
   B01 += valueP;
   B00 += valueQ;
  }
 }

 pvalueIS[i][j] = 1.0/(1.0 + exp(B11-B01) + exp(B10-B01) + exp(B00-B01));

 double p01 = 0.0;
 p01 = log(1.0/(exp(B11-B01) + exp(B10-B01) + exp(B00-B01)));
 return p01;
}


void component::GibbsRefine(component & ExampleMode)
{
 int i, j, k, SampleID, Gibbs;
 double BestLL = 0;
 //int updated = 0; //indicates whether peak status(data[][][]) has been updated
 //printf("GibbsNum: %d", NumGibbs);
 //printf("\n");

 // keep record of four nuisance parameters in the Gibbs sample  
 vector< vector< vector<double> > > Record_Gibbs;
 for(i=0; i<NumTissue; i++)
 {
  Record_Gibbs.push_back( vector< vector<double> >() );
  for(j=0; j<NumGibbs; j++)
  {
   Record_Gibbs[i].push_back( vector<double>() );
   for(k=0; k<4; k++) Record_Gibbs[i][j].push_back(0); // 4 columns: Theta, Phi, Gamma, Psi;
  } 
 }

 //keep record of loglikelihood in the Gibbs sample
 vector<double> Record_LL;
 for(i=0; i<NumGibbs; i++) Record_LL.push_back(0);

 //get summary statistics for the input observed PLAC-Seq Contact matrix
 vector< vector<double> > statx;
 vector< vector< vector<double> > > StatRecord;

 for(j=0; j<9; j++)
 {
  statx.push_back( vector<double>() );
  for(SampleID=0; SampleID<NumTissue; SampleID++)	statx[j].push_back(0);
 }
 for(i=0; i< (int)(NumGibbs*0.8)+1; i++)
 {
  StatRecord.push_back( vector< vector<double> >() );
  for(j=0; j<9; j++)
  {
   StatRecord[i].push_back( vector<double>() );
   for(SampleID=0; SampleID<NumTissue; SampleID++)	StatRecord[i][j].push_back(0);
  }
 }

 SummaryStat(x, statx);
 for(SampleID=0; SampleID<NumTissue; SampleID++)
 {
  printf("%.01f\t%.01f\t%.01f\t%.01f\t%.01f\t%.01f\t%.01f\t%.01f\t%.01f\n",
         statx[0][SampleID], statx[1][SampleID], statx[2][SampleID], 
         statx[3][SampleID], statx[4][SampleID], statx[5][SampleID], 
         statx[6][SampleID], statx[7][SampleID], statx[8][SampleID]);
  for(j=0; j<9; j++)	StatRecord[0][j][SampleID] = statx[j][SampleID];	//first row is the real data
 }
 
 vector< vector< vector<int> > > xsim;
 vector< vector<double> > statxsim;

 double loglikefull_initial = 0;
 for(SampleID = 0; SampleID<NumTissue; SampleID++)
 {
  loglikefull_initial += GetNB(data, theta[SampleID], phi[SampleID], SampleID);
  loglikefull_initial += GetLogLikelihood(data, alpha, IsingGamma[SampleID], IsingPsi[SampleID], SampleID);
 }
 printf("initial loglikefull = %.4lf\n", loglikefull_initial); 
 BestLL = loglikefull_initial;

 //mcmc
 for(Gibbs=0; Gibbs<NumGibbs; Gibbs++) //do Gibbs sample, j -> GibbsID
 {
  loglikefull = 0;
  
   //update Z 
   for(int i=0; i<NumPoint-1; i++)
   {
    for(int j=i+1; j<NumPoint; j++){
      //double ratio = x[i][j][SampleID]/e[i][j][SampleID];
      //if(ratio>1.1 && ratio<2.0){
      GibbsUpdate(data, alpha, i, j);
      //updated = 1;
      //}
    }
   }
   //printf("Best Theta = %.4lf\n", theta[SampleID]);
   
   //update four nuisance parameters
  // if (updated == 1)
  // {
  if(Gibbs%100==0)
  {
   for(SampleID = 0; SampleID<NumTissue; SampleID++){
   theta[SampleID] = UpdateTheta(SampleID);
   phi[SampleID] = UpdatePhi(SampleID);
   IsingGamma[SampleID] = UpdateGamma(SampleID);
   IsingPsi[SampleID] = UpdatePsi(SampleID);
   }
  }
  // updated = 0;

  for(SampleID = 0; SampleID<NumTissue; SampleID++){
   loglikefull += GetNB(data, theta[SampleID], phi[SampleID], SampleID);
   loglikefull += GetLogLikelihood(data, alpha, IsingGamma[SampleID], IsingPsi[SampleID], SampleID);
  }
//   printf("check point %d	%.4lf	%.4lf\n", GibbsID, loglikealpha, loglikefull);
  
  for(SampleID=0; SampleID<NumTissue; SampleID++)
  {
   Record_Gibbs[SampleID][Gibbs][0] = theta[SampleID];
   Record_Gibbs[SampleID][Gibbs][1] = phi[SampleID]; 
   Record_Gibbs[SampleID][Gibbs][2] = IsingGamma[SampleID];
   Record_Gibbs[SampleID][Gibbs][3] = IsingPsi[SampleID];
  }
  
  Record_LL[Gibbs] = loglikefull;

//printf("%d\t%d\n", Gibbs+1, (int)(NumGibbs*0.2));
  //if(Gibbs == 0) BestLL = loglikefull;
  //else if(Gibbs >= 1)
  //{
   if( loglikefull > BestLL ) 
   {
    BestLL = loglikefull;
    for(i=0; i<NumPoint-1; i++)
    {
     for(j=i+1; j<NumPoint; j++)
     {
      for(SampleID=0; SampleID<NumTissue; SampleID++)
      {
       dataMode[i][j][SampleID] = data[i][j][SampleID];
       dataMode[j][i][SampleID] = dataMode[i][j][SampleID];
      }
     } 
    }
   }
  //}
  
  //do posterior predictive check
  if(Gibbs+1 > NumGibbs*0.2)
  {
   //simulate PLAC-seq contact frequency
   SimCount(xsim);

   //get summary statistics for the simulated count data
   SummaryStat(xsim, statxsim);
   
   //keep record
   for(SampleID=0; SampleID<NumTissue; SampleID++)
   {
    for(j=0; j<9; j++) 
    {
     StatRecord[Gibbs-int(NumGibbs*0.2)+1][j][SampleID] = statxsim[j][SampleID];
    }
   }
  }

  if(Gibbs%1000==0) printf("Gibbs = %d, prop = %.4lf,  BestLL = %.4lf\n", Gibbs, (Gibbs+0.0)/NumGibbs, BestLL);
 }

 for(SampleID=0; SampleID<NumTissue; SampleID++)
 {
  for(i=0; i<NumPoint-1; i++)
  {
   for(j=i+1; j<NumPoint; j++){
     p[i][j][SampleID] = CalcPr(dataMode, alpha, i, j, SampleID);
     p[j][i][SampleID] = p[i][j][SampleID];
   }
  }
 }

 strcpy( StrOutCur, StrOutDir );
 FILE * s = fopen(strcat(StrOutCur, "/Record_Gibbs.txt"), "w");
 for(Gibbs=0; Gibbs<NumGibbs; Gibbs++)
 {
  for(SampleID=0; SampleID<NumTissue; SampleID++)
  {
   for(k=0; k<4; k++) fprintf(s, "%.4lf\t", Record_Gibbs[SampleID][Gibbs][k]); 
  }
  fprintf(s, "%.4lf\n", Record_LL[Gibbs]);
 }
 fclose(s);

 //output peak status
 strcpy( StrOutCur, StrOutDir );
 FILE * sz = fopen(strcat(StrOutCur, "/Record_long_format.txt"), "w");
 for(SampleID=0; SampleID<NumTissue; SampleID++)
 {
  for(i=0; i<NumPoint-1; i++)
  {
   for(j=i+1; j<NumPoint; j++)	fprintf(sz, "%d\t%d\t%d\t%d\t%.8lf\t%d\t%.8lf\t\n", SampleID, i*binsize+bininitial, j*binsize+bininitial, x[i][j][SampleID], e[i][j][SampleID], dataMode[i][j][SampleID], p[i][j][SampleID]);
 }
}
 fclose(sz);

 //calculate peak likelihood
 for(i=0; i<NumPoint-1; i++)
 {
  for(j=i+1; j<NumPoint; j++)
  {
   loglikeSP[i][j] = PvalueCalSP(dataMode, alpha, i, j);
   loglikeGS[i][j] = PvalueCalFstS(dataMode, alpha, i, j);
   loglikeIS[i][j] = PvalueCalSecS(dataMode, alpha, i, j);
  }
 }
 
 //output peak likelihood
 strcpy( StrOutCur, StrOutDir );
 FILE * sp = fopen(strcat(StrOutCur, "/PP.txt"), "w");
 for(i=0; i<NumPoint-1; i++)
 {
  for(j=i+1; j<NumPoint; j++)  fprintf(sp, "%d\t%d\t%d\t%.4lf\t%d\t%1.6e\t%d\t%.4lf\t%d\t%1.6e\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\n", i*binsize+bininitial, j*binsize+bininitial, x[i][j][0], e[i][j][0], dataMode[i][j][0], pvalue[i][j][0], x[i][j][1], e[i][j][1], dataMode[i][j][1], pvalue[i][j][1], pvalueSP[i][j], pvalueGS[i][j], pvalueIS[i][j], loglikeSP[i][j], loglikeGS[i][j], loglikeIS[i][j]);
 }
 fclose(sp);

 strcpy( StrOutCur, StrOutDir );
 FILE * sest = fopen(strcat(StrOutCur, "/Record_Para.txt"), "w"); 
 fprintf(sest, "Best LogLike = %.4lf\n", BestLL);
 for(SampleID=0; SampleID<NumTissue; SampleID++)
 {
  fprintf(sest, "Tissue = %d\n", SampleID);
  fprintf(sest, "Best Theta = %.4lf\n", theta[SampleID]);
  fprintf(sest, "Best Phi = %.4lf\n", phi[SampleID]);
  fprintf(sest, "Best Gamma = %.4lf\n", IsingGamma[SampleID]);
  fprintf(sest, "Best Psi = %.4lf\n", IsingPsi[SampleID]);
 }
 fprintf(sest, "NumGibbs = %d\n", NumGibbs);
 fprintf(sest, "SEED = %d\n", seed);
 fclose(sest);
 
 //output posterior predictive check p-values
 vector< vector<double> > StatCount, StatCountTie;
 for(i=0; i<9; i++)
 {
  StatCount.push_back( vector<double>() );
  StatCountTie.push_back( vector<double>() );
  for(SampleID=0; SampleID<NumTissue; SampleID++)
  {
   StatCount[i].push_back(0);
   StatCountTie[i].push_back(0);
  }
 }

 for(SampleID=0; SampleID<NumTissue; SampleID++)
 {
  for(i=1; i< (int)StatRecord.size(); i++)
  {
   for(j=0; j<9; j++)
   {
    if( StatRecord[0][j][SampleID] >= StatRecord[i][j][SampleID] ) StatCountTie[j][SampleID]++;
    if( StatRecord[0][j][SampleID] >  StatRecord[i][j][SampleID] ) StatCount[j][SampleID]++;
   }
  }
 }

 //output posterior predictive check p-values
 strcpy( StrOutCur, StrOutDir );
 FILE * sppc = fopen(strcat(StrOutCur, "/Record_PPC_Pvalue.txt"),"w");
 for(SampleID=0; SampleID<NumTissue; SampleID++)
 {
  fprintf(sppc, "Tissue = %d\n", SampleID);
  fprintf(sppc, "p-value_of_min = %.4lf\n",             ( StatCountTie[0][SampleID] + StatCount[0][SampleID] )/( (StatRecord.size()-1)*2 ) );
  fprintf(sppc, "p-value_of_10_percentile = %.4lf\n",   ( StatCountTie[1][SampleID] + StatCount[1][SampleID] )/( (StatRecord.size()-1)*2 ) );
  fprintf(sppc, "p-value_of_25_percentile = %.4lf\n",   ( StatCountTie[2][SampleID] + StatCount[2][SampleID] )/( (StatRecord.size()-1)*2 ) );
  fprintf(sppc, "p-value_of_median = %.4lf\n",          ( StatCountTie[3][SampleID] + StatCount[3][SampleID] )/( (StatRecord.size()-1)*2 ) );
  fprintf(sppc, "p-value_of_75_percentile = %.4lf\n",   ( StatCountTie[4][SampleID] + StatCount[4][SampleID] )/( (StatRecord.size()-1)*2 ) );
  fprintf(sppc, "p-value_of_90_percentile = %.4lf\n",   ( StatCountTie[5][SampleID] + StatCount[5][SampleID] )/( (StatRecord.size()-1)*2 ) );
  fprintf(sppc, "p-value_of_max = %.4lf\n",             ( StatCountTie[6][SampleID] + StatCount[6][SampleID] )/( (StatRecord.size()-1)*2 ) );
  fprintf(sppc, "p-value_of_mean = %.4lf\n",            ( StatCountTie[7][SampleID] + StatCount[7][SampleID] )/( (StatRecord.size()-1)*2 ) );
  fprintf(sppc, "p-value_of_variance = %.4lf\n",        ( StatCountTie[8][SampleID] + StatCount[8][SampleID] )/( (StatRecord.size()-1)*2 ) );
 }
 fclose(sppc);    

 strcpy( StrOutCur, StrOutDir );
 sppc = fopen(strcat(StrOutCur,"/Record_PPC_stat.txt"),"w");
 for(SampleID=0; SampleID<NumTissue; SampleID++)
 {
  for(i=0; i< (int)StatRecord.size(); i++)
  {
   fprintf(sppc, "Tissue = %d\n", SampleID);
   for(j=0; j<7; j++) fprintf(sppc, "%.0lf\t", StatRecord[i][j][SampleID]);
   fprintf(sppc, "%.4lf    ", StatRecord[i][7][SampleID]);
   fprintf(sppc, "%.4lf\n",  StatRecord[i][8][SampleID]);
  }
 }
 fclose(sppc);

 //QC datamode
 vector<int> countmode;
 for(i=0; i<(int)alpha.size(); i++) countmode.push_back(0); // count the frequency of each configuration
 
 for(i=0; i<NumPoint-1; i++)
 {
  for(j=i+1; j<NumPoint; j++)
  {
   int numpeak=0;
   for(k=0; k<NumTissue; k++)
   {
    if(dataMode[i][j][k]==1) numpeak += (int) pow(2, k);
   }
   countmode[numpeak]++;
  }
 }
 
 strcpy( StrOutCur, StrOutDir );
 FILE *sprob = fopen(strcat(StrOutCur, "/Record_prob.txt"), "w");
 for(i=0; i<(int)alpha.size(); i++) fprintf(sprob, "%d\t%d\t%.8lf\n", i, countmode[i], countmode[i]/( NumPoint*(NumPoint-1)/2.0 ));
 fclose(sprob);
}
