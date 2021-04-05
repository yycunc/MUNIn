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

  theta.push_back(0);
  phi.push_back(0);
  IsingGamma.push_back(0);
  IsingPsi.push_back(0); 
 }
 
 for(i=0; i<NumTissue; i++){
  theta[i] = theta_orig[i];
  phi[i] = phi_orig[i];
  IsingGamma[i] = IsingGamma_orig[i];
  IsingPsi[i] = IsingPsi_orig[i];
 }
 
  
 loglikefull=0;
}

void component::InitialZ()
{
 int i, j, SampleID;

 for(i=0; i<NumPoint; i++)
 {
  data.push_back( vector< vector<int> >() );
  for(j=0; j<NumPoint; j++)
  {
   data[i].push_back( vector<int>() );
   for(SampleID=0; SampleID<NumTissue; SampleID++) data[i][j].push_back(0);
  }
 }

 for(i=0; i<NumPoint; i++)
 {
  for(j=0; j<NumPoint; j++)
  {
   for(SampleID=0; SampleID<NumTissue; SampleID++) 
   {
    data[i][j][SampleID] = dataInitial[i][j][SampleID];
    data[j][i][SampleID] = data[i][j][SampleID];
   }
  }
 }
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

double component::GetNB_single( vector< vector < vector<int> > > & zinput, double ValueTheta, double ValuePhi, int i, int j, int SampleID )
{
 double sum = 0;

 double nbmean = e[i][j][SampleID] * exp( ValueTheta * (zinput[i][j][SampleID]+1)/2 );
 sum = gsl_sf_lngamma(x[i][j][SampleID]+ValuePhi) - gsl_sf_lngamma(ValuePhi) + ValuePhi*log(ValuePhi) + x[i][j][SampleID]*log(nbmean) - (x[i][j][SampleID]+ValuePhi)*log(nbmean+ValuePhi);

 return sum;
}

double component::GetLogLikelihood_single( vector< vector < vector<int> > > & zinput, vector<double> & a, double ValueGamma, double ValuePsi, int i, int j, int SampleID )
{
 double sum = 0;
 int ii, jj;
 int i_start = 0;
 int i_end = NumPoint-1;
 int j_start = 0;
 int j_end = NumPoint-1;
 
 if(i-1 >= 0)	i_start = i-1;
 if(i+1 <= NumPoint-1)	i_end = i+1;
 if(j-1 >= 0)	j_start = j-1;
 if(j+1 <= NumPoint-1)	j_end = j+1;

 for(ii=i_start; ii<i_end; ii++)
 {
  for(jj=j_start; jj<j_end; jj++)
  {
   int count = 0;
   if(ii-1 >= 0)                         count += zinput[ii-1][jj][SampleID];
   if(ii+1 <= NumPoint-1)        count += zinput[ii+1][jj][SampleID];
   if(jj-1 >= 0)                         count += zinput[ii][jj-1][SampleID];
   if(jj+1 <= NumPoint-1)        count += zinput[ii][jj+1][SampleID];

   vector<int> peak; // index vector, the k th element is peak
   vector<int> back; // index vector, the k th element is background

   peak = zinput[ii][jj];
   back = zinput[ii][jj];

   peak[SampleID] = 1;
   back[SampleID] = -1;

   int id, CountPeak=0, CountBack=0;
   for(id=0; id<NumTissue; id++)
   {
    if(peak[id] == 1)     CountPeak += (int) pow(2, id);
    if(back[id] == 1)     CountBack += (int) pow(2, id);
   }
   sum += -log( 1.0 + exp(-2 * ValueGamma * zinput[ii][jj][SampleID] - 2 * ValuePsi * zinput[ii][jj][SampleID] * count) * pow(a[CountBack]/a[CountPeak], zinput[ii][jj][SampleID]) );
  }
 }

 return sum;
}

void component::GibbsUpdate(vector< vector< vector<int> > > & zinput, vector<double> & a, int i, int j)
{
 int k, l;
/*
 double B11 = 0.0;
 double B01 = 0.0;
 double B10 = 0.0;
 double B00 = 0.0;
 int k, l;

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
*/
 vector<int> TmpValue;
 vector<double> LL;
 vector<double> PP;
 vector< vector<int> > TmpValueRecord;
 for(k=0; k<PPSize; k++)
 {
  TmpValueRecord.push_back( vector<int>() );
  for(l=0; l<NumTissue; l++)
  {
   TmpValueRecord[k].push_back(0);
  }
 }

 for(l=0; l<NumTissue; l++)   TmpValue.push_back(0);
 for(k=0; k<PPSize; k++)
 {
  LL.push_back(0);
  PP.push_back(0);
 }
  
 for(k=0; k<PPSize; k++)
 {
  //convert decimal to binary, from i (decimal) -> TmpValue (binary, length 4)
  int Value = k;
  for(l=NumTissue-1; l>=0; l--)
  {
   TmpValue[l] = (int) floor( Value / pow(2, l) );
   Value = Value - TmpValue[l]*(int)pow(2, l);
  }
  //fprintf(salpha, "%d\t", i);
  for(l=0; l<NumTissue; l++)  TmpValueRecord[k][l] = TmpValue[l];

  double CountLL = 0;
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
   
   CountLL += TmpValue[SampleID] * valueP;
   CountLL += (1.0 - TmpValue[SampleID]) * valueQ; //count the number of tissues with peak, defined as S
   //printf("%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n", TmpValue[SampleID], valueP, valueQ, TmpValue[SampleID] * valueP, CountLL); 
  }

  LL[k] = CountLL;
  //printf("%.4lf\n", LL[k]);
 }

 //double CountP = 0;
 for(k=0; k<PPSize; k++)
 {
  double CountP = 0;
  for(l=0; l<PPSize; l++)
  {
   if(l != k) CountP += exp(LL[l] - LL[k]);
  }

  PP[k] = 1.0/(1.0 + CountP);
 }

 double value = gsl_rng_uniform(rng);
 double start = 0;
 for(k=0; k<PPSize; k++)
 {
  if(k == 0){
    if(value >= start && value <= start + PP[k]){
      for(l=0; l<NumTissue; l++)
      {
       zinput[i][j][l] = TmpValueRecord[k][l]*2-1;
       zinput[j][i][l] = zinput[i][j][l];
      }
      //printf("%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%d\t%d\t%d\t%d\n", start, start+PP[k], value, PP[0], PP[1], PP[2], PP[3], TmpValueRecord[k][0], TmpValueRecord[k][1], zinput[i][j][0], zinput[i][j][1]);
      break; 
    } else{
      start = start + PP[k];
    }
  } else{
    if(value > start && value <= start + PP[k]){
      for(l=0; l<NumTissue; l++)
      {
       zinput[i][j][l] = TmpValueRecord[k][l]*2-1;
       zinput[j][i][l] = zinput[i][j][l];
      }
      //printf("%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%d\t%d\t%d\t%d\n", start, start+PP[k], value, PP[0], PP[1], PP[2], PP[3], TmpValueRecord[k][0], TmpValueRecord[k][1], zinput[i][j][0], zinput[i][j][1]);
      break; 
    } else{
      start = start + PP[k];
    }
  }
 }    

/*
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
  //double value = gsl_rng_uniform(rng);
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
*/
 //printf("%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%d\t%d\n", PP[0], PP[1], PP[2], PP[3], pvalue11, pvalue10, pvalue01, pvalue00, value, zinput[i][j][0], zinput[i][j][1]);

}

double component::UpdateTheta(int SampleID)
{
 double Theta = theta[SampleID];
 double ValueOld = GetNB(data, Theta, phi[SampleID], SampleID);
 double range = theta_orig[SampleID]*0.2;
 double stepsize = range/20;
 int i = 1;

 double ThetaNew = theta_orig[SampleID] - range;
 if(ThetaNew < 0)
 {
  ThetaNew = -ThetaNew;
 }
 double ValueNew = GetNB(data, ThetaNew, phi[SampleID], SampleID);
 while(i<=40)
 {
  double ThetaUpdated = theta_orig[SampleID] - range + i * stepsize;
  if(ThetaUpdated < 0)
  {
   ThetaUpdated = -ThetaUpdated;
  }
  double ValueUpdated = GetNB(data, ThetaUpdated, phi[SampleID], SampleID);
  if(ValueUpdated >= ValueNew)
  {
   ThetaNew = ThetaUpdated;
   ValueNew = ValueUpdated;
  }
  i++;
 }

 if(ValueNew >= ValueOld)
 {
  Theta = ThetaNew;
 }
 else if(ValueNew < ValueOld)
 {
  double value = gsl_rng_uniform(rng);
  if(value < exp(ValueNew - ValueOld))
  {
   Theta = ThetaNew;
  }
 }

 return Theta;
}

double component::UpdatePhi(int SampleID)
{
 double Phi = phi[SampleID];
 double ValueOld = GetNB(data, theta[SampleID], Phi, SampleID);
 double range = phi_orig[SampleID]*0.2;
 double stepsize = range/20;
 int i = 1;

 double PhiNew = phi_orig[SampleID] - range;
 if(PhiNew < 0)
 {
  PhiNew = -PhiNew;
 }
 double ValueNew = GetNB(data, theta[SampleID], PhiNew, SampleID);
 while(i<=40)
 {
  double PhiUpdated = phi_orig[SampleID] - range + i * stepsize;
  if(PhiNew < 0)
  {
   PhiUpdated = -PhiUpdated;
  }
  double ValueUpdated = GetNB(data, theta[SampleID], PhiUpdated, SampleID);
  if(ValueUpdated >= ValueNew)
  {
   PhiNew = PhiUpdated;
   ValueNew = ValueUpdated;
  }
  i++;
 }

 if(ValueNew >= ValueOld)
 {
  Phi = PhiNew;
 }
 else if(ValueNew < ValueOld)
 {
  double value = gsl_rng_uniform(rng);
  if(value < exp(ValueNew - ValueOld))
  {
   Phi = PhiNew;
  }
 }

 return Phi;
}

double component::UpdateGamma(int SampleID)
{
 double Gamma = IsingGamma[SampleID];
 double ValueOld = GetLogLikelihood(data, alpha, Gamma, IsingPsi[SampleID], SampleID);
 double range = -IsingGamma_orig[SampleID]*0.2;
 double stepsize = range/20;
 int i = 1;

 double GammaNew = IsingGamma_orig[SampleID] - range;
 if(GammaNew > 0)
 {
  GammaNew = -GammaNew;
 }
 double ValueNew = GetLogLikelihood(data, alpha, GammaNew, IsingPsi[SampleID], SampleID);
 while(i<=40)
 {
  double GammaUpdated = IsingGamma_orig[SampleID] - range + i * stepsize;
  if(GammaNew > 0)
  {
   GammaUpdated = -GammaUpdated;
  }
  double ValueUpdated = GetLogLikelihood(data, alpha, GammaUpdated, IsingPsi[SampleID], SampleID);
  if(ValueUpdated >= ValueNew)
  {
   GammaNew = GammaUpdated;
   ValueNew = ValueUpdated;
  }
  i++;
 }

 if(ValueNew >= ValueOld)
 {
  Gamma = GammaNew;
 }
 else if(ValueNew < ValueOld)
 {
  double value = gsl_rng_uniform(rng);
  if(value < exp(ValueNew - ValueOld))
  {
   Gamma = GammaNew;
  }
 }

 return Gamma;
}

double component::UpdatePsi(int SampleID)
{
 double Psi = IsingPsi[SampleID];
 double ValueOld = GetLogLikelihood(data, alpha, IsingGamma[SampleID], Psi, SampleID);
 double range = IsingPsi_orig[SampleID]*0.2;
 double stepsize = range/20;
 int i = 1;

 double PsiNew = IsingPsi_orig[SampleID] - range;
 if(PsiNew < 0)
 {
  PsiNew = -PsiNew;
 }
 double ValueNew = GetLogLikelihood(data, alpha, IsingGamma[SampleID], PsiNew, SampleID);
 while(i<=40)
 {
  double PsiUpdated = IsingPsi_orig[SampleID] - range + i * stepsize;
  if(PsiNew < 0)
  {
   PsiUpdated = -PsiUpdated;
  }
  double ValueUpdated = GetLogLikelihood(data, alpha, IsingGamma[SampleID], PsiUpdated, SampleID);
  if(ValueUpdated >= ValueNew)
  {
   PsiNew = PsiUpdated;
   ValueNew = ValueUpdated;
  }
  i++;
 }
 
 if(ValueNew >= ValueOld)
 {
  Psi = PsiNew;
 }
 else if(ValueNew < ValueOld)
 { 
  double value = gsl_rng_uniform(rng);
  if(value < exp(ValueNew - ValueOld))
  {
   Psi = PsiNew;
  } 
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

void component::PPModeCar(vector< vector< vector<int> > > & zinput, vector<double> & a, int i, int j)
{
 int k, l;
 //initialize PP vector
 vector<int> TmpValue;
 vector<double> LL;
 vector<double> PP;
 vector< vector<int> > TmpValueRecord;
 for(k=0; k<PPSize; k++)
 {
  TmpValueRecord.push_back( vector<int>() );
  for(l=0; l<NumTissue; l++)   TmpValueRecord[k].push_back(0);
  LL.push_back(0);
  PP.push_back(0);
 }
 for(l=0; l<NumTissue; l++)   TmpValue.push_back(0);

 for(k=0; k<PPSize; k++)
 {
  //convert decimal to binary, from i (decimal) -> TmpValue (binary, length 4)
  int Value = k;
  for(l=NumTissue-1; l>=0; l--)
  {
   TmpValue[l] = (int) floor( Value / pow(2, l) );
   Value = Value - TmpValue[l]*(int)pow(2, l);
  }
  for(l=0; l<NumTissue; l++)  TmpValueRecord[k][l] = TmpValue[l];
  double CountLL = 0;
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
   
   CountLL += TmpValue[SampleID] * valueP;
   CountLL += (1.0 - TmpValue[SampleID]) * valueQ; //count the number of tissues with peak, defined as S
   //printf("%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n", TmpValue[SampleID], valueP, valueQ, TmpValue[SampleID] * valueP, CountLL);
  }

  LL[k] = CountLL;
 }

 for(k=0; k<PPSize; k++)
 {
  double CountP = 0;
  for(l=0; l<PPSize; l++)
  {
   if(l != k) CountP += exp(LL[l] - LL[k]);
  }
  PP[k] = 1.0/(1.0 + CountP);
  PPMode[k][i][j] = PP[k];
  loglikeMode[k][i][j] = log(1.0/CountP);
 }
}

void component::GibbsRefine(component & ExampleMode)
{
 int i, j, k, m, SampleID, Gibbs;
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
   for(int j=i+1; j<NumPoint; j++)
   {      
    GibbsUpdate(data, alpha, i, j);
   }
  }
  //printf("Best Theta = %.4lf\n", theta[SampleID]);
  
  if(Gibbs%100==0)
  {
   for(SampleID = 0; SampleID<NumTissue; SampleID++)
   {
    theta[SampleID] = UpdateTheta(SampleID);
    phi[SampleID] = UpdatePhi(SampleID);
    IsingGamma[SampleID] = UpdateGamma(SampleID);
    IsingPsi[SampleID] = UpdatePsi(SampleID);
   }
  }
  
  for(SampleID = 0; SampleID<NumTissue; SampleID++)
  {
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
   ExampleMode = *this;
  }
  else if( loglikefull < BestLL )
  {
   double value = gsl_rng_uniform(rng);
   if(value < exp(loglikefull - BestLL))
   {
    BestLL = loglikefull;
    ExampleMode = *this;
   }
  }

  for(m=0; m<num_lowObs; m++)
  {
   i = lowObs[m][0];
   j = lowObs[m][1];
   SampleID = lowObs[m][2];
   //if( x[i][j][lowObs[i][2]] > e[lowObs[i][0]][lowObs[i][1]][lowObs[i][2]] )
   //{
   // printf ("x = %d, e = %.4lf\n", x[lowObs[i][0]][lowObs[i][1]][lowObs[i][2]], e[lowObs[i][0]][lowObs[i][1]][lowObs[i][2]]);
   //}
   if( data[i][j][SampleID] == 1 )
   {
    double loglikefull_orig = 0;
    for(int k = 0; k<NumTissue; k++)
    {
     loglikefull_orig += GetNB_single(data, theta[k], phi[k], i, j, k);
     loglikefull_orig += GetLogLikelihood_single(data, alpha, IsingGamma[k], IsingPsi[k], i, j, k);
    }

    data[i][j][SampleID] = -1;
    data[j][i][SampleID] = data[i][j][SampleID];
    double loglikefull_after = 0;
    for(int k = 0; k<NumTissue; k++)
    {
     loglikefull_after += GetNB_single(data, theta[k], phi[k], i, j, k);
     loglikefull_after += GetLogLikelihood_single(data, alpha, IsingGamma[k], IsingPsi[k], i, j, k);
    }
    if( loglikefull_after > loglikefull_orig )
    {
     loglikefull += loglikefull_after - loglikefull_orig; 
     if( loglikefull > BestLL)
     {
      BestLL = loglikefull;
      ExampleMode = *this;
     }
     else if( loglikefull < BestLL )
     {
      double value = gsl_rng_uniform(rng);
      if(value < exp(loglikefull - BestLL))
      {
       BestLL = loglikefull;
       ExampleMode = *this;
      }
     }
     else
     {
      data[i][j][SampleID] = 1;
      data[j][i][SampleID] = data[i][j][SampleID];
     }
    }
    else
    {
     data[i][j][SampleID] = 1;
     data[j][i][SampleID] = data[i][j][SampleID];
    }
   }
  }

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
     p[i][j][SampleID] = CalcPr(ExampleMode.data, alpha, i, j, SampleID);
     p[j][i][SampleID] = p[i][j][SampleID];
   }
  }
 }

 //calculate posterior probabilty
 for(i=0; i<NumPoint-1; i++)
 {
  for(j=i+1; j<NumPoint; j++)   PPModeCar(ExampleMode.data, alpha, i, j);
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
   for(j=i+1; j<NumPoint; j++)	fprintf(sz, "%d\t%d\t%d\t%d\t%.8lf\t%d\t%.8lf\t\n", SampleID, i*binsize+bininitial, j*binsize+bininitial, x[i][j][SampleID], e[i][j][SampleID], ExampleMode.data[i][j][SampleID], p[i][j][SampleID]);
  }
 }
 fclose(sz);
 
 //output peak likelihood
 strcpy( StrOutCur, StrOutDir );
 FILE * sp = fopen(strcat(StrOutCur, "/PP.txt"), "w");
 for(i=0; i<NumPoint-1; i++)
 {
  for(j=i+1; j<NumPoint; j++)
  {
   fprintf(sp, "%d\t%d", i*binsize+bininitial, j*binsize+bininitial);
   for(SampleID=0; SampleID<NumTissue; SampleID++)  fprintf(sp, "\t%d\t%.4lf\t%d\t%1.6e", x[i][j][SampleID], e[i][j][SampleID], ExampleMode.data[i][j][SampleID], pvalue[i][j][SampleID]);
   for(k=0; k<PPSize; k++)  fprintf(sp, "\t%.8lf", PPMode[k][i][j]);
   for(k=0; k<PPSize; k++)  fprintf(sp, "\t%.8lf", loglikeMode[k][i][j]);
   fprintf(sp, "\n");
  }
 }
 fclose(sp);

 strcpy( StrOutCur, StrOutDir );
 FILE * sest = fopen(strcat(StrOutCur, "/Record_Para.txt"), "w"); 
 fprintf(sest, "Best LogLike = %.4lf\n", BestLL);
 for(SampleID=0; SampleID<NumTissue; SampleID++)
 {
  fprintf(sest, "Tissue = %d\n", SampleID);
  fprintf(sest, "Best Theta = %.16lf\n", ExampleMode.theta[SampleID]);
  fprintf(sest, "Best Phi = %.16lf\n", ExampleMode.phi[SampleID]);
  fprintf(sest, "Best Gamma = %.16lf\n", ExampleMode.IsingGamma[SampleID]);
  fprintf(sest, "Best Psi = %.16lf\n", ExampleMode.IsingPsi[SampleID]);
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
    if(ExampleMode.data[i][j][k]==1) numpeak += (int) pow(2, k);
   }
   countmode[numpeak]++;
  }
 }
 
 strcpy( StrOutCur, StrOutDir );
 FILE *sprob = fopen(strcat(StrOutCur, "/Record_prob.txt"), "w");
 for(i=0; i<(int)alpha.size(); i++) fprintf(sprob, "%d\t%d\t%.8lf\n", i, countmode[i], countmode[i]/( NumPoint*(NumPoint-1)/2.0 ));
 fclose(sprob);
}
