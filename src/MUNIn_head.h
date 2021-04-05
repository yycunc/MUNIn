
class TuningPara
{
public: 
double step, large, small; //tuning parameters in HMC for position

TuningPara & operator=(const TuningPara & rhs); 

void UpdateTuningPara(const vector<int> & acceptance, int j, int size, double lowrate, double highrate);

//constructor and destructor
TuningPara(void):step(1),large(10),small(0.000001){}
~TuningPara(){}
};


void TuningPara::UpdateTuningPara(const vector<int> & acceptance, int current_j, int size, double lowrate, double highrate)
{
 int k;
 double sum=0;
 for(k=current_j-size+1; k<=current_j; k++) sum += acceptance[k];
 
 if(sum<lowrate*size)
   {
    large = step;
	step = (step+small)/2;
   }
 if(sum>highrate*size)
   {
    small = step;
	step = (step+large)/2;
   }
 if(large-small< large*0.1)
   {
    large *= 10;
	small *= 0.1;
   }
 if(large > 10)
   {
    large = 10;
   }
 if(small < 0.01)
   {
    small = 0.01;
   }
}

TuningPara & TuningPara::operator=(const TuningPara & rhs)
{
 step = rhs.step;
 large = rhs.large;
 small = rhs.small;

 return *this;
}

class component{
public:	
    double loglikefull; // full likelihood 
    vector<double> theta;
    vector<double> phi;
    vector<double> IsingGamma;
    vector<double> IsingPsi;
    //vector< vector< vector<int> > > data; //peak indicator: z=1 peak, z=-1 background
    
    //keep track of Metropolis-Hasting acceptance rate
    vector<int> CountTheta;
    vector<int> CountPhi;
    vector<int> CountGamma;
    vector<int> CountPsi;

    //tune RW MH parametres, acceptance ratio 20% ~ 40%
    vector< TuningPara > TuneTheta; 
    vector< TuningPara > TunePhi;
    vector< TuningPara > TuneGamma;	
    vector< TuningPara > TunePsi;	

    vector< vector< vector<int> > > data;

    component & operator=(const component & rhs);	 
    
    //Calculating the median value of the x/e ratio
    double CalcMedian(vector<double> & scores);

    //initialization of peak indicator z
    void InitialZ();
	
    //initialization of other nuisance parameters
    void InitialPara();
		
    //update peak indicator z, binary sample Z, for loci pair (j,k) of each tissue i
    void GibbsUpdate(vector< vector< vector<int> > > & zinput, vector<double> & a, int i, int j);
	
    //calculate log NB likelihood
    double GetNB(vector< vector < vector<int> > > & zinput, double ValueTheta, double ValuePhi, int SampleID);

    //calculate log Ising prior
    double GetLogLikelihood(vector< vector < vector<int> > > & zinput, vector<double> & a, double ValueGamma, double ValuePsi, int SampleID);

    //calculate log NB likelihood
    double GetNB_single(vector< vector < vector<int> > > & zinput, double ValueTheta, double ValuePhi, int i, int j, int SampleID);
    
    //calculate log Ising prior
    double GetLogLikelihood_single(vector< vector < vector<int> > > & zinput, vector<double> & a, double ValueGamma, double ValuePsi, int i, int j, int SampleID);

    //RW MH update Theta
    double UpdateTheta(int SampleID);

    //RW MH update Phi
    double UpdatePhi(int SampleID);

    //RW MH update Gamma
    double UpdateGamma(int SampleID);

    //RW MH update Psi
    double UpdatePsi(int SampleID);
	
    //Gibbs sample
    void GibbsRefine(component & ExampleMode);	

    //Posterior predictive check
    void SimCount(vector<vector< vector<int> > > & xsim);

    //caulate the peak probabilty of each bin
    double CalcPr(vector< vector< vector<int> > > & zinput, vector<double> & a, int i, int j, int SampleID);

    //caculate the pvalue of shared peaks
    void PPModeCar(vector< vector< vector<int> > > & zinput, vector<double> & a, int i, int j);
			
    //constructor and destructor
    component(void){}
    ~component(){}
};

component & component::operator=(const component & rhs)
{
 loglikefull = rhs.loglikefull;
 theta = rhs.theta;
 phi = rhs.phi;
 IsingGamma = rhs.IsingGamma;
 IsingPsi = rhs.IsingPsi;

 data = rhs.data; 
 
 CountTheta = rhs.CountTheta;
 CountPhi = rhs.CountPhi;
 CountGamma = rhs.CountGamma; 
 CountPsi = rhs.CountPsi;

 TuneTheta = rhs.TuneTheta;
 TunePhi = rhs.TunePhi;
 TuneGamma = rhs.TuneGamma;
 TunePsi = rhs.TunePsi;
 
 return *this;
}
