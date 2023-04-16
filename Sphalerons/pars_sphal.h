//////////////////////////////

//In this header file, the quantities needed are defined
//The input file and the parameters are read and defined.

//Parameters
int Nt, S(0.), NCool(20.), beta, trash(0.);
double Estar;
double dLim(0.), uLim(0.), apar(0.);
Real sigma;



//Reading of the input file parameters and assignation
void parameters(int argc, char *argv[]){
  while( argc > 1 ) {
    
    switch(argv[1][0]) {
    case 'L':
      NCool = atoi( &argv[1][1] );
      break;
    case 'S':
      S = atoi( &argv[1][1] );
      break;
    case 'N':
      Nt = atoi( &argv[1][1] );
      break;
    case 'B':
      beta = atoi( &argv[1][1] );
      break;
    case 'O':
      Estar = atof( &argv[1][1] );
      break;
    case 'D':
      dLim = atof( &argv[1][1] );
      break;
    case 'U':
      uLim = atof( &argv[1][1] );
      break;   
    case 'T':
      trash = atoi( &argv[1][1] );
      break;
    case 'E':
      sigma = atof( &argv[1][1] );
      break;
    case 'A':
      apar = atof( &argv[1][1] );
      break;
    default:
      cerr << "Unlucky: Retry input values ens\n";
      exit (8);
    }
    ++argv;
    --argc;
  }
}



//Other useful quantities
static Real E0=0.0; 
static const Real alpha=0; 
//static const Real lambda=0.9999999;
static const int Nboot = 1000, Nlambda=260;
static const int Ncool = 1;
int EO = 2;
/////////////////////////////////////////                                 

string ens;

//Global quantities
PrecVec lambda(1000);
Real lambda1, lambda2;
int ilambda1, ilambda2;


//SETTA LA PRECISIONE DESIDERATA IN BITS                                  
const int P = 1024;
struct Initer
{ 
  Initer()
  {
    mpfr_class::set_dprec(P);
  }
};
  
Initer initer;
/////////////////////////////////////////

/////////////////////////////////////////                                   
//SETTA LIMITE INFERIORE E SUPERIORE INTEGRAZIONE NUMERICA
#if defined(BG)
const Real infLimit=0;
#endif
#if defined(HLN)
const Real infLimit=E0;
#endif
const Real supLimit=inf;			
/////////////////////////////////////////
 
 
