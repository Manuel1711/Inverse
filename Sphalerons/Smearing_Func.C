#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <cmath> 
#include "mp.h"
#include "smear.h"



int main(int argc, char *argv[]){
  parameters(argc, argv);
  cout << "********* " << Estar  << " *********" << endl;
  cout << "L: " << S << " #Points: " << Nt <<  " Nt: " << beta << " dLim: " << dLim << " uLim: " << uLim << " trash: " << trash << " sigma: " << sigma << " NCool: " << NCool << " Estar: " << Estar << endl;
  
  //Valori tau
  Real t[Nt];
  for(int i=1; i<=Nt; i++){
    //t[i-1] = 0.2 + (i-1)*0.05;
    t[i-1] = i;
    cout << "t: " << t[i-1] << endl;
  }
  
  //Valori lambda
  for(int ilambda=0; ilambda<Nlambda; ilambda++){
    if(ilambda < 20) lambda(ilambda) = conv(to_string(ilambda))/1000000000;
    else if(ilambda >= 20 and ilambda < 30) lambda(ilambda) = conv(to_string(ilambda-19))/100000000;
    else if(ilambda >= 30 and ilambda < 50) lambda(ilambda) = conv(to_string(ilambda-29))/10000000;
    else if(ilambda >= 50 and ilambda < 70) lambda(ilambda) = conv(to_string(ilambda-49))/1000000;
    else if(ilambda >= 70 and ilambda < 90) lambda(ilambda) = conv(to_string(ilambda-69))/100000;
    else if(ilambda >=90 and ilambda <110) lambda(ilambda) = conv(to_string(ilambda-89))/10000;
    else if(ilambda>=110 and ilambda<200) lambda(ilambda) = conv(to_string(ilambda-109))/1000;
    else if(ilambda>=200 and ilambda < 260) lambda(ilambda) = 0.1 + conv(to_string(ilambda-199))/100;
  }
    
  
  
  //Mean and Sigma Values
  FILE *Input_Corr_M;
  char open_Input_Corr_M[1024];
  if(trash==0) sprintf(open_Input_Corr_M, "Mean_Double.dat");
  else sprintf(open_Input_Corr_M, "bootstrap_tcorr_pureYM_L_%d_T_%d/Mean_topcharge_tcorr_ncool_%d.dat", S, beta, NCool);
  if ((Input_Corr_M = fopen(open_Input_Corr_M, "r")) == NULL ){
    printf("Error opening the input file: %s\n",open_Input_Corr_M);
    exit(EXIT_FAILURE);
  }
  
  PrecVec Corr_err(Nt+trash), Corr_Mu(Nt+trash);
  for(int i=0; i<Nt+trash; i++){
    double a;
    char App3[1024], App4[1024];
    fscanf(Input_Corr_M, "%lf" "%s" "%s", &a, App3, App4);
    Corr_Mu(i) = conv(App3);
    Corr_err(i) = conv(App4);
    cout << "KKK " << a << " " << Corr_Mu(i) << "  " <<  Corr_err(i) << endl;
  }
  fclose(Input_Corr_M);
  

  //Matrice di covarianza
  PrecMatr Cov(Nt,Nt);
  for(int i=0; i<Nt; i++){
    for(int j=0; j<Nt; j++){
      if(i==j){
	Cov(i,j)=Corr_err(i+trash)*Corr_err(i+trash);
      }
      else Cov(i,j)=0;
      cout << "Cov: " << Cov(i,j) << endl;
    }
  }
  
  
  //A,R,f
  PrecMatr A=A_Comp(t,Estar, apar);
  PrecVec R=R_Comp(t);
  PrecVec f;
#if defined(HLN)
  f = f_func(t, sigma, Estar, apar);
  cout << "f(" << Estar << ")= " << f << " f(0) = " << f_func(t,sigma,0,apar) << endl;
#endif
  
  //Coeff
  PrecMatr gl(Nlambda, Nt);
  for(int ilambda=0; ilambda<Nlambda; ilambda++){
    gl.row(ilambda) = g_comp(lambda(ilambda), A, Cov, f, R, Corr_Mu(0));
    cout << "gl: " << gl.row(ilambda) << endl;
  }
  
  //Smear
  char open_Delta_S[1024], open_Diff[1024];
  sprintf(open_Delta_S, "Output/Delta_Smear.out");
  Smear_output(open_Delta_S, Estar, t, gl.row(0));

  Real Ver;
  for(int i=0; i<Nt; i++){
    Ver += gl(0,i)*Corr_Mu(i+trash);
  }
  cout << "EEE: " << Ver << endl;
  
  
  //Calcolo integrale quadrato target per funzionale W
  Real TSq = Target_Int(apar);
  
  
  
  
  
  //BOOTSTRAP
  FILE *Input_Corrs;
  char open_Input_Corrs[1024];
  if(trash==0) sprintf(open_Input_Corrs, "final_corr_BOOTSTRAP_SAMPLES.dat");
  else sprintf(open_Input_Corrs, "bootstrap_tcorr_pureYM_L_%d_T_%d/bootstrap_topcharge_tcorr_ncool_%d.dat", S, beta, NCool);
  if ((Input_Corrs = fopen(open_Input_Corrs, "r")) == NULL ){
    printf("Error opening the input file: %s\n",open_Input_Corrs);
    exit(EXIT_FAILURE);
  } 
  
  PrecMatr Corr_Boot(Nboot, Nt+trash);
  for(int iboot=0; iboot<Nboot; iboot++){
    for(int i=0; i<Nt+trash; i++){
      
      int App1, App2;
      char App3[1024], App4[1024], App5[1024];
      if(Nboot==1) fscanf(Input_Corrs,  "%s" "%s" "%s", App3, App4, App5);
      else fscanf(Input_Corrs, "%d" "%s" "%s" "%s", &App1, App3, App4, App5);
      Corr_Boot(iboot,i) = conv(App4);
      cout << "iboot " << iboot << "  "  << i << "  " <<  Corr_Boot(iboot,i) <<  endl;
      
    }
  }
  

  //Calcolo Densità
  PrecMatr Dens(Nboot, Nlambda);
  cout << "HERE" << endl;
  for(int iboot=0; iboot<Nboot; iboot++){
    for(int ilambda=0; ilambda<Nlambda; ilambda++){
      for(int it=0; it<Nt; it++){
       	Dens(iboot, ilambda) += Corr_Boot(iboot, it+trash)*gl(ilambda,it);
      }
    }
  }
  PrecVec Dens_Mu(Nlambda), Dens_S(Nlambda);
  
  
  for(int ilambda=0; ilambda<Nlambda; ilambda++){
    Dens_Mu(ilambda) = 2*Pi*Boot_Mean(Dens.col(ilambda), Nboot)/(beta*beta);
    Dens_S(ilambda) = 2*Pi*Boot_Sigma(Dens.col(ilambda), Nboot)/(beta*beta);
  }
  
  Residual_Study(Corr_Mu(0), Cov, Estar, gl, t, f, A);
  
  Real SigmaF = Sigma(Dens_S(ilambda2), Dens_S(ilambda1), Dens_Mu(ilambda1), Dens_Mu(ilambda2));
  
  //Output
  ofstream file;
  char fout[1024];
  if(trash==0) sprintf(fout,"Output/Double.out");
  else sprintf(fout,"Output/Alpha/Ns%d_Nt%d/Output_ens.txt", S, beta);
  file.open(fout, std::ios::app);
  if (!file.is_open()) {
    std::cerr << "Errore: non è stato possibile aprire il file." << std::endl;
    return 1;
  } 
  file << NCool << "  " << Dens_Mu(ilambda2) << "  " << SigmaF << endl;
  //file << Estar << "  " << Dens_Mu(ilambda2) << "  " << SigmaF << endl;
  file.close();

  Print_Study_Lambda(Dens_Mu, Dens_S, Corr_Mu(0), Cov, Estar, gl, t, SigmaF);    
  

  
  
  return 0;
  
} 
 
