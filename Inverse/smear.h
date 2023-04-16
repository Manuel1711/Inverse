//#include "pars.h"
#include "/Users/manuel/Documents/GitHub/Sphelareons/Programmi/pars_sphal.h"
#include "statistical.h"
 
 


// DEFINITIONS OF THE FUNCTIONS USED FOR THE RESOLUTION METHOD



//Target function
Real Z( Real s, Real Es){
  return (1+erf(Es/(sqrt(2)*s)))/2;
}
Real Target_F(Real E, Real Es, Real s){
  //return exp(-pow(E-Es,2)/(2*s*s))/(sqrt(2*Pi)*s*Z(s,Es));
  return 1/pow((s*Pi/2),2)*(E-Estar)/(sinh((E)/s)); 
}


//Basis function
Real K(Real omega, Real t, Real beta){
  Real ret;
#if defined(EXP)
  ret=exp(-omega*t); //+ exp(-(beta-t)*omega);
#endif 
#if defined(COS)
  Real A=cosh(omega*(t - beta/2));
  Real B=sinh(beta*omega/2);
  //cout << "B: " << beta << endl;
  ret = A/B*omega; //Se moltiplico per omega alla fine ottengo rho/omega che è quello che mi interessa 
#endif
#if defined(COS_SPHAL)
  Real A=cosh(omega*(t-beta/2));
  Real B=sinh(beta*omega/2);
  ret = A/B*omega;
#endif
  return ret; 
}


//Smearing function
Real Delta_Smear(Real omega, PrecVec q, Real t_in[]){  
  Real D;
  for(int i=0; i<Nt; i++){
    if(omega==0) D += q(i)*2/beta;
    else D += q(i)*K(omega, t_in[i], beta);
  }
  return D; 
}

//Target funct integral
Real Target_Int(double ap){
  const auto f_cs=
    [=](const Real& E) -> Real
    {
      return Target_F(E, 0, sigma)*Target_F(E, 0, sigma)*exp(ap*E);
    };
  return bq::gauss_kronrod<Real,61>::integrate(f_cs,infLimit,supLimit,5,1e-16);
}



//Coefficients computation
PrecVec Coeff(PrecVec R, PrecMatr Winv, PrecVec f, Real l){
  Real den =  R.transpose()*Winv*R;
#if defined(HLN)
  Real numA = R.transpose()*Winv*(1-l)*f;
  Real num = 1-numA;
  return Winv*(1-l)*f+ Winv*R*num/den;
#endif
#if defined(BG)
  return Winv*R/den;
#endif
}


//Spectral function from coefficients
Real spectral(PrecVec q, PrecVec C){
  Real rho=0;
  for(int i=0; i<Nt; i++){
    rho += q(i)*C(i);
  }
  return Pi*rho;
}


//Naive statistical uncertainty
Real stat_unc(PrecVec q, PrecVec dC){
  Real err=0;
  for(int i=0; i<Nt; i++){
    err += q(i)*dC(i);
    cout << "q: " << q(i) << "  "  << dC(i) << "  " << err << endl;
  }
  
  return Pi*err;
}




// USEFUL OUTPUT FUNCTIONS

//Output Smearing function
void Smear_output(char open_Delta_S[], double Estar, Real t_a[], PrecVec g){

  FILE *Delta_S;
  if ((Delta_S = fopen(open_Delta_S, "w")) == NULL ){
    printf("Error opening the input file: %s\n",open_Delta_S);
    exit(EXIT_FAILURE);
  }
  
  
#if defined(BG)
  E0=0;
#endif
  fprintf(Delta_S, "@type xy\n");
  for(double i=0; i<300; i++)
    fprintf(Delta_S, "%s " "%s\n", conv(E0 +0.0001 + i/100).c_str(), conv(Delta_Smear(E0 + 0.0001 + i/100, g, t_a)).c_str());
  
#if defined(HLN)
  
  fprintf(Delta_S, "\n \n @type xy \n");
  for(double i=0; i<300; i++){
    fprintf(Delta_S, "%s " "%s\n", conv(E0 + 0.0001 +i/100).c_str(), conv(Target_F(E0 + 0.0001 + i/100, Estar, sigma)).c_str());
  }
  fprintf(Delta_S, "\n \n @type xy \n");
  for(double i=0; i<300; i++){
    Real df=Target_F(E0 + 0.0001 + i/100, Estar, sigma) - Delta_Smear(E0 + 0.0001 + i/100, g, t_a);
    fprintf(Delta_S, "%s " "%s\n", conv(E0 + 0.0001 + i/100).c_str(), conv(df).c_str());
  }
  
  fclose(Delta_S);
  
  
#endif
  
}





//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Computation A, R, f (Independent from \lambda and correlators)
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


//Computation A
PrecMatr A_Comp(Real t_a[], double Estar, double ap){
  PrecMatr A(Nt,Nt);
  for(int i=0; i<Nt; i++){
    for(int j=0; j<Nt; j++){

      //With exp basis 
#if defined(EXP)
#if defined(HLN)
      A(i,j) = exp(-(ti+tj-alpha)*E0)/(ti+tj-alpha);
#endif
#if defined(BG)
      A(i,j) = -(-2+2*Estar*(ti+tj)-pow(Estar,2)*pow(ti+tj,2))/pow(ti+tj,3);
#endif      
#endif

      //With cos basis integrate
#if defined(COS) || (COS_SPHAL)
      const auto f_cs=
	[=](const Real& E) -> Real
	{
#if defined(BG)
	  return K(E, t_a[i], beta)*(E - Estar)*(E - Estar)*K(E,t_a[j],beta);
#endif
#if defined(HLN)
	  return K(E, t_a[i], beta)*K(E,t_a[j],beta)*exp(ap*E); 
#endif
	};
      A(i,j) = bq::gauss_kronrod<Real,61>::integrate(f_cs,infLimit,supLimit,5,1e-16);
#endif
      
    }//j
  }//i
  cout << "A: " << A << endl;
  return A;
} 

//Computation R
PrecVec R_Comp(Real t_a[]){
  PrecVec R(Nt);
  for(int i=0; i<Nt; i++){

#if defined(EXP)
#if defined(HLN)
    R(i) = 1/(t_a[i])*exp(-E0*t_a[i]);
#endif
#if defined(BG)
    R(i) = 1/(t_a[i]);
#endif
#endif

#if defined(COS) || (COS_SPHAL)
    const auto f_cs=
    [=](const Real& E) -> Real
    {
     return K(E, t_a[i], beta);
    };
    R(i) = bq::gauss_kronrod<Real,61>::integrate(f_cs,infLimit,supLimit,5,1e-16);
#endif
  }
  return R;
}


//Computation f
Real f_NInt(Real infLimit, Real supLimit, Real ti, Real s, Real Es, double ap){
  const auto f_cs=
    [=](const Real& E) -> Real
    {
     return K(E, ti, beta)*Target_F(E,Es,s)*exp(ap*E);
    };
  return
    bq::gauss_kronrod<Real,61>::integrate(f_cs,infLimit,supLimit,5,1e-16);
  
}
Real N(Real t, Real s, Real Es){
  return 1/(2*Z(s,Es))*exp(((alpha-t)*((alpha-t)*pow(s,2)+2*Es))/2);
}
Real D(Real t, Real s, Real Es){
  return  1+erf(((alpha-t)*pow(s,2)+Es-E0)/(sqrt(2)*s)); 
}
PrecVec f_func(Real t_a[], Real s, Real Es, double ap){
  PrecVec f(Nt);
  for(int i=0; i<Nt; i++){
#if defined(EXP)
    f(i) = N(t_a[i], s, Es)*D(t_a[i], s, Es);
#endif
#if defined(COS) || (COS_SPHAL)
    f(i) = f_NInt(infLimit, supLimit,t_a[i], s, Es, ap);
#endif
    //cout << "f: " << f(i) << endl;
  }
  return f;
}
 



 
 


//Computation W e coefficients
PrecVec g_comp(Real lambda, PrecMatr A, PrecMatr Cov, PrecVec f, PrecVec R, Real Corr0){
  PrecMatr W_Mat(Nt,Nt);
  W_Mat = (1-lambda)*A + lambda*Cov/(pow(Corr0,2));
  //Inversion matrix W
  const auto Winv=W_Mat.inverse();
  //Computation g
  return Coeff(R,Winv,f,lambda);
}


//Computation functional a posteriori (global coefficients)
Real W_func_comp(Real lambda, Real Corr_zero, PrecMatr Cov, double Estar, PrecVec g, double ap, PrecMatr A, PrecVec f){
  Real A1=g.transpose()*A*g;
  Real A2=-2*f.transpose()*g;
  Real A3=Target_Int(ap);
  Real B=lambda/(Corr_zero*Corr_zero)*g.transpose()*Cov*g;
  return (1-lambda)*(A1+A2+A3)+B;
}


void Residual_Study(Real Corr_Z, PrecMatr Cov, double Estar, PrecMatr gl, Real t[], PrecVec f, PrecMatr A){
  

  FILE *RS;
  char open_RS[1024];
  sprintf(open_RS, "Output/RS.out");
  if ((RS = fopen(open_RS, "w")) == NULL ){
    printf("Error opening the input file: %s\n",open_RS);
    exit(EXIT_FAILURE);
  } 
  
  int a=0;
  
  for(int ilambda=0; ilambda<Nlambda; ilambda++){
    fprintf(RS, "%s " "%s " "%s\n", conv(lambda(ilambda)).c_str(), conv(W_func_comp(0, Corr_Z, Cov, Estar, gl.row(ilambda), apar, A, f)/Target_Int(apar)).c_str(), conv(W_func_comp(1, Corr_Z, Cov, Estar, gl.row(ilambda), apar, A, f)).c_str());
    if(dLim*W_func_comp(0, Corr_Z, Cov, Estar, gl.row(ilambda), apar, A, f)/Target_Int(apar)>W_func_comp(1, Corr_Z, Cov, Estar, gl.row(ilambda), apar, A, f) and a==0){
      lambda1 = lambda(ilambda-1);
      ilambda1=ilambda;
      a++;
    }
    if(uLim*W_func_comp(0, Corr_Z, Cov, Estar, gl.row(ilambda), apar, A, f)/Target_Int(apar)>W_func_comp(1, Corr_Z, Cov, Estar, gl.row(ilambda), apar, A, f) and a==1){
      lambda2 = lambda(ilambda-1);
      ilambda2=ilambda;
      a++;
    }
    //cout << " ilambda : " << ilambda << endl;
  }

  cout << "PPP: " << lambda1 << "  " << lambda2 << endl;
 
  
  fclose(RS);
  
}

Real Sigma(Real Sigma_Stat2, Real Sigma_Stat1, Real M_1, Real M_2){
  cout << "BOOO " << endl;
  Real Ps = M_1-M_2/Sigma_Stat1;
  Real syst = abs(M_1-M_2)*erf(abs(Ps)/sqrt(2));
  return sqrt(Sigma_Stat2*Sigma_Stat2 + syst*syst);
  
}

void Print_Study_Lambda(PrecVec Dens_Mu, PrecVec Dens_S, Real Corr_Z, PrecMatr Cov, double Estar, PrecMatr gl, Real t[], Real SigmaF){ 
  
  
  //File di output e bin in lambda
  FILE *Lambda_Shape_out, *Lambda_A;
  char open_Lambda_Shape_out[1024], open_Lambda_A[1024];
  sprintf(open_Lambda_A, "Output/Lambda_A.out");
  //if(trash == 0) sprintf(open_Lambda_Shape_out, "Output/Lambda_Kotov.out");
  if(trash == 0) sprintf(open_Lambda_Shape_out, "Output/Lambda_Double.out");
  else sprintf(open_Lambda_Shape_out, "Output/Alpha/Ns%d_Nt%d/N_Cool_%d/Lambda_ShapeD%lfU%lf_%lf.out", S, beta, NCool, dLim, uLim, Estar);
  if ((Lambda_Shape_out = fopen(open_Lambda_Shape_out, "w")) == NULL ){
    printf("Error opening the input file: %s\n",open_Lambda_Shape_out);
    exit(EXIT_FAILURE);
    fprintf(Lambda_Shape_out, "@target G0.S0\n@type xydy\n");
  }
  if ((Lambda_A = fopen(open_Lambda_A, "w")) == NULL ){
    printf("Error opening the input file: %s\n",open_Lambda_A);
    exit(EXIT_FAILURE);
  } 
  
  //Output rho in funzione di lambda
  PrecMatr A0 = A_Comp(t, Estar, 0);
  PrecVec f0 = f_func(t, sigma, Estar, 0);
  fprintf(Lambda_Shape_out, "@type xydy \n");

  for(int ilambda=0; ilambda<Nlambda; ilambda++){


    cout << "OOO: " << lambda(ilambda) << "  " << sqrt(W_func_comp(0,Corr_Z, Cov, Estar, gl.row(ilambda), 0, A0, f0)/Target_Int(0)) << endl;
    
    fprintf(Lambda_Shape_out, "%s " "%s " "%s\n", conv(sqrt(W_func_comp(0,Corr_Z, Cov, Estar, gl.row(ilambda), 0, A0, f0)/Target_Int(0))).c_str(), conv(Dens_Mu(ilambda)).c_str(), conv(Dens_S(ilambda)).c_str());
    //fprintf(Lambda_A, "%s " "%s\n", conv(lambda).c_str(), conv(W_func_comp(0,Corr_zero,Cov,Estar,TSq)).c_str());
  }//lambda
  
  char eo;
  if(EO==0) eo='e';
  if(EO==1) eo='o';
  if(EO==2) eo='s';

  fprintf(Lambda_Shape_out, "\n & \n");
  fprintf(Lambda_Shape_out, "%s " "%s " "%s\n & \n", conv(sqrt(W_func_comp(0,Corr_Z, Cov, Estar, gl.row(ilambda1), 0, A0, f0)/Target_Int(0))).c_str(), conv(Dens_Mu(ilambda1)).c_str(), conv(Dens_S(ilambda1)).c_str());
  fprintf(Lambda_Shape_out, "%s " "%s " "%s\n & \n", conv(sqrt(W_func_comp(0,Corr_Z, Cov, Estar, gl.row(ilambda2), 0, A0, f0)/Target_Int(0))).c_str(), conv(Dens_Mu(ilambda2)).c_str(), conv(Dens_S(ilambda2)).c_str());
  fprintf(Lambda_Shape_out, "%s " "%s " "%s\n & \n", conv(sqrt(W_func_comp(0,Corr_Z, Cov, Estar, gl.row(ilambda2), 0, A0, f0)/Target_Int(0))).c_str(), conv(Dens_Mu(ilambda2)).c_str(), conv(SigmaF).c_str());

  //fprintf(Lambda_Shape_out, "%s  %s\n%s  %s\n", conv(Max_Lambda).c_str(), conv(-rho_m_max-rho_s_max).c_str(), conv(Max_Lambda).c_str(), conv(2*(rho_m_max+rho_s_max)).c_str());
  //fprintf(Lambda_Shape_out, "&\n@    s1 legend  \"\\xl\\f{}\\smax\\N=%s\"\n", conv(Max_Lambda).c_str());
  fprintf(Lambda_Shape_out, "@    xaxis  label \"d(g)\"\n@    yaxis  label \"-(\\xG\\f{}\\N\\s sphal\\N)/( T\\S4\\N)\"\n");
  fprintf(Lambda_Shape_out, "@version 50125\n@world xmin 0.005\n@world ymin %s\n@world xmax 0.15\n@world ymax %s", conv(Dens_Mu(0)-Dens_S(0)).c_str(), conv(Dens_Mu(0)+Dens_S(0)).c_str());
  fprintf(Lambda_Shape_out, "\n@    xaxis  tick major 1.5\n@    yaxis  tick major 0.0005\n");
  fprintf(Lambda_Shape_out, "@    s1 errorbar linewidth 3.0 \n@    s1 errorbar riser linewidth 3.0 \n@    s2 errorbar linewidth 3.0 \n@    s2 errorbar riser linewidth 3.0 \n@  s3 errorbar linewidth 2.0 \n@ s3 errorbar riser linewidth 2.0 \n@  s0 symbol 1 \n@    s0 line type 0 \n@    s1 symbol 2 \n@    s1 symbol fill pattern 1 \n@    s2 symbol 2 \n@    s2 symbol fill pattern 1 \n@    legend 0.65, 0.8 \n@    s2 legend  \"%f d\\s2-\\N[g**]\\S\\S2\\N = B[g**]\" \n@    s1 legend  \"%f d\\s2-\\N[g*]\\S2\\N = B[g*]\" \n@g0 hidden false \n@g0 type XY \n@g0 stacked false \n@g0 bar hgap 0.000000 \n@g0 fixedpoint off \n@g0 fixedpoint type 0 \n@g0 fixedpoint xy 0.000000, 0.000000 \n@g0 fixedpoint format general general \n@g0 fixedpoint prec 6, 6 \n@version 50122\n", uLim, dLim );
  
  
  fclose(Lambda_Shape_out);
  
}
