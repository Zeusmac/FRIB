#include <iostream>

double Gaus(double *dim, double *par)
{
  // 5 paremeters used in function: a, x0, sigma, slope, intercept
  double a = par[0];
  double x0 = par[1];
  double sigma = par[2];
  double bkgrd = par[3] * dim[0] + par[4];
  double x = dim[0];

  return a * TMath::Exp(-((x - x0) * (x - x0)) / (2 * sigma * sigma)) + bkgrd;
}

double LinFit(double *dim, double *par)
{
  return par[0] + par[1] * dim[0];
}


double decayE(double *dim, double *par)
{
  /* function for fitting decay curves 
  a is amplitude of source,
  b is tau of the source,
  c is tau of the daughter,
  bkgrd is backround,
  this first term is the exponetial decay the second term is batemans equation.
  */
  double a = par[0];
  double b = par[1];
  double c = par[2];
  double bkgrd = par[3];
  double x = dim[0];

  if( x < 0 ) return bkgrd;
  return a * TMath::Exp(-(b) * (x)) + ((b*a)/(c-b))*(TMath::Exp(-(b)*(x))-TMath::Exp(-(c)*(x)))+ bkgrd;
}

void Decayfit(TH1 *hist,double lower, double upper)
{
  //function to apply the fit DecayE to the root interface
  TF1 *decay = new TF1("decay", decayE, lower, upper,4);
  decay->SetParName(0, "source #");
  decay->SetParName(1, "Source lambda");
  decay->SetParName(2, "Daughter lambda");
  decay->SetParName(3, "bkground");
  decay->SetParameters(500, 5.2511150042*1e-8,8.0598509367*1e-9,1400);
  // decay->SetParLimits(0, 0, 100000);
  // decay->SetParLimits(1, 0, .1);
  //decay->SetParLimits(2, 8.0598509367*1e-9,8.0598509367*1e-9);
  // decay->SetParLimits(3, 0, 4800);
  // hist->GetListOfFunctions()->Add(decay);
  hist->Fit(decay, "+", "", lower, upper);
  Double_t decayhl = TMath::Log(2)/(decay ->GetParameter(1));
  Double_t decay_error = (TMath::Log(2)/TMath::Power((decay -> GetParameter(1)),2))*(decay -> GetParError(1));
  printf("Decay: %7fms +_ %7fms \n", decayhl*1e-6, decay_error*1e-6);
  Double_t decay_daughter = TMath::Log(2)/(decay -> GetParameter(2));
  Double_t decay_daughter_error = (TMath::Log(2)/TMath::Power((decay ->GetParameter(2)),2))*(decay -> GetParError(2));
  printf("Decay_daughter: %7fms +_ %7fms \n", decay_daughter*1e-6, decay_daughter_error*1e-6);
  gPad->Modified();
  gPad->Update();
  
}

void MultGausFit(TH1 *hist, int PeakNo){
     double lower,upper;
     for (int i = 0; i < PeakNo; i++) {
      cout<<"Write lower and upper" <<endl;
      cin>>lower >>upper;
      TF1 *F = new TF1("F", Gaus, lower, upper, 5);
      F->SetParName(0, "a");
      F->SetParName(1, "x0");
      F->SetParName(2, "sigma");
      F->SetParName(3, "slope");
      F->SetParName(4, "intercept");
      F->SetParameters(hist->GetMaximum(), ((lower+upper)/2), 2, 0.05, hist->GetMinimum());
      F->SetParLimits(0, 0, 1E6);
      F->SetParLimits(1, lower, upper);
      F->SetParLimits(4, 0, 1E6);
      //hist->GetListOfFunctions()->Add(F);
      hist->Fit(F, "+","", lower, upper);

      gPad->Modified();
      gPad->Update();
  }
}