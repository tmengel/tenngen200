#include "TennGenFunctions.h"
#include <bits/stdc++.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>

TENNGEN_BEGIN_NAMESPACE

Double_t Multiplicity(Double_t* X, Double_t* par){
    Double_t x      = X[0];
    return 0.0569104-0.00129485*x+1.45861e-05*x*x-8.83828e-08*TMath::Power(x,3)+3.10488e-10*TMath::Power(x,4)-6.54566e-13*TMath::Power(x,5)+8.25612e-16*TMath::Power(x,6)-5.9263e-19*TMath::Power(x,7)+2.10337e-22*TMath::Power(x,8)-2.40376e-26*TMath::Power(x,9);
}

Double_t getdNdptOverptIntegrand(Double_t* rad, Double_t* par){

	Double_t r 			= rad[0];
	// parameters used to fit: mass, beta, temp, n, norm
	Double_t mass 		= par[0];// not mT
	Double_t pt			= par[1];
	Double_t betaMax 	= par[2];
	Double_t temp 		= par[3];
	Double_t n 			= par[4];
	

	Double_t beta = betaMax*TMath::Power(r,n);
	if(beta > 0.99999999999999999999) beta = 0.99999999999999999999;

	Double_t mT 	= TMath::Sqrt(mass*mass + pt*pt);
	Double_t rho0 	= TMath::ATanH(beta);
	Double_t avoidFPE = pt*TMath::SinH(rho0)/temp;
	if(avoidFPE > 700.) avoidFPE = 700.;
	Double_t bk1arg = mT*TMath::CosH(rho0)/temp;
	Double_t integrand = r*mT*TMath::BesselI0(avoidFPE)*TMath::BesselK1(bk1arg);

	return integrand;
}

Double_t getdNdpt(Double_t* pT, Double_t* params){

	TF1* dNdptOverptIntegrandFunc = new TF1("integrandFunc", getdNdptOverptIntegrand, 0, 1, 5 );
	
	Double_t pt		= pT[0];
	Double_t mass 	= params[0];
	Double_t beta 	= params[1];
	Double_t temp 	= params[2];
	Double_t n 		= params[3];
	Double_t norm	= params[4];
	
	dNdptOverptIntegrandFunc->SetParameters(mass,pt,beta,temp,n);

	Double_t dNdptOverpt 	= dNdptOverptIntegrandFunc->Integral(0,1);
	
	Double_t dNdpt_normalized			= pt* norm * dNdptOverpt;
	gSystem->ProcessEvents();
	gROOT->Reset();
	return dNdpt_normalized;
}

float MASS[4][6] ={{ 0.13957 , 0.13957 , 0.49368 , 0.49368 , 0.93827 , 0.93827 }  , { 0.13957 , 0.13957 , 0.49368 , 0.49368 , 0.93827 , 0.93827 }  , { 0.13957 , 0.13957 , 0.49368 , 0.49368 , 0.93827 , 0.93827 }  , { 0.13957 , 0.13957 , 0.49368 , 0.49368 , 0.93827 , 0.93827 }  };
float BETA[4][6] ={{ 0.759913 , 0.760621 , 0.7319 , 0.732047 , 0.648344 , 0.680417 }  , { 0.765472 , 0.763144 , 0.739863 , 0.739893 , 0.641706 , 0.684773 }  , { 0.772956 , 0.773497 , 0.754914 , 0.756607 , 0.661272 , 0.692553 }  , { 0.913187 , 0.788666 , 0.780603 , 0.780076 , 0.681797 , 0.696818 }  };
float TEMP[4][6] ={{ 0.175254 , 0.17647 , 0.163823 , 0.167908 , 0.25008 , 0.256333 }  , { 0.175317 , 0.17444 , 0.166239 , 0.170082 , 0.247217 , 0.256199 }  , { 0.173122 , 0.17437 , 0.164637 , 0.168712 , 0.240833 , 0.249546 }  , { 0.0835104 , 0.169124 , 0.16093 , 0.163796 , 0.218115 , 0.22822 }  };
float NPAR[4][6] ={{ 32.371 , 33.1419 , 11.3702 , 12.1617 , 44.6847 , 74.0065 }  , { 30.7 , 28.7665 , 11.8083 , 12.3836 , 35.0342 , 65.2486 }  , { 28.408 , 29.0054 , 12.4268 , 13.1964 , 34.7224 , 54.8703 }  , { 27.5124 , 27.7284 , 14.806 , 14.6286 , 27.4966 , 36.6379 }  };
float NORM[4][6] ={{ 10211.8 , 9733.9 , 7669.71 , 6958.68 , 714.188 , 833.707 }  , { 6923.66 , 7024.51 , 4679.75 , 4275.71 , 521.372 , 559.871 }  , { 4030.29 , 3828.84 , 2667.37 , 2413.47 , 341.159 , 360.574 }  , { 109206 , 1663.87 , 1117.2 , 1030.07 , 254.529 , 241.505 }  };

std::string partcode(const int j){

	std::string particleID;

	if(j == 0) {particleID="pi-";}
	else if(j == 1) {particleID="pi+";}
	else if(j == 2) {particleID="ka-";}
	else if(j == 3) {particleID="ka+";}
	else if(j == 4) {particleID="pba";}
	else if(j == 5) {particleID="pro";}

	return particleID; 
}

void makeDistros200(std::string rootfilename){

    string outFile = rootfilename;
    TFile* outfile = new TFile(outFile.c_str(), "RECREATE");
	string PNGdir = "distro-PNG";
  if (mkdir(PNGdir.c_str(), 0777) == -1) std::cerr << "Error :  " << strerror(errno) << endl;

	int nCents = 4;
	int numparts = 6;
	TCanvas* c1;
	TF1* funcBGBW;
    TF1* multiplicity;

	for(int i =0; i<nCents; i++){
		for(int j =0; j< numparts; j++){

		c1 = new TCanvas(); // a la Rademakers
		gPad->SetLogy();
		gROOT-> SetBatch(kTRUE);// save canvases without displaying them
		c1->Update();


		funcBGBW = new TF1("funcBGBW",getdNdpt,0.00000000000001,10.,5); // actually has 5 parameters

		funcBGBW-> SetParameters(MASS[i][j],BETA[i][j],TEMP[i][j],NPAR[i][j],NORM[i][j]);
		funcBGBW->SetParNames("mass","beta (c)","temp","n","norm");

		string centrality = "cent"+std::to_string(i);
		string particleID = partcode(j);
		string histoFileName=particleID+"_pt_"+centrality;
		const char* funcName = histoFileName.c_str();

		TH1F* h1 = new TH1F(funcName,funcName, 1000, 0.,10.);

		h1->FillRandom("funcBGBW",1000000);
		h1->Scale(1.0/h1->Integral());

		TString xlabel = "p_{T}";
		TString ylabel = "#frac{d^{2}N}{dydp_{T}}";
		h1-> SetXTitle(xlabel);
		h1-> SetYTitle(ylabel);

		string imgPathAndName = "./distro-PNG/"+histoFileName+".png";
		TImage *png = TImage::Create();// FIXME try to use canvas method instead of png object
		png->FromPad(c1);
		const char* imgPathAndNameConstCharPtr = imgPathAndName.c_str();
		png->WriteImage(imgPathAndNameConstCharPtr);

		delete c1;
		delete funcBGBW;


		}

	}


    c1 = new TCanvas();
    string multName = "MultiplicityDistro";
    TH1F* h1 = new TH1F("MultiplicityDistro",multName.c_str(), 100, 0.,750.);
    multiplicity = new TF1("multiplicity",Multiplicity,0.0,750.,1);
    h1->FillRandom("multiplicity",1000000);
    h1->Scale(1.0/h1->Integral());
    
    string imgPathAndName = "./distro-PNG/"+multName+".png";
	TImage *png = TImage::Create();// FIXME try to use canvas method instead of png object
	png->FromPad(c1);
	const char* imgPathAndNameConstCharPtr = imgPathAndName.c_str();
	png->WriteImage(imgPathAndNameConstCharPtr);
	delete multiplicity;
    delete c1;
           
    outfile->Write();
    
}


TENNGEN_END_NAMESPACE
