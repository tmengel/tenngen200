#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

#include "TRandom3.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TROOT.h"
#include "Riostream.h"
#include <TSystem.h>

using namespace std;


int getHarmonicCent(int multiplicity);
int getPhenixCent(int multiplicity);
int getCorrectedMultiplicity(int raw_multiplicity);
//float getMass(int K_F);
int PartInx(float mass);
Double_t getdNdptOverptIntegrand(Double_t* rad, Double_t* par);
Double_t getdNdpt(Double_t* pT, Double_t* params);
Double_t getdNdyIntegrand(Double_t* myPt, Double_t* par);
Double_t Multiplicity(Double_t* X, Double_t* par);
Double_t EtaFunction(Double_t* Eta, Double_t* params);
double RatioFunction(float pT, int cent, int type);
double OmegaFunction(float pT, int cent);
float HarmonicFunction(float pT, int cent, int type, int part);
float dNdPhi(float Phi, float V1, float V2, float V3, float V4, float PSI1, float PSI2, float PSI3, float PSI4);
int PartInd(int K_F);
string getPartString(int partInd);
int getKFpartType(const int partInd );

const int raw_high[5]={608,430,311,145,55};
const int raw_low[5]={431,312,146,56,0};
const int raw_bin[5]= {177,118,165,89,55};

const float correction_low[4] = {82,36,-10,-22};
const float correction_high[4]= {138,64,46,32};

const float piplusRatios[4] = {0.396738386,0.400412797,0.400958581,0.40230616};
const float piminusRatios[4] = {0.396738386,0.400412797,0.400958581,0.40230616};
const float kaplusRatios[4] = {0.063108127,0.061919505,0.060984334,0.059731231};
const float kaminusRatios[4] = {0.06118953,0.059236326,	0.058838257,0.057484421};
const float proRatios[4] = {0.043099904,0.041486068,0.042385006,0.042604604};
const float pbarRatios[4] = {0.03295875,0.032404541,0.033371486,0.034295646};

struct TGParticle{

        float pt,eta,phi;
        int KF;

        ~TGParticle() {};
        TGParticle(const float PT, const float ETA, const float PHI, const int K_F) :
                pt(PT), eta(ETA), phi(PHI), KF(K_F) {}
        float px =  pt*TMath::Cos(phi);
        float py =  pt*TMath::Sin(phi);
        float pz = pt*TMath::SinH(eta);
        



        float mass = getMass(KF) ;
        
        float p0 = TMath::Sqrt(pt*pt*TMath::CosH(eta)*TMath::CosH(eta) + mass*mass); 

        float Pt(){
                return pt;
        }
        float Eta(){
                return eta;
        }
        float Phi(){
                return phi;
        }
        float M(){
                return mass;
        }
        float Px(){
                return px;
        }
        float Py(){
                return py;
        }
        float Pz(){
                return pz;
        }
        float E(){
                return p0;
        }
        int Kf(){
                return KF;
        }
        
        float getMass(int kf);



};
float TGParticle::getMass(int K_F){
     float mass;
     if(K_F==211||K_F==-211){ mass= 0.13957;}
     else if(K_F==321||K_F==-321){mass= 0.49368;}
     else if(K_F==2212||K_F==-2212){mass=0.93827;}
     else {mass= 0.13957;} 
     return mass;
}

struct TGEvent{
        std::vector<TGParticle> particles;
        int centrality, multiplicity;
        ~TGEvent() {};
        TGEvent(std::vector<TGParticle> PARTS, const int CENT, const int MULT) :
                particles(PARTS), centrality(CENT), multiplicity(MULT) {}
        
        std::vector<TGParticle> parts(){
                return particles;
        }
        int cent(){
                return centrality;
        }
        int mult(){
                return multiplicity;
        }
};

std::vector<TGEvent> TannGen(const int nEvent, TRandom3* fRandom, const int setCent = 0, const float etaRange =1.1, const int doV1=1, const int doV2=1,const int doV3=1,const int doV4=1){
    

    TFile* myFile = new TFile("AuAu200Distros.root");  
    string centString;
    string partString;
    string HarmonicString = "_"+std::to_string(doV1)+std::to_string(doV2)+std::to_string(doV3)+std::to_string(doV4);
    if(setCent ==-1) centString = "allCents";
    else{ centString = "Cent"+std::to_string(setCent);}
    partString = getPartString(6);
    string tFileOUT="./TG_"+std::to_string(nEvent)+partString+"_eta"+std::to_string(etaRange)+"_"+centString+HarmonicString+".root";


    //TRandom3* fRandom= new TRandom3(3);
    
    int nRandCalls = 0;
    float eta_range = etaRange;

    float psi1,psi2,psi3,psi4;
    float ptPart, etaPart,phiPart, v1, v2, v3, v4, AbsVn, testValue, mass, dndphi,tempphi;
    double OMEGA, piplus,piminus,kaplus,kaminus,pbar,proton, Norm, partTEST;
    int rawmultiplicity, numparts, multiplicity, harmonicCent,spectraCent, K_F, partcode, phiCHECK, multiCheck,tmpMult, centCheck;
    string ptHistoPip,ptHistoPim,ptHistoKap,ptHistoKam,ptHistoPro,ptHistoPbar;
    int numpip,numpim,numkap,numkam,numpro,numpba;
    TH1F* ptDistroPip;
    TH1F* ptDistroPim;
    TH1F* ptDistroKap;
    TH1F* ptDistroKam;
    TH1F* ptDistroPro;
    TH1F* ptDistroPbar;
    TH1F* MultiDistro;

    TFile *ff = new TFile(tFileOUT.c_str(),"RECREATE"); 


    /*-----------------------------------------HISTOS-------------------------------------------*/ 

        char expression2[128];
        char expression3[128];
        char expression4[128];
        char expression5[128];
        char expression6[128];
        char expression7[128];
        char expression8[128];
        char expression12[128];
        char expression13[128];
        char expression14[128];
        char expression15[128];
        char expression16[128];
        char expression17[128];
        char expression18[128];
        char expression22[128];
        char expression23[128];
        char expression24[128];
        char expression25[128];
        char expression26[128];
        char expression27[128];
        char expression28[128];
        char expression32[128];
        char expression33[128];
        char expression34[128];
        char expression35[128];
        char expression36[128];
        char expression37[128];
        char expression38[128];
        char expression42[128];
        char expression43[128];
        char expression44[128];
        char expression45[128];
        char expression46[128];
        char expression47[128];
        char expression48[128];
        char expression52[128];
        char expression53[128];
        char expression54[128];
        char expression55[128];
        char expression56[128];
        char expression57[128];
        char expression58[128];
        char expression60[128];
        char expression61[128];

        sprintf(expression2 , "p_{T} Distribution for #pi^{+} cent 0 (0-10 )");
        sprintf(expression3 , "p_{T} Distribution for #pi^{-} cent 0 (0-10 )");
        sprintf(expression4 , "p_{T} Distribution for K^{+} cent 0 (0-10 )");
        sprintf(expression5 , "p_{T} Distribution for K^{-} cent 0 (0-10 )");
        sprintf(expression6 , "p_{T} Distribution for p^{+} cent 0 (0-10 )");
        sprintf(expression7 , "p_{T} Distribution for p^{-} cent 0 (0-10 )");
        sprintf(expression8 , "p_{T} Distribution for all particles cent 0 (0-10 )");
        sprintf(expression12 , "p_{T} Distribution for #pi^{+} cent 1 (10-20 )");
        sprintf(expression13 , "p_{T} Distribution for #pi^{-} cent 1 (10-20 )");
        sprintf(expression14 , "p_{T} Distribution for K^{+} cent 1 (10-20 )");
        sprintf(expression15 , "p_{T} Distribution for K^{-} cent 1 (10-20 )");
        sprintf(expression16 , "p_{T} Distribution for p^{+} cent 1 (10-20 )");
        sprintf(expression17 , "p_{T} Distribution for p^{-} cent 1 (10-20 )");
        sprintf(expression18 , "p_{T} Distribution for all particles cent 1 (10-20 )");
        sprintf(expression22 , "p_{T} Distribution for #pi^{+} cent 2 (40-60 )");
        sprintf(expression23 , "p_{T} Distribution for #pi^{-} cent 2 (40-60 )");
        sprintf(expression24 , "p_{T} Distribution for K^{+} cent 2 (40-60 )");
        sprintf(expression25 , "p_{T} Distribution for K^{-} cent 2 (40-60 )");
        sprintf(expression26 , "p_{T} Distribution for p^{+} cent 2 (40-60 )");
        sprintf(expression27 , "p_{T} Distribution for p^{-} cent 2 (40-60 )");
        sprintf(expression28 , "p_{T} Distribution for all particles cent 2 (40-60 )");
        sprintf(expression32 , "p_{T} Distribution for #pi^{+} cent 3 (40-60 )");
        sprintf(expression33 , "p_{T} Distribution for #pi^{-} cent 3 (40-60 )");
        sprintf(expression34 , "p_{T} Distribution for K^{+} cent 3 (40-60 )");
        sprintf(expression35 , "p_{T} Distribution for K^{-} cent 3 (40-60 )");
        sprintf(expression36 , "p_{T} Distribution for p^{+} cent 3 (40-60 )");
        sprintf(expression37 , "p_{T} Distribution for p^{-} cent 3 (40-60 )");
        sprintf(expression38 , "p_{T} Distribution for all particles cent 3 (40-60 )");
        sprintf(expression48 , "#eta Distribution for all particles");
        sprintf(expression58 , "#phi Distribution for all particles");
        sprintf(expression60 , "#phi vs #eta distribution of all particles");
        sprintf(expression61 , "#phi vs #eta distribution of all particles weighted by p_{T}");
        

        TH1D *histpT_piPlus_cent0 = new TH1D("histpT_piPlus_cent0", expression2,200,0,10); //for piPlus
        histpT_piPlus_cent0-> Sumw2();
        histpT_piPlus_cent0->SetXTitle("p_{T} (GeV/c)");
        histpT_piPlus_cent0->SetYTitle("dN/dp_{T}");
        
        TH1D *histpT_piMinus_cent0 = new TH1D("histpT_piMinus_cent0", expression3,200,0,10); //for piMinus
        histpT_piMinus_cent0-> Sumw2();
        histpT_piMinus_cent0->SetXTitle("p_{T} (GeV/c)");
        histpT_piMinus_cent0->SetYTitle("dN/dp_{T} ");
        
        TH1D *histpT_kPlus_cent0 = new TH1D("histpT_kPlus_cent0", expression4,200,0,10); //for kPlus
        histpT_kPlus_cent0-> Sumw2();
        histpT_kPlus_cent0->SetXTitle("p_{T} (GeV/c)");
        histpT_kPlus_cent0->SetYTitle("dN/dp_{T} ");
        
        TH1D *histpT_kMinus_cent0 = new TH1D("histpT_kMinus_cent0", expression5,200,0,10); //for kMinus
        histpT_kMinus_cent0-> Sumw2();
        histpT_kMinus_cent0->SetXTitle("p_{T} (GeV/c)");
        histpT_kMinus_cent0->SetYTitle("dN/dp_{T} ");
        
        TH1D *histpT_p_cent0 = new TH1D("histpT_p_cent0", expression6,200,0,10); //for proton
        histpT_p_cent0-> Sumw2();
        histpT_p_cent0->SetXTitle("p_{T} (GeV/c)");
        histpT_p_cent0->SetYTitle("dN/dp_{T} ");
        
        TH1D *histpT_pbar_cent0 = new TH1D("histpT_pbar_cent0", expression7,200,0,10); //for anti-proton
        histpT_pbar_cent0-> Sumw2();
        histpT_pbar_cent0->SetXTitle("p_{T} (GeV/c)");
        histpT_pbar_cent0->SetYTitle("dN/dp_{T} ");

        TH1D *histpT_all_cent0 = new TH1D("histpT_all_cent0", expression8,200,0,10); //for all particles (charged + neutral)
        histpT_all_cent0-> Sumw2();
        histpT_all_cent0->SetXTitle("p_{T} (GeV/c)");
        histpT_all_cent0->SetYTitle("dN/dp_{T} ");


        TH1D *histpT_piPlus_cent1 = new TH1D("histpT_piPlus_cent1", expression12,200,0,10); //for piPlus
        histpT_piPlus_cent1-> Sumw2();
        histpT_piPlus_cent1->SetXTitle("p_{T} (GeV/c)");
        histpT_piPlus_cent1->SetYTitle("dN/dp_{T}");
        
        TH1D *histpT_piMinus_cent1 = new TH1D("histpT_piMinus_cent1", expression13,200,0,10); //for piMinus
        histpT_piMinus_cent1-> Sumw2();
        histpT_piMinus_cent1->SetXTitle("p_{T} (GeV/c)");
        histpT_piMinus_cent1->SetYTitle("dN/dp_{T} ");
        
        TH1D *histpT_kPlus_cent1 = new TH1D("histpT_kPlus_cent1", expression14,200,0,10); //for kPlus
        histpT_kPlus_cent1-> Sumw2();
        histpT_kPlus_cent1->SetXTitle("p_{T} (GeV/c)");
        histpT_kPlus_cent1->SetYTitle("dN/dp_{T} ");
        
        TH1D *histpT_kMinus_cent1 = new TH1D("histpT_kMinus_cent1", expression15,200,0,10); //for kMinus
        histpT_kMinus_cent1-> Sumw2();
        histpT_kMinus_cent1->SetXTitle("p_{T} (GeV/c)");
        histpT_kMinus_cent1->SetYTitle("dN/dp_{T} ");
        
        TH1D *histpT_p_cent1 = new TH1D("histpT_p_cent1", expression16,200,0,10); //for proton
        histpT_p_cent1-> Sumw2();
        histpT_p_cent1->SetXTitle("p_{T} (GeV/c)");
        histpT_p_cent1->SetYTitle("dN/dp_{T} ");
        
        TH1D *histpT_pbar_cent1 = new TH1D("histpT_pbar_cent1", expression17,200,0,10); //for anti-proton
        histpT_pbar_cent1-> Sumw2();
        histpT_pbar_cent1->SetXTitle("p_{T} (GeV/c)");
        histpT_pbar_cent1->SetYTitle("dN/dp_{T} ");

        TH1D *histpT_all_cent1 = new TH1D("histpT_all_cent1", expression18,200,0,10); //for all particles (charged + neutral)
        histpT_all_cent1-> Sumw2();
        histpT_all_cent1->SetXTitle("p_{T} (GeV/c)");
        histpT_all_cent1->SetYTitle("dN/dp_{T} ");

        TH1D *histpT_piPlus_cent2 = new TH1D("histpT_piPlus_cent2", expression22,200,0,10); //for piPlus
        histpT_piPlus_cent2-> Sumw2();
        histpT_piPlus_cent2->SetXTitle("p_{T} (GeV/c)");
        histpT_piPlus_cent2->SetYTitle("dN/dp_{T}");
        
        TH1D *histpT_piMinus_cent2 = new TH1D("histpT_piMinus_cent2", expression23,200,0,10); //for piMinus
        histpT_piMinus_cent2-> Sumw2();
        histpT_piMinus_cent2->SetXTitle("p_{T} (GeV/c)");
        histpT_piMinus_cent2->SetYTitle("dN/dp_{T} ");
        
        TH1D *histpT_kPlus_cent2 = new TH1D("histpT_kPlus_cent2", expression24,200,0,10); //for kPlus
        histpT_kPlus_cent2-> Sumw2();
        histpT_kPlus_cent2->SetXTitle("p_{T} (GeV/c)");
        histpT_kPlus_cent2->SetYTitle("dN/dp_{T} ");
        
        TH1D *histpT_kMinus_cent2 = new TH1D("histpT_kMinus_cent2", expression25,200,0,10); //for kMinus
        histpT_kMinus_cent2-> Sumw2();
        histpT_kMinus_cent2->SetXTitle("p_{T} (GeV/c)");
        histpT_kMinus_cent2->SetYTitle("dN/dp_{T} ");
        
        TH1D *histpT_p_cent2 = new TH1D("histpT_p_cent2", expression26,200,0,10); //for proton
        histpT_p_cent2-> Sumw2();
        histpT_p_cent2->SetXTitle("p_{T} (GeV/c)");
        histpT_p_cent2->SetYTitle("dN/dp_{T} ");
        
        TH1D *histpT_pbar_cent2 = new TH1D("histpT_pbar_cent2", expression27,200,0,10); //for anti-proton
        histpT_pbar_cent2-> Sumw2();
        histpT_pbar_cent2->SetXTitle("p_{T} (GeV/c)");
        histpT_pbar_cent2->SetYTitle("dN/dp_{T} ");

        TH1D *histpT_all_cent2 = new TH1D("histpT_all_cent2", expression28,200,0,10); //for all particles (charged + neutral)
        histpT_all_cent2-> Sumw2();
        histpT_all_cent2->SetXTitle("p_{T} (GeV/c)");
        histpT_all_cent2->SetYTitle("dN/dp_{T} ");

        TH1D *histpT_piPlus_cent3 = new TH1D("histpT_piPlus_cent3", expression32,200,0,10); //for piPlus
        histpT_piPlus_cent3-> Sumw2();
        histpT_piPlus_cent3->SetXTitle("p_{T} (GeV/c)");
        histpT_piPlus_cent3->SetYTitle("dN/dp_{T}");
        
        TH1D *histpT_piMinus_cent3 = new TH1D("histpT_piMinus_cent3", expression33,200,0,10); //for piMinus
        histpT_piMinus_cent3-> Sumw2();
        histpT_piMinus_cent3->SetXTitle("p_{T} (GeV/c)");
        histpT_piMinus_cent3->SetYTitle("dN/dp_{T} ");
        
        TH1D *histpT_kPlus_cent3 = new TH1D("histpT_kPlus_cent3", expression34,200,0,10); //for kPlus
        histpT_kPlus_cent3-> Sumw2();
        histpT_kPlus_cent3->SetXTitle("p_{T} (GeV/c)");
        histpT_kPlus_cent3->SetYTitle("dN/dp_{T} ");
        
        TH1D *histpT_kMinus_cent3 = new TH1D("histpT_kMinus_cent3", expression35,200,0,10); //for kMinus
        histpT_kMinus_cent3-> Sumw2();
        histpT_kMinus_cent3->SetXTitle("p_{T} (GeV/c)");
        histpT_kMinus_cent3->SetYTitle("dN/dp_{T} ");
        
        TH1D *histpT_p_cent3 = new TH1D("histpT_p_cent3", expression36,200,0,10); //for proton
        histpT_p_cent3-> Sumw2();
        histpT_p_cent3->SetXTitle("p_{T} (GeV/c)");
        histpT_p_cent3->SetYTitle("dN/dp_{T} ");
        
        TH1D *histpT_pbar_cent3 = new TH1D("histpT_pbar_cent3", expression37,200,0,10); //for anti-proton
        histpT_pbar_cent3-> Sumw2();
        histpT_pbar_cent3->SetXTitle("p_{T} (GeV/c)");
        histpT_pbar_cent3->SetYTitle("dN/dp_{T} ");

        TH1D *histpT_all_cent3 = new TH1D("histpT_all_cent3", expression38,200,0,10); //for all particles (charged + neutral)
        histpT_all_cent3-> Sumw2();
        histpT_all_cent3->SetXTitle("p_{T} (GeV/c)");
        histpT_all_cent3->SetYTitle("dN/dp_{T} ");


        TH1D *histeta_all = new TH1D("histeta_all", expression48,200,-1.5,1.5); //for all charged particle
        histeta_all -> Sumw2();
        histeta_all->SetXTitle("eta");
        histeta_all->SetYTitle("dN/d#eta");


        TH1D *histphi_all = new TH1D("histphi_all", expression58,200,0,2*TMath::Pi()); //for all charged particles
        histphi_all -> Sumw2();
        histphi_all->SetXTitle("phi");
        histphi_all->SetYTitle("dN/d#phi");


        TH2D *histphi_eta_all = new TH2D("histphi_eta_all",expression60,200,0,2*TMath::Pi(),200,-1.1,1.1); //for all particles
        histphi_eta_all -> Sumw2();
        histphi_eta_all->SetXTitle("#phi (radians)");
        histphi_eta_all->SetYTitle("#eta (pseudo-rapidity)");
        histphi_eta_all->SetZTitle("dN^{ch.}/d#phid#eta");

        TH2D *hist_pT_phi_eta_all = new TH2D("hist_pT_phi_eta_all",expression61,200,0,2*TMath::Pi(),200,-1.1,1.1); //for all particles
        hist_pT_phi_eta_all -> Sumw2();
        hist_pT_phi_eta_all->SetXTitle("#phi (radians)");
        hist_pT_phi_eta_all->SetYTitle("#eta (pseudo-rapidity)");
        hist_pT_phi_eta_all->SetZTitle("dp_{T}^{ch.}/d#phid#eta");

    /*----------------------------------------END-HISTOS---------------------------------------*/
   
   
        MultiDistro = (TH1F*)myFile->Get("MultiplicityDistro");
        ptHistoPip = "pi+_pt_cent"+std::to_string(spectraCent);
        ptHistoPim = "pi-_pt_cent"+std::to_string(spectraCent);
        ptHistoKap = "ka+_pt_cent"+std::to_string(spectraCent);
        ptHistoKam = "ka-_pt_cent"+std::to_string(spectraCent);
        ptHistoPro = "pro_pt_cent"+std::to_string(spectraCent);
        ptHistoPbar = "pba_pt_cent"+std::to_string(spectraCent);

        ptDistroPip = (TH1F*)myFile->Get(ptHistoPip.c_str());
        ptDistroPim = (TH1F*)myFile->Get(ptHistoPim.c_str());
        ptDistroKap = (TH1F*)myFile->Get(ptHistoKap.c_str());
        ptDistroKam = (TH1F*)myFile->Get(ptHistoKam.c_str());
        ptDistroPro = (TH1F*)myFile->Get(ptHistoPro.c_str());
        ptDistroPbar = (TH1F*)myFile->Get(ptHistoPbar.c_str());


    std::vector<TGEvent> myEvents;
    
     
    cout << "STARTING TANNGEN: nEvents = " << nEvent << " Eta = " << eta_range << " Cent = " << setCent <<endl; 
 
    while(myEvents.size()!=nEvent){

        numparts = 0;
        std:vector<TGParticle> myparticles;
        
        psi1 = fRandom->Uniform(0,2*TMath::Pi());
        nRandCalls++;
        psi2 = 0;//fRandom->Uniform(0,2*TMath::Pi());
        psi3 = fRandom->Uniform(0,2*TMath::Pi());
        nRandCalls++;
        psi4 = 0; //fRandom->Uniform(0,2*TMath::Pi());

      
        tmpMult =0;
        multiCheck=0;
        while(multiCheck!=1){
                tmpMult = (int)MultiDistro->GetRandom(fRandom); 
                nRandCalls++;
                if(getPhenixCent(tmpMult) == setCent){
                        multiCheck=1;
                        rawmultiplicity = tmpMult; 
                }   
        }
                        
                
        harmonicCent = getHarmonicCent(rawmultiplicity);
        spectraCent = getPhenixCent(rawmultiplicity);
        
        float tmphigh = correction_high[spectraCent];
        float tmplow = correction_low[spectraCent];
        
        multiplicity = int(2.0*eta_range*(rawmultiplicity+int(fRandom->Uniform(tmplow,tmphigh))));
        nRandCalls++;

        numpip = int(piplusRatios[spectraCent]*multiplicity);
        numpim = int(piminusRatios[spectraCent]*multiplicity);
        numkap = int(kaplusRatios[spectraCent]*multiplicity);
        numkam = int(kaminusRatios[spectraCent]*multiplicity);
        numpro = int(proRatios[spectraCent]*multiplicity);
        numpba = int(pbarRatios[spectraCent]*multiplicity);
        nRandCalls+=2*multiplicity;

         for(int i=0; i<numpip;i++){

                ptPart =ptDistroPip->GetRandom(fRandom);

                etaPart = fRandom->Uniform(-eta_range,eta_range);

                K_F = 211;                        
                mass = 0.13957;
                partcode= PartInx(mass);

                v2 = doV2*HarmonicFunction(ptPart, harmonicCent,partcode, 0);
                v3 = doV3*HarmonicFunction(ptPart, harmonicCent,partcode, 1);
                v4 = doV4*HarmonicFunction(ptPart, harmonicCent,partcode, 2);
                v1 = doV1*0.02*v2;

                AbsVn = 2.0*(1.02*TMath::Abs(v2)+TMath::Abs(v3)+TMath::Abs(v4))/TMath::Pi();
                
                phiCHECK =0;
                while(phiCHECK!=1){
                        dndphi = AbsVn*fRandom->Uniform(0.0,1.0);
                        nRandCalls+=2;
                        tempphi = fRandom->Uniform(0,2.0*TMath::Pi());
                        if(dndphi < dNdPhi(tempphi,v1,v2,v3,v4,psi1,psi2,psi3,psi4)){
                                phiPart = tempphi;
                                phiCHECK = 1;
                        }
                }
                numparts++;

                histphi_all->Fill(phiPart); 
                histeta_all->Fill(etaPart); 
                histphi_eta_all->Fill(phiPart,etaPart);
                hist_pT_phi_eta_all->Fill(phiPart,etaPart,ptPart);

                if(spectraCent==0) histpT_all_cent0->Fill(ptPart);
                else if(spectraCent==1) histpT_all_cent1->Fill(ptPart);
                else if(spectraCent==2)  histpT_all_cent2->Fill(ptPart);
                else if(spectraCent ==3) histpT_all_cent3->Fill(ptPart);

                if(spectraCent == 0){ histpT_piPlus_cent0->Fill( ptPart );} 
                else if(spectraCent == 1){ histpT_piPlus_cent1->Fill( ptPart );} 
                else if(spectraCent == 2){ histpT_piPlus_cent2->Fill( ptPart );} 
                else if(spectraCent == 3){ histpT_piPlus_cent3->Fill( ptPart );} 

                myparticles.push_back(TGParticle(ptPart,etaPart,phiPart,K_F));

        }
        for(int i=0; i<numpim;i++){
                ptPart =ptDistroPim->GetRandom(fRandom);
                etaPart = fRandom->Uniform(-eta_range,eta_range);

                K_F = -211;                        
                mass = 0.13957;
                partcode= PartInx(mass);

                v2 = doV2*HarmonicFunction(ptPart, harmonicCent,partcode, 0);
                v3 = doV3*HarmonicFunction(ptPart, harmonicCent,partcode, 1);
                v4 = doV4*HarmonicFunction(ptPart, harmonicCent,partcode, 2);
                v1 = doV1*0.02*v2;

                AbsVn = 2.0*(1.02*TMath::Abs(v2)+TMath::Abs(v3)+TMath::Abs(v4))/TMath::Pi();
                
                phiCHECK =0;
                while(phiCHECK!=1){
                        dndphi = AbsVn*fRandom->Uniform(0.0,1.0);
                        tempphi = fRandom->Uniform(0,2.0*TMath::Pi());
                        nRandCalls+=2;
                        if(dndphi < dNdPhi(tempphi,v1,v2,v3,v4,psi1,psi2,psi3,psi4)){
                                phiPart = tempphi;
                                phiCHECK = 1;
                        }
                }
                numparts++;

                histphi_all->Fill(phiPart); 
                histeta_all->Fill(etaPart); 
                histphi_eta_all->Fill(phiPart,etaPart);
                hist_pT_phi_eta_all->Fill(phiPart,etaPart,ptPart);

                if(spectraCent==0) histpT_all_cent0->Fill(ptPart);
                else if(spectraCent==1) histpT_all_cent1->Fill(ptPart);
                else if(spectraCent==2)  histpT_all_cent2->Fill(ptPart);
                else if(spectraCent ==3) histpT_all_cent3->Fill(ptPart);

                if(spectraCent == 0){ histpT_piMinus_cent0->Fill( ptPart );} 
                else if(spectraCent == 1){ histpT_piMinus_cent1->Fill( ptPart );} 
                else if(spectraCent == 2){ histpT_piMinus_cent2->Fill( ptPart );} 
                else if(spectraCent == 3){ histpT_piMinus_cent3->Fill( ptPart );}

                myparticles.push_back(TGParticle(ptPart,etaPart,phiPart,K_F));
                
        }
        for(int i=0; i<numkap;i++){
                ptPart =ptDistroKap->GetRandom(fRandom);
                etaPart = fRandom->Uniform(-eta_range,eta_range);

                K_F = 321;                        
                mass = 0.49368;
                partcode= PartInx(mass);

                v2 = doV2*HarmonicFunction(ptPart, harmonicCent,partcode, 0);
                v3 = doV3*HarmonicFunction(ptPart, harmonicCent,partcode, 1);
                v4 = doV4*HarmonicFunction(ptPart, harmonicCent,partcode, 2);
                v1 = doV1*0.02*v2;

                AbsVn = 2.0*(1.02*TMath::Abs(v2)+TMath::Abs(v3)+TMath::Abs(v4))/TMath::Pi();
                
                phiCHECK =0;
                while(phiCHECK!=1){
                        dndphi = AbsVn*fRandom->Uniform(0.0,1.0);
                        tempphi = fRandom->Uniform(0,2.0*TMath::Pi());
                        nRandCalls+=2;
                        if(dndphi < dNdPhi(tempphi,v1,v2,v3,v4,psi1,psi2,psi3,psi4)){
                                phiPart = tempphi;
                                phiCHECK = 1;
                        }
                }
                numparts++;

                histphi_all->Fill(phiPart); 
                histeta_all->Fill(etaPart); 
                histphi_eta_all->Fill(phiPart,etaPart);
                hist_pT_phi_eta_all->Fill(phiPart,etaPart,ptPart);

                if(spectraCent==0) histpT_all_cent0->Fill(ptPart);
                else if(spectraCent==1) histpT_all_cent1->Fill(ptPart);
                else if(spectraCent==2)  histpT_all_cent2->Fill(ptPart);
                else if(spectraCent ==3) histpT_all_cent3->Fill(ptPart);

                if(spectraCent == 0){ histpT_kPlus_cent0->Fill( ptPart );}  
                else if(spectraCent == 1){ histpT_kPlus_cent1->Fill( ptPart );} 
                else if(spectraCent == 2){ histpT_kPlus_cent2->Fill( ptPart );} 
                else if(spectraCent == 3){ histpT_kPlus_cent3->Fill( ptPart );}
                
                myparticles.push_back(TGParticle(ptPart,etaPart,phiPart,K_F));
                
        }
        for(int i=0; i<numkam;i++){
                ptPart =ptDistroKam->GetRandom(fRandom);
                etaPart = fRandom->Uniform(-eta_range,eta_range);

                K_F = -321;                        
                mass = 0.49368;
                partcode= PartInx(mass);

                v2 = doV2*HarmonicFunction(ptPart, harmonicCent,partcode, 0);
                v3 = doV3*HarmonicFunction(ptPart, harmonicCent,partcode, 1);
                v4 = doV4*HarmonicFunction(ptPart, harmonicCent,partcode, 2);
                v1 = doV1*0.02*v2;

                AbsVn = 2.0*(1.02*TMath::Abs(v2)+TMath::Abs(v3)+TMath::Abs(v4))/TMath::Pi();
                
                phiCHECK =0;
                while(phiCHECK!=1){
                        dndphi = AbsVn*fRandom->Uniform(0.0,1.0);
                        tempphi = fRandom->Uniform(0,2.0*TMath::Pi());
                        nRandCalls+=2;
                        if(dndphi < dNdPhi(tempphi,v1,v2,v3,v4,psi1,psi2,psi3,psi4)){
                                phiPart = tempphi;
                                phiCHECK = 1;
                        }
                }
                numparts++;

                histphi_all->Fill(phiPart); 
                histeta_all->Fill(etaPart); 
                histphi_eta_all->Fill(phiPart,etaPart);
                hist_pT_phi_eta_all->Fill(phiPart,etaPart,ptPart);

                if(spectraCent==0) histpT_all_cent0->Fill(ptPart);
                else if(spectraCent==1) histpT_all_cent1->Fill(ptPart);
                else if(spectraCent==2)  histpT_all_cent2->Fill(ptPart);
                else if(spectraCent ==3) histpT_all_cent3->Fill(ptPart);

                if(spectraCent == 0){ histpT_kMinus_cent0->Fill( ptPart );}  
                else if(spectraCent == 1){ histpT_kMinus_cent1->Fill( ptPart );} 
                else if(spectraCent == 2){ histpT_kMinus_cent2->Fill( ptPart );} 
                else if(spectraCent == 3){ histpT_kMinus_cent3->Fill( ptPart );}
                
                myparticles.push_back(TGParticle(ptPart,etaPart,phiPart,K_F));
        }
        for(int i=0; i<numpro;i++){
                ptPart =ptDistroPro->GetRandom(fRandom);
                etaPart = fRandom->Uniform(-eta_range,eta_range);

                K_F = 2212;                        
                mass = 0.93827;
                partcode= PartInx(mass);

                v2 = doV2*HarmonicFunction(ptPart, harmonicCent,partcode, 0);
                v3 = doV3*HarmonicFunction(ptPart, harmonicCent,partcode, 1);
                v4 = doV4*HarmonicFunction(ptPart, harmonicCent,partcode, 2);
                v1 = doV1*0.02*v2;

                AbsVn = 2.0*(1.02*TMath::Abs(v2)+TMath::Abs(v3)+TMath::Abs(v4))/TMath::Pi();
                
                phiCHECK =0;
                while(phiCHECK!=1){
                        dndphi = AbsVn*fRandom->Uniform(0.0,1.0);
                        tempphi = fRandom->Uniform(0,2.0*TMath::Pi());
                        nRandCalls+=2;
                        if(dndphi < dNdPhi(tempphi,v1,v2,v3,v4,psi1,psi2,psi3,psi4)){
                                phiPart = tempphi;
                                phiCHECK = 1;
                        }
                }
                numparts++;

                histphi_all->Fill(phiPart); 
                histeta_all->Fill(etaPart); 
                histphi_eta_all->Fill(phiPart,etaPart);
                hist_pT_phi_eta_all->Fill(phiPart,etaPart,ptPart);

                if(spectraCent==0) histpT_all_cent0->Fill(ptPart);
                else if(spectraCent==1) histpT_all_cent1->Fill(ptPart);
                else if(spectraCent==2)  histpT_all_cent2->Fill(ptPart);
                else if(spectraCent ==3) histpT_all_cent3->Fill(ptPart);

                if(spectraCent == 0){ histpT_p_cent0->Fill( ptPart );}  
                else if(spectraCent == 1){ histpT_p_cent1->Fill( ptPart );} 
                else if(spectraCent == 2){ histpT_p_cent2->Fill( ptPart );} 
                else if(spectraCent == 3){ histpT_p_cent3->Fill( ptPart );}
                
                myparticles.push_back(TGParticle(ptPart,etaPart,phiPart,K_F));
                
        }
        for(int i=0; i<numpba;i++){

                ptPart =ptDistroPbar->GetRandom(fRandom);
                etaPart = fRandom->Uniform(-eta_range,eta_range);

                K_F = -2212;                        
                mass = 0.93827;
                partcode= PartInx(mass);

                v2 = doV2*HarmonicFunction(ptPart, harmonicCent,partcode, 0);
                v3 = doV3*HarmonicFunction(ptPart, harmonicCent,partcode, 1);
                v4 = doV4*HarmonicFunction(ptPart, harmonicCent,partcode, 2);
                v1 = doV1*0.02*v2;

                AbsVn = 2.0*(1.02*TMath::Abs(v2)+TMath::Abs(v3)+TMath::Abs(v4))/TMath::Pi();
                
                phiCHECK =0;
                while(phiCHECK!=1){
                        dndphi = AbsVn*fRandom->Uniform(0.0,1.0);
                        tempphi = fRandom->Uniform(0,2.0*TMath::Pi());
                        nRandCalls+=2;
                        if(dndphi < dNdPhi(tempphi,v1,v2,v3,v4,psi1,psi2,psi3,psi4)){
                                phiPart = tempphi;
                                phiCHECK = 1;
                        }
                }
                numparts++;

                histphi_all->Fill(phiPart); 
                histeta_all->Fill(etaPart); 
                histphi_eta_all->Fill(phiPart,etaPart);
                hist_pT_phi_eta_all->Fill(phiPart,etaPart,ptPart);

                if(spectraCent==0) histpT_all_cent0->Fill(ptPart);
                else if(spectraCent==1) histpT_all_cent1->Fill(ptPart);
                else if(spectraCent==2)  histpT_all_cent2->Fill(ptPart);
                else if(spectraCent ==3) histpT_all_cent3->Fill(ptPart);

                if(spectraCent == 0){ histpT_pbar_cent0->Fill( ptPart );}  
                else if(spectraCent == 1){ histpT_pbar_cent1->Fill( ptPart );} 
                else if(spectraCent == 2){ histpT_pbar_cent2->Fill( ptPart );} 
                else if(spectraCent == 3){ histpT_pbar_cent3->Fill( ptPart );}
                
                myparticles.push_back(TGParticle(ptPart,etaPart,phiPart,K_F));
        
                
        }
        myEvents.push_back(TGEvent(myparticles, spectraCent, multiplicity));
        gSystem->ProcessEvents();
        gROOT->Reset();
        if(myEvents.size()%(nEvent-1)==0){
                cout << "event #: " << myEvents.size() << " complete" <<endl; 
        }

    }
        
   
    //nRandCalls = nRandCalls/nEvent;
    //cout << "Ending: " << "RandCalls/Event: " << nRandCalls << endl;
    
   
    delete MultiDistro; 
    delete ptDistroPip,ptDistroPim,ptDistroKap,ptDistroKam,ptDistroPro,ptDistroPbar; 
    //delete fRandom;
    delete myFile;
    

    TDirectory *all_particles_etaphi = ff->mkdir("all_particles_etaphi");
    
    if(setCent ==0){

        TDirectory *pt_cent0 = ff->mkdir("pt_cent0");
        pt_cent0->cd();

        histpT_piPlus_cent0->Scale(1.0/histpT_piPlus_cent0->Integral());
        histpT_piMinus_cent0->Scale(1.0/histpT_piMinus_cent0->Integral());
        histpT_kPlus_cent0->Scale(1.0/histpT_kPlus_cent0->Integral());
        histpT_kMinus_cent0->Scale(1.0/histpT_kMinus_cent0->Integral());
        histpT_p_cent0->Scale(1.0/histpT_p_cent0->Integral());
        histpT_pbar_cent0->Scale(1.0/histpT_pbar_cent0->Integral());
        histpT_all_cent0->Scale(1.0/histpT_all_cent0->Integral());

        histpT_piPlus_cent0->Write();
        histpT_piMinus_cent0->Write();
        histpT_kPlus_cent0->Write();
        histpT_kMinus_cent0->Write();
        histpT_p_cent0->Write();
        histpT_pbar_cent0->Write(); 
        histpT_all_cent0->Write();

    }
    if(setCent ==1){
        TDirectory *pt_cent1 = ff->mkdir("pt_cent1");
        pt_cent1->cd();

        histpT_piPlus_cent1->Scale(1.0/histpT_piPlus_cent1->Integral());
        histpT_piMinus_cent1->Scale(1.0/histpT_piMinus_cent1->Integral());
        histpT_kPlus_cent1->Scale(1.0/histpT_kPlus_cent1->Integral());
        histpT_kMinus_cent1->Scale(1.0/histpT_kMinus_cent1->Integral());
        histpT_p_cent1->Scale(1.0/histpT_p_cent1->Integral());
        histpT_pbar_cent1->Scale(1.0/histpT_pbar_cent1->Integral());
        histpT_all_cent1->Scale(1.0/histpT_all_cent1->Integral());

        histpT_piPlus_cent1->Write();
        histpT_piMinus_cent1->Write();
        histpT_kPlus_cent1->Write();
        histpT_kMinus_cent1->Write();
        histpT_p_cent1->Write();
        histpT_pbar_cent1->Write(); 
        histpT_all_cent1->Write();
    }
    if(setCent==2){
        TDirectory *pt_cent2 = ff->mkdir("pt_cent2");
        pt_cent2->cd();

        histpT_piPlus_cent2->Scale(1.0/histpT_piPlus_cent2->Integral());
        histpT_piMinus_cent2->Scale(1.0/histpT_piMinus_cent2->Integral());
        histpT_kPlus_cent2->Scale(1.0/histpT_kPlus_cent2->Integral());
        histpT_kMinus_cent2->Scale(1.0/histpT_kMinus_cent2->Integral());
        histpT_p_cent2->Scale(1.0/histpT_p_cent2->Integral());
        histpT_pbar_cent2->Scale(1.0/histpT_pbar_cent2->Integral());
        histpT_all_cent2->Scale(1.0/histpT_all_cent2->Integral());

        histpT_piPlus_cent2->Write();
        histpT_piMinus_cent2->Write();
        histpT_kPlus_cent2->Write();
        histpT_kMinus_cent2->Write();
        histpT_p_cent2->Write();
        histpT_pbar_cent2->Write(); 
        histpT_all_cent2->Write();
    }
    if(setCent==3){
        TDirectory *pt_cent3 = ff->mkdir("pt_cent3");
        pt_cent3->cd();

        histpT_piPlus_cent3->Scale(1.0/histpT_piPlus_cent3->Integral());
        histpT_piMinus_cent3->Scale(1.0/histpT_piMinus_cent3->Integral());
        histpT_kPlus_cent3->Scale(1.0/histpT_kPlus_cent3->Integral());
        histpT_kMinus_cent3->Scale(1.0/histpT_kMinus_cent3->Integral());
        histpT_p_cent3->Scale(1.0/histpT_p_cent3->Integral());
        histpT_pbar_cent3->Scale(1.0/histpT_pbar_cent3->Integral());
        histpT_all_cent3->Scale(1.0/histpT_all_cent3->Integral());

        histpT_piPlus_cent3->Write();
        histpT_piMinus_cent3->Write();
        histpT_kPlus_cent3->Write();
        histpT_kMinus_cent3->Write();
        histpT_p_cent3->Write();
        histpT_pbar_cent3->Write(); 
        histpT_all_cent3->Write();
    }

    all_particles_etaphi->cd();

    histeta_all->Scale(1.0/histeta_all->Integral());
    hist_pT_phi_eta_all->Scale(1.0/hist_pT_phi_eta_all->Integral());
    histphi_eta_all->Scale(1.0/histphi_eta_all->Integral());
    histphi_all->Scale(1.0/histphi_all->Integral());

    histphi_eta_all->Write();
    hist_pT_phi_eta_all->Write();
    histeta_all->Write();
    histphi_all->Write();
    ff->Close();


    return myEvents; 
    
}

int getHarmonicCent(int multiplicity){
        if(multiplicity>=431){return 0;}
	else if(multiplicity>=312){return 1;}
	else if(multiplicity>=217){return 2;}
	else if(multiplicity>=146){return 3;}
	else if(multiplicity>=94){return 4;}
	else if(multiplicity>=56){return 5;}
	else{return 6;}
}
int getPhenixCent(int multiplicity){
    if(multiplicity>=431){return 0;}
    else if(multiplicity>=312){return 1;}
    else if(multiplicity>=146){return 2;}
    else if(multiplicity>=56){return 3;}
    else{return 6;}
}
// float getMass(int K_F){
//      float mass;
//      if(K_F==211||K_F==-211){ mass= 0.13957;}
//      else if(K_F==321||K_F==-321){mass= 0.49368;}
//      else if(K_F==2212||K_F==-2212){mass=0.93827;}
//      else {mass= 0.13957;} 
//      return mass;
// }
int PartInx(float mass){
    if(mass==0.13957){return 0;}
    else if(mass==0.49368){return 1;}
    else if(mass==0.93827){return 2;}
    else{return 0;}
}
int getCorrectedMultiplicity(int raw_multiplicity){
		double correct_multiplicity;

		int raw_high_bins[7]={608,430,311,216,145,93,55};
		int raw_low_bins[7]={431,312,217,146,94,56,0};
		int raw_bin_width[7]= {177,118,94,70,51,37,55};



		int corrected_high_bins[7]={775,509,332,241,148,103,52};
		int corrected_low_bins[7]={510,333,242,149,104,53,0};
		int corrected_bin_width[7]= {265,176,90,92,44,52};


		int multiplicity_bin = getHarmonicCent(raw_multiplicity);
		double binFrac = (raw_multiplicity-raw_low_bins[multiplicity_bin])/raw_bin_width[multiplicity_bin];

		correct_multiplicity = int(corrected_low_bins[multiplicity_bin] + binFrac*corrected_bin_width[multiplicity_bin]);

		return correct_multiplicity;
}
Double_t getdNdptOverptIntegrand(Double_t* rad, Double_t* par){
	// (dN/dpt)/pt= r*dr*mt*I0((pt*sinh(rho))/T)*K1((mt*cosh(rho))/T)
	// rho=arctanh(beta); beta=betaMax*(r/R)^n
	//Double_t pT		= pt[0];
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
Double_t getdNdyIntegrand(Double_t* myPt, Double_t* par){
	Double_t pt   = myPt[0]; // x-axis of integration
	Double_t mass = par[0];
	
	return getdNdpt(myPt,par);
}
Double_t Multiplicity(Double_t* X, Double_t* par){
    Double_t x      = X[0];
    return 0.0569104-0.00129485*x+1.45861e-05*x*x-8.83828e-08*TMath::Power(x,3)+3.10488e-10*TMath::Power(x,4)-6.54566e-13*TMath::Power(x,5)+8.25612e-16*TMath::Power(x,6)-5.9263e-19*TMath::Power(x,7)+2.10337e-22*TMath::Power(x,8)-2.40376e-26*TMath::Power(x,9);
}
Double_t EtaFunction(Double_t* Eta, Double_t* params){
    Double_t eta    = Eta[0];
    Double_t A      = params[0]; //centraility
    if(A==0){  
        return 643.536-0.770784*eta+11.6018*eta*eta+0.458893*TMath::Power(eta,3)-4.61068*TMath::Power(eta,4)-0.0315425*TMath::Power(eta,5)+0.194857*TMath::Power(eta,6)+0.000584403*TMath::Power(eta,7)-0.00253618*TMath::Power(eta,8);
    }
    else if(A==1){
        return 496.741-1.37939*eta+8.1113*eta*eta+0.554655*TMath::Power(eta,3)-3.332*TMath::Power(eta,4)-0.0385251*TMath::Power(eta,5)+0.138938*TMath::Power(eta,6)+0.000740975*TMath::Power(eta,7)-0.00178451*TMath::Power(eta,8);
    }
    else if(A==2){
        return 348.744-1.43291*eta+4.94739*eta*eta+0.418635*TMath::Power(eta,3)-2.0966*TMath::Power(eta,4)-0.0278022*TMath::Power(eta,5)+0.0833634*TMath::Power(eta,6)+0.000535907*TMath::Power(eta,7)-0.00101071*TMath::Power(eta,8);
    }
    else if(A==3){  
        return 234.335-1.60194*eta+3.23287*eta*eta+0.433013*TMath::Power(eta,3)-1.3691*TMath::Power(eta,4)-0.0275018*TMath::Power(eta,5)+0.0549573*TMath::Power(eta,6)+0.00050864*TMath::Power(eta,7)-0.000677341*TMath::Power(eta,8);
    }
    else if(A==4){  
        return 152.768-0.838729*eta+1.41065*eta*eta+0.243627*TMath::Power(eta,3)-0.775913*TMath::Power(eta,4)-0.01596*TMath::Power(eta,5)+0.0306733*TMath::Power(eta,6)+0.000297723*TMath::Power(eta,7)-0.000370548*TMath::Power(eta,8);
    }
    else if(A==5){  
        return 91.7547-0.435019*eta+1.00857*eta*eta+0.134504*TMath::Power(eta,3)-0.457398*TMath::Power(eta,4)-0.00935008*TMath::Power(eta,5)+0.0175202*TMath::Power(eta,6)+0.00018388*TMath::Power(eta,7)-0.000205433*TMath::Power(eta,8);
    }
    else{return 0.0;}
}
double RatioFunction(float pT, int cent, int type){
        float pt     = pT;
        int A      = cent; //centraility
        int B      = type; // ratiotype
        if(B ==1 || B ==3  || B ==4){
                if(pt > 4.0){ pt = 4.0;}
        }
        if(A==0){
                if(B==0){ return 1.05232*TMath::Power(pt,0.0444398)*TMath::Exp(-0.0391436*pt);}
                else if(B==1){ return 0.973497*TMath::Power(pt,0.00746373)*TMath::Exp(-0.0231327*pt);}
                else if(B==2){ return 0.714843*TMath::Power(pt,-0.0765232)*TMath::Exp(0.0302362*pt);}
                else if(B==3){ return 0.368629*TMath::Power(pt,0.943152)*TMath::Exp(-0.172881*pt);}
                else if(B==4){ return 0.378045*TMath::Power(pt,0.902791)*TMath::Exp(-0.139721*pt);}
                else if(B==5){ return 0.430421*TMath::Power(pt,3.51419)*TMath::Exp(-1.18616*pt);}
                else if(B==6){ return 0.623285*TMath::Power(pt,3.61225)*TMath::Exp(-1.24021*pt);}
        }
        else if(A==1){
                if(B==0){ return 0.999801*TMath::Power(pt,-0.0249828)*TMath::Exp(0.0081137*pt);}
                else if(B==1){ return 0.963395*TMath::Power(pt,-0.0157204)*TMath::Exp(-0.0115677*pt);}
                else if(B==2){ return 0.746027*TMath::Power(pt,-0.0585174)*TMath::Exp(0.0057852*pt);}
                else if(B==3){ return 0.357072*TMath::Power(pt,0.933124)*TMath::Exp(-0.164785*pt);}
                else if(B==4){ return 0.369955*TMath::Power(pt,0.918328)*TMath::Exp(-0.143286*pt);}
                else if(B==5){ return 0.418842*TMath::Power(pt,3.43009)*TMath::Exp(-1.16914*pt);}
                else if(B==6){ return 0.571988*TMath::Power(pt,3.49185)*TMath::Exp(-1.18507*pt);}
        }
        else if(A==2){
                if(B==0){ return 1.03665*TMath::Power(pt,0.0206189)*TMath::Exp(-0.0233864*pt);}
                else if(B==1){ return 0.979407*TMath::Power(pt,0.00336944)*TMath::Exp(-0.0266132*pt);}
                else if(B==2){ return 0.766358*TMath::Power(pt,-0.0469923)*TMath::Exp(-0.00814672*pt);}
                else if(B==3){ return 0.341224*TMath::Power(pt,0.885281)*TMath::Exp(-0.144538*pt);}
                else if(B==4){ return 0.346016*TMath::Power(pt,0.842213)*TMath::Exp(-0.101821*pt);}
                else if(B==5){ return 0.392903*TMath::Power(pt,3.21868)*TMath::Exp(-1.10371*pt);}
                else if(B==6){ return 0.536303*TMath::Power(pt,3.29701)*TMath::Exp(-1.12682*pt);}
        }
        else if(A==3){
                if(B==0){ return 1.00893*TMath::Power(pt,-0.0147369)*TMath::Exp(0.00145578*pt);}
                else if(B==1){ return 0.961856*TMath::Power(pt,-0.0420376)*TMath::Exp(-0.00230289*pt);}
                else if(B==2){ return 0.777819*TMath::Power(pt,-0.0822537)*TMath::Exp(0.000508055*pt);}
                else if(B==3){ return 0.309516*TMath::Power(pt,0.775202)*TMath::Exp(-0.0930582*pt);}
                else if(B==4){ return 0.319566*TMath::Power(pt,0.776558)*TMath::Exp(-0.0740102*pt);}
                else if(B==5){ return 0.395983*TMath::Power(pt,2.91526)*TMath::Exp(-1.06911*pt);}
                else if(B==6){ return 0.517985*TMath::Power(pt,2.98984)*TMath::Exp(-1.07494*pt);}
        }
        //return 0.0;
}
double OmegaFunction(float pT, int cent){
        float pt     = pT;
        int A   = cent; //centrality
        if(A==0){
         return 1+0.368629*TMath::Power(pt,0.943152)*TMath::Exp(-0.172881*pt)+0.430421*TMath::Power(pt,3.51419)*TMath::Exp(-1.18616*pt)+(1.05232*TMath::Power(pt,0.0444398)*TMath::Exp(-0.0391436*pt))*(1+0.378045*TMath::Power(pt,0.902791)*TMath::Exp(-0.139721*pt)+0.623285*TMath::Power(pt,3.61225)*TMath::Exp(-1.24021*pt));
        }
        else if(A==1){
            return 1+0.357072*TMath::Power(pt,0.933124)*TMath::Exp(-0.164785*pt)+0.418842*TMath::Power(pt,3.43009)*TMath::Exp(-1.16914*pt)+(0.999801*TMath::Power(pt,-0.0249828)*TMath::Exp(0.0081137*pt))*(1+0.369955*TMath::Power(pt,0.918328)*TMath::Exp(-0.143286*pt)+0.571988*TMath::Power(pt,3.49185)*TMath::Exp(-1.18507*pt));
        }
        else if(A==2){
            return 1+0.341224*TMath::Power(pt,0.885281)*TMath::Exp(-0.144538*pt)+0.392903*TMath::Power(pt,3.21868)*TMath::Exp(-1.10371*pt)+(1.03665*TMath::Power(pt,0.0206189)*TMath::Exp(-0.0233864*pt))*(1+0.346016*TMath::Power(pt,0.842213)*TMath::Exp(-0.101821*pt)+0.536303*TMath::Power(pt,3.29701)*TMath::Exp(-1.12682*pt));
        }
        else if(A==3){
            return 1+0.309516*TMath::Power(pt,0.775202)*TMath::Exp(-0.0930582*pt)+0.395983*TMath::Power(pt,2.91526)*TMath::Exp(-1.06911*pt)+(1.00893*TMath::Power(pt,-0.0147369)*TMath::Exp(0.00145578*pt))*(1+0.319566*TMath::Power(pt,0.776558)*TMath::Exp(-0.0740102*pt)+0.517985*TMath::Power(pt,2.98984)*TMath::Exp(-1.07494*pt));
        }
       // return 0.0;
}
float HarmonicFunction(float pT, int cent, int type, int part){
        float pt;
        if(pT > 4.0){  pt = 4.0;  }
        else{ pt = pT; }
        int A   = cent;
        int B   = part;
        int C = type;
        if(pt < 0.5){ return 0.0; }
        else{
                if(A==0){
                        if(B==0){
                                if(C==0){
                                        return -0.0143392+0.0759773*pt-0.0154024*pt*pt-0.000915513*pt*pt*pt+0.000264875*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                        return -0.00558048+0.0136122*pt+0.0356842*pt*pt-0.0164884*pt*pt*pt+0.00195171*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return -0.00256332-0.00483424*pt+0.0390473*pt*pt-0.016438*pt*pt*pt+0.00205353*pt*pt*pt*pt;
                                }
                        }
                        else if(B==1){
                                if(C==0){
                                        return -0.0287711+0.0777838*pt-0.0110908*pt*pt-0.00251265*pt*pt*pt+0.000435634*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                        return -0.00843333+0.00413273*pt+0.03708*pt*pt-0.0140296*pt*pt*pt+0.00138853*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return -0.0147847+0.0238129*pt+0.00149144*pt*pt+0.00122299*pt*pt*pt-0.000583605*pt*pt*pt*pt;
                                }
                        }
                        else if(B==2){
                                if(C==0){
                                        return 0.0139405-0.0778337*pt+0.122313*pt*pt-0.0441081*pt*pt*pt+0.00502149*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                        return 0.0204243-0.0952528*pt+0.118049*pt*pt-0.0375978*pt*pt*pt+0.00388916*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return 0.0368641-0.132059*pt+0.140577*pt*pt-0.0465704*pt*pt*pt+0.00527664*pt*pt*pt*pt;
                                }
                        }
                }

                else if(A==1){
                        if(B==0){
                                if(C==0){
                                        return -0.0164604+0.119236*pt-0.0140501*pt*pt-0.00715798*pt*pt*pt+0.00137041*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                        return -0.00615225+0.0205719*pt+0.0373663*pt*pt-0.0176272*pt*pt*pt+0.00201875*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return 0.00343526-0.0156758*pt+0.058631*pt*pt-0.0234438*pt*pt*pt+0.00258536*pt*pt*pt*pt;
                                }
                        }
                        else if(B==1){
                                if(C==0){
                                        return -0.0335949+0.108473*pt+0.00189956*pt*pt-0.0119077*pt*pt*pt+0.00177362*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                        return -0.0141866+0.0239085*pt+0.0233996*pt*pt-0.0081926*pt*pt*pt+0.000430152*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return 0.0178734-0.0725294*pt+0.105302*pt*pt-0.0391775*pt*pt*pt+0.00466704*pt*pt*pt*pt;
                                }
                        }
                        else if(B==2){
                                if(C==0){
                                        return 0.0147481-0.0885341*pt+0.16892*pt*pt-0.0604128*pt*pt*pt+0.00661366*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                        return 0.020801-0.0910493*pt+0.118184*pt*pt-0.0365487*pt*pt*pt+0.0035618*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return 0.0300511-0.108966*pt+0.122315*pt*pt-0.0365423*pt*pt*pt+0.00350489*pt*pt*pt*pt;
                                }
                        }
                }

                else if(A==2){
                        if(B==0){
                                if(C==0){
                                        return -0.0220529+0.172125*pt-0.0353618*pt*pt-0.003559*pt*pt*pt+0.00113968*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                        return -0.0066372+0.0262161*pt+0.0372216*pt*pt-0.0187145*pt*pt*pt+0.00228567*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return 0.000152642+0.00135534*pt+0.0523496*pt*pt-0.0225954*pt*pt*pt+0.0025451*pt*pt*pt*pt;
                                }
                        }
                        else if(B==1){
                                if(C==0){
                                        return -0.0424241+0.152629*pt-0.00506494*pt*pt-0.0151633*pt*pt*pt+0.00254353*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                        return -0.0149325+0.0253627*pt+0.0329371*pt*pt-0.0153877*pt*pt*pt+0.00170996*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return -0.0171898+0.0261749*pt+0.032913*pt*pt-0.0180592*pt*pt*pt+0.00240376*pt*pt*pt*pt;
                                }
                        }
                        else if(B==2){
                                if(C==0){
                                        return 0.0128407-0.0812974*pt+0.196424*pt*pt-0.0729275*pt*pt*pt+0.0081403*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                        return 0.0216277-0.0905268*pt+0.125852*pt*pt-0.0410326*pt*pt*pt+0.00433817*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return 0.0296393-0.113592*pt+0.137947*pt*pt-0.0424535*pt*pt*pt+0.00422479*pt*pt*pt*pt;
                                }
                        }
                }

                else if(A==3){
                        if(B==0){
                                if(C==0){
                                        return -0.0273469+0.215291*pt-0.0580156*pt*pt+0.0015503*pt*pt*pt+0.00068957*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                        return -0.00634738+0.0244379*pt+0.0472794*pt*pt-0.0265474*pt*pt*pt+0.00383202*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return 0.00529299-0.0155944*pt+0.0851034*pt*pt-0.0399046*pt*pt*pt+0.00537977*pt*pt*pt*pt;
                                }
                        }
                        else if(B==1){
                                if(C==0){
                                        return -0.0457415+0.184931*pt-0.0184578*pt*pt-0.0130774*pt*pt*pt+0.00241422*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                        return 0.00059818-0.0174573*pt+0.0752039*pt*pt-0.0318181*pt*pt*pt+0.00386052*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return 0.00319935-0.0357498*pt+0.0956003*pt*pt-0.0389201*pt*pt*pt+0.0046787*pt*pt*pt*pt;
                                }
                        }
                        else if(B==2){
                                if(C==0){
                                        return 0.00914554-0.0597874*pt+0.203465*pt*pt-0.0797661*pt*pt*pt+0.00929514*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                        return 0.0304227-0.111558*pt+0.150866*pt*pt-0.0511995*pt*pt*pt+0.00556649*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return -0.0025491-0.0227755*pt+0.0628781*pt*pt-0.0165041*pt*pt*pt+0.00111185*pt*pt*pt*pt;
                                }
                        }
                }

                else if(A==4){
                        if(B==0){
                                if(C==0){
                                        return -0.0300557+0.23959*pt-0.0712208*pt*pt+0.004233*pt*pt*pt+0.000504197*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                        return -0.0047109+0.0195728*pt+0.0522525*pt*pt-0.0282469*pt*pt*pt+0.00377098*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return -0.0132215+0.0468*pt+0.0341852*pt*pt-0.0206421*pt*pt*pt+0.00294137*pt*pt*pt*pt;
                                }
                        }
                        else if(B==1){
                        if(C==0){
                                        return -0.0425067+0.191418*pt-0.0147714*pt*pt-0.0177701*pt*pt*pt+0.00341417*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                        return -0.000136675-0.0175618*pt+0.0863983*pt*pt-0.0430817*pt*pt*pt+0.00620464*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return -0.0570229+0.17675*pt-0.123802*pt*pt+0.0478088*pt*pt*pt-0.00638515*pt*pt*pt*pt;
                                }
                        }
                        else if(B==2){
                        if(C==0){
                                        return 0.0054852-0.0327023*pt+0.19693*pt*pt-0.0815048*pt*pt*pt+0.0098101*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                        return 0.0109575-0.0600514*pt+0.115052*pt*pt-0.0418587*pt*pt*pt+0.00470501*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return 0.04566-0.148386*pt+0.193706*pt*pt-0.0675996*pt*pt*pt+0.00792379*pt*pt*pt*pt;
                                }
                        }
                }

                else if(A==5){
                        if(B==0){
                                if(C==0){
                                        return -0.0327924+0.250176*pt-0.0765101*pt*pt+0.00390845*pt*pt*pt+0.000819225*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                        return -0.003819+0.0112537*pt+0.0588299*pt*pt-0.0333377*pt*pt*pt+0.0046983*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return -0.00266408+0.00640717*pt+0.023585*pt*pt-0.0121294*pt*pt*pt+0.00172664*pt*pt*pt*pt;
                                }
                        }
                        else if(B==1){
                                if(C==0){
                                        return -0.0131631+0.0325158*pt-0.00707803*pt*pt+0.000728541*pt*pt*pt-8.91768e-05*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                        return -0.0571195+0.249406*pt-0.074045*pt*pt+0.00463722*pt*pt*pt+0.000540562*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return -0.0143528+0.0402737*pt+0.0106742*pt*pt-0.00873702*pt*pt*pt+0.000978765*pt*pt*pt*pt;
                                }
                        }
                        else if(B==2){
                                if(C==0){
                                return -0.00897794+0.0202506*pt+0.159824*pt*pt-0.0719297*pt*pt*pt+0.00894275*pt*pt*pt*pt;
                                }
                                else if(C==1){
                                return -0.00854808+0.0237419*pt+0.00737678*pt*pt+0.00711372*pt*pt*pt-0.00275382*pt*pt*pt*pt;
                                }
                                else if(C==2){
                                        return 0.0;//0+1.764e-315*pt+2.122e-314*pt*pt+4.68595e-310*pt*pt*pt+1.94786e+118*pt*pt*pt*pt;
                                }
                        }
                }
        }
        //return 0.0;
}
float dNdPhi(float Phi, float V1, float V2, float V3, float V4, float PSI1, float PSI2, float PSI3, float PSI4){
        return (1.0 + 2.0*(V1*TMath::Cos(Phi-PSI1)+V2*TMath::Cos(2.0*(Phi-PSI2))+V3*TMath::Cos(3.0*(Phi-PSI3))+V4*TMath::Cos(4.0*(Phi-PSI4))))/(2.0*TMath::Pi());
}
int PartInd(int K_F){
     if(K_F==211||K_F==-211){return 0;}
     else if(K_F==321||K_F==-321){return 1;}
     else if(K_F==2212||K_F==-2212){return 2;}
     else {return 0;} 
}
string getPartString(int partInd){
        string particleID;
    	if(partInd ==0){particleID="pi-";}
	else if(partInd ==1){particleID="pi+";}
        else if(partInd ==2){particleID="ka-";}
        else if(partInd ==3){particleID="ka+";}
        else if(partInd ==4){particleID="pba";}
        else if(partInd ==5){particleID="pro";}
        else if(partInd ==6){particleID="ALL";}
		else{particleID = "ERROR_PARTINPUT";}
        return particleID;	
}
int getKFpartType(const int partInd ){
	if(partInd==0) return -211;
	else if(partInd == 1) return 211;
	else if(partInd == 2) return -321;
	else if(partInd == 3) return 321;
	else if(partInd == 4) return -2212;
	else if(partInd == 5) return 2212;
	else return 865;
}
