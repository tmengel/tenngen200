#include <iostream>
#include <string>
#include <math.h>

#include "Pythia8/Pythia.h"

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"



using namespace Pythia8;


int PythiaGen(const int nEvent, const double pthardmin, const double pthardmax,const int ptbinnum,string pythiaRootDir) {


  string TfoutStr=pythiaRootDir+"/pythiaPtHardBin"+std::to_string(ptbinnum)+".root";
 // string TfoutStr="./"+pythiaRootDir+"/pythiaPtHardBinMerged.root";//+std::to_string(ptbinnum)+".root";
  TFile *Tfout = new TFile(TfoutStr.c_str(),"RECREATE"); 
  TTree* outTree = new TTree("particleTree", "");
  TTree* eventTree = new TTree("EventInfo", "");
  
  float scalingconst;
  int nevents, eventMultiplicity;
  int MAXTRACKS=2000;

  float etarange = 1.1;

  float TrackPt[MAXTRACKS];
  float TrackEta[MAXTRACKS];
  float TrackPhi[MAXTRACKS];
  float TrackE[MAXTRACKS];

  eventTree->Branch("numEvents",&nevents,"numEvents/I");
  eventTree->Branch("scalingconst",&scalingconst,"scalingconst/F");

  outTree->Branch("multiplicity", &eventMultiplicity, "multiplicity/I");
  outTree->Branch("particlePt", TrackPt, "particlePt[multiplicity]/F");
  outTree->Branch("particleEta", TrackEta, "particleEta[multiplicity]/F");
  outTree->Branch("particlePhi", TrackPhi, "particlePhi[multiplicity]/F");
  outTree->Branch("particleE", TrackE, "particleE[multiplicity]/F");

  double nRange= 1000;
  double pTrange = 100;

  // One does not need complete events to study pThard spectrum only.
  bool completeEvents = true;

  // Optionally minimize output (almost) to final results.
  bool smallOutput = true;

  string NumEventString = "Main::numberOfEvents =" +std::to_string(nEvent);

  Pythia pythia;

  // Shorthand for some public members of pythia (also static ones).
  Settings& settings = pythia.settings;
  const Info& info = pythia.info;
  Event& event = pythia.event;

  // Optionally limit output to minimal one.
  if (smallOutput) {
    pythia.readString("Init:showProcesses = off");
    pythia.readString("Init:showMultipartonInteractions = off");
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Next:numberCount = 1000000000");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
  }

    pythia.readString("Beams:idA =  2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 200.");
    pythia.readString("Tune:pp = 14 ");
    pythia.readString(NumEventString);
    pythia.readString("HardQCD:all = on");
    settings.parm("PhaseSpace:pTHatMin", pthardmin);
    settings.parm("PhaseSpace:pTHatMax", pthardmax);
    if (ptbinnum < 10) {
      if (!completeEvents) {
      pythia.readString("PartonLevel:all = on");
      pythia.readString("PartonLevel:ISR = off");
      pythia.readString("PartonLevel:FSR = off");
      pythia.readString("HadronLevel:all = off");
      }
    } 
    else {
      if (!completeEvents) pythia.readString("PartonLevel:all = off");
    }
   
    


    pythia.init();
    int nAbort = 5;
    int iAbort = 0;

    double weight;
    int j;
      
    // Begin event loop.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      weight = info.weight();
      // Generate events. Skip if failure.
      if (!pythia.next()) { /*Check to make sure pythia didn't crash or output junk*/
            if (++iAbort < nAbort) continue;
            cout << " Event generation aborted prematurely, owing to error!\n";
            break;
      }
      j=0;
      for (int i = 0; i < event.size(); ++i){
             if (event[i].isFinal() && event[i].isCharged()) {
                if(abs(event[i].eta())<=etarange){

                    TrackPt[j]= event[i].pT();
                    TrackEta[j]= event[i].eta();
                    TrackPhi[j]= event[i].phi()+M_PI;
                    TrackE[j]=  event[i].e();
                    j++;
                }
            }
       }
       eventMultiplicity =j;
       outTree->Fill();
       eventMultiplicity =0;
       for(int i=0;i<j;i++){
         TrackPt[i]= 0;
         TrackEta[i]= 0;
         TrackPhi[i]= 0;
         TrackE[i]= 0;
       }


    // End of event loop. Statistics.
    }
    if (!smallOutput) pythia.stat();

    // Normalize to cross section for each case, and add to sum.
    double sigmaNorm = (info.sigmaGen() / info.weightSum()) * (nRange / pTrange);

    scalingconst = sigmaNorm;
    nevents = nEvent;
    eventTree->Fill();

    Tfout->Write();
    Tfout->Close();

  
  return 0;
}
