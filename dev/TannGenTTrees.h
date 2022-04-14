#include "TannGen.h"
#include "TTree.h"

using namespace std;

int TannGenTTrees(const int nEvent, TRandom3* fRandom, string OutDir, const int nTrees, const int setCent = 0, const float etaRange = 1.1, const int doV1=1, const int doV2=1,const int doV3=1,const int doV4=1){

    int MAXTRACKS = 2000;
    int nevents, eventcent, eventMultiplicity;
    float eventeta;
    float TrackPt[MAXTRACKS];
    float TrackEta[MAXTRACKS];
    float TrackPhi[MAXTRACKS];
    int TrackKF[MAXTRACKS];

    string centString = "Cent"+std::to_string(setCent);

    string tFileOUT="./"+OutDir+"/TG_"+std::to_string(nEvent)+centString+"TTree"+std::to_string(nTrees)+".root";

    TFile* outFile = new TFile(tFileOUT.c_str(),"RECREATE");
    TTree* outTree = new TTree("particleTree", "");

    TTree* eventTree = new TTree("EventInfo", "");

    eventTree->Branch("numEvents",&nevents,"numEvents/I");
    eventTree->Branch("centBin",&eventcent,"cent/I");
    eventTree->Branch("etaRange",&eventeta,"etaRange/F");

    outTree->Branch("multiplicity", &eventMultiplicity, "multiplicity/I");
    outTree->Branch("particlePt", TrackPt, "particlePt[multiplicity]/F");
    outTree->Branch("particleEta", TrackEta, "particleEta[multiplicity]/F");
    outTree->Branch("particlePhi", TrackPhi, "particlePhi[multiplicity]/F");
    outTree->Branch("particleKF", TrackKF, "particleKF[multiplicity]/I");

     std::vector<TGEvent> myEvents =  TannGen(nEvent,fRandom,setCent,etaRange,doV1,doV2,doV3,doV4);

     nevents = myEvents.size();
     eventcent = setCent;
     eventeta = etaRange;

     eventTree->Fill();

     std::vector<TGParticle> EventParts;
     
     for(int i =0;i<nEvent;i++){

        EventParts=myEvents[i].parts();
        eventMultiplicity = EventParts.size();

        for(int j =0;j<EventParts.size(); j++){
                TrackPt[j]= EventParts[j].Pt();
                TrackEta[j]= EventParts[j].Eta();
                TrackPhi[j]= EventParts[j].Phi();
                TrackKF[j]=  EventParts[j].Kf();
        }

        outTree->Fill();
        eventMultiplicity = 0;

        for(int j =0;j<EventParts.size(); j++){
                TrackPt[j]= 0;
                TrackEta[j]= 0;
                TrackPhi[j]= 0;
                TrackKF[j]=  0;
        }

     }

     outFile->Write();
     delete eventTree;
     delete outTree;
     delete outFile;

     return 0;
}