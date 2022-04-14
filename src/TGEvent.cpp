#include "TGEvent.h"

TENNGEN_BEGIN_NAMESPACE

TGEvent& TGEvent::operator+=( const TGEvent& addEvent) {
    if((addEvent.centrality == centrality)&& (addEvent.eta == eta)){
        for (int i = 0; i < addEvent.size(); ++i) {
            add( addEvent[i] );
           
        }
        multiplicity = particles.size();
    }
  return *this;
}

TGEvent& TGEvent::operator+=( const TGPartList& addList) {
 for (int i = 0; i < addList.size(); ++i) {
    add( addList[i] );
    
  }
  multiplicity = particles.size();
  return *this;
}

TGEvent& TGEvent::operator+=(const TGParticle& addParticle) {
    add( addParticle );
    multiplicity = particles.size();
    return *this;
}

TGEvent operator+(const TGEvent &a, const TGEvent &b){
    TGEvent sumEvent(a);
    if(a.centrality ==b.centrality && a.eta == b.eta){
        for(int i = 0; i<b.size(); i++){
            sumEvent.add(b[i]);
        }
        sumEvent.multiplicity = sumEvent.size();
    }
    return sumEvent;
}


TGEvent& TGEvent::operator=( const TGEvent& event) {
  // Do not copy if same.
  if (this != &event) {
    // Reset all current info in the list.
    clear();
    // Copy all the particles one by one.
    for (int i = 0; i < event.size(); ++i) add( event[i] );
    multiplicity=event.multiplicity;
    centrality=event.centrality;
    eta=event.eta;
  }
  return *this;

}

void TGEvent::setEventUserIndex(int in){
    for(int i=0;i<this->size();i++){
        particles.setListUserIndex(in);
    }
}

ostream& operator<<(ostream& os, const TGEvent& v) {
  os << fixed << setprecision(3) << v.multiplicity << " "
     << v.eta << " " <<  v.centrality <<  " ";
  return os;
}


float TGEvent::recoPsiN(const int n){
    float psiNum,psiDem;

    for(int i=0; i< particles.size(); i++ ){
      
        psiNum+=particles[i].Pt()*sin( 1.0*n*(M_PI + atan2( -particles[i].Py(),-particles[i].Py())));
        psiDem+=particles[i].Pt()*cos( 1.0*n*(M_PI + atan2( -particles[i].Py(),-particles[i].Py())));
    }
    return (M_PI+atan2(-psiNum,-psiDem))/(1.0*n);   
}

float TGEvent::avgPt(){
    float sum=0;
    TGParticle Tmp;
    for(int i=0; i< particles.size(); i++ ){
        Tmp = particles[i];
        sum+=Tmp.Pt();
    }
    return sum/multiplicity;
}

TGEventList& TGEventList::operator=( const TGEventList& eventlist) {
  // Do not copy if same.
  if (this != &eventlist) {
    // Reset all current info in the list.
    clear();
    // Copy all the particles one by one.
    for (int i = 0; i < eventlist.size(); ++i) add( eventlist[i] );
  }
  return *this;

}

void TGEventList::setEventListUserIndex(int in) {
    for(int i=0; i< this->size();i++) events[i].setEventUserIndex(in);
}



TGEventList& TGEventList::operator+=(const TGEventList& addList){
        for (int i = 0; i < addList.size(); ++i) {
            add( addList[i] );
        }
  return *this;
}
TGEventList& TGEventList::operator+=(const TGEvent& addEvent){
    add(addEvent);
    return *this;
}


TGEventList operator+(const TGEventList &a, const TGEventList &b){
    TGEventList sumList;
    if(a.size()>=b.size()) for(int i=0;i<b.size();i++) sumList+=(a[i]+b[i]);
    if(a.size()<b.size()) for(int i = 0; i< a.size();i++) sumList+=(a[i]+b[i]);
    return sumList;
}


void TGEventList::writeTTree(std::string& outdir){

    int nevents, eventcent;
    float eventeta;
    int eventMultiplicity;
    float TrackPt[2000];
    float TrackEta[2000];
    float TrackPhi[2000];
    int TrackKF[2000];

    TFile *tt = new TFile(outdir.c_str(),"RECREATE");
    TTree* outTree = new TTree("particleTree", "");
    TTree* eventTree = new TTree("EventInfo", "");

    eventTree->Branch("numEvents",&nevents,"numEvents/I");
    eventTree->Branch("centBin",&eventcent,"cent/I");
    eventTree->Branch("etaRange",&eventeta,"etaRange/F");

    outTree->Branch("multiplicity", &eventMultiplicity, "multiplicity/I");
    outTree->Branch("particlePt", TrackPt, "particlePt[multiplicity]/F");
    outTree->Branch("particleEta", TrackEta, "particleEta[multiplicity]/F");
    outTree->Branch("particlePhi", TrackPhi, "particlePhi[multiplicity]/F");
    outTree->Branch("particleKF", TrackKF, "particleE[multiplicity]/I");

    nevents = size();
    eventcent = Centrality(0);
    eventeta = Eta(0);

    eventTree->Fill();

    for(int i =0; i< nevents; i++){
    
        eventMultiplicity = (events[i]).size();
        for(int j=0;j<(events[i]).size();j++){
                TrackPt[j]= events[i][j].Pt();
                TrackEta[j]= events[i][j].Eta();
                TrackPhi[j]= events[i][j].Phi();
                TrackKF[j]=  events[i][j].Kf();
        }

        outTree->Fill();
        eventMultiplicity =0;
        for(int j =0;j<(events[i]).size(); j++){
                TrackPt[j]= 0;
                TrackEta[j]= 0;
                TrackPhi[j]= 0;
                TrackKF[j]=  0;
        }

    }
    tt->Write();
    tt->Close();

}

void TGEventList::writeHistos(std::string& outdir){

    //int setCent = Centrality(0);
    TFile *qa = new TFile(outdir.c_str(),"RECREATE");
    char expression1[128];
    char expression2[128];
    char expression3[128];
    char expression4[128];
    char expression5[128];
    char expression6[128];
    char expression7[128];
    char expression8[128];
    char expression9[128];
    char expression10[128];
    char expression11[128];
    char expression12[128];
    char expression13[128];
    char expression14[128];

    sprintf(expression1 , "p_{T} Distribution for #pi^{+} cent %d", Centrality(0));
    sprintf(expression2 , "p_{T} Distribution for #pi^{-} cent %d", Centrality(0));
    sprintf(expression3 , "p_{T} Distribution for K^{+} cent %d" , Centrality(0));
    sprintf(expression4 , "p_{T} Distribution for K^{-} cent %d" , Centrality(0));
    sprintf(expression5 , "p_{T} Distribution for p^{+} cent %d", Centrality(0));
    sprintf(expression6 , "p_{T} Distribution for p^{-} cent %d", Centrality(0));
    sprintf(expression7 , "p_{T} Distribution for all particles cent %d", Centrality(0));
    sprintf(expression8 , "#eta Distribution for all particles");
    sprintf(expression9 , "#phi Distribution for all particles");
    sprintf(expression10 , "#phi vs #eta distribution of all particles weighted by p_{T}");

    sprintf(expression11 , "#Psi_{Reco,1} Distribution for all particles");
    sprintf(expression12 , "#Psi_{Reco,2} Distribution for all particles");
    sprintf(expression13 , "#Psi_{Reco,3} Distribution for all particles");
    sprintf(expression14 , "#Psi_{Reco,4} Distribution for all particles");

    
    TH1D *histpT_piPlus = new TH1D("histpT_piPlus", expression1,200,0,10); 
    histpT_piPlus->SetXTitle("p_{T} (GeV/c)");
    histpT_piPlus->SetYTitle("dN/dp_{T}");

    TH1D *histpT_piMinus = new TH1D("histpT_piMinus", expression2,200,0,10); 
    histpT_piMinus->SetXTitle("p_{T} (GeV/c)");
    histpT_piMinus->SetYTitle("dN/dp_{T} ");

    TH1D *histpT_kPlus = new TH1D("histpT_kPlus", expression3,200,0,10);
    histpT_kPlus->SetXTitle("p_{T} (GeV/c)");
    histpT_kPlus->SetYTitle("dN/dp_{T} ");

    TH1D *histpT_kMinus = new TH1D("histpT_kMinus", expression4,200,0,10);
    histpT_kMinus->SetXTitle("p_{T} (GeV/c)");
    histpT_kMinus->SetYTitle("dN/dp_{T} ");

    TH1D *histpT_p = new TH1D("histpT_p", expression5,200,0,10); 
    histpT_p->SetYTitle("dN/dp_{T} ");
    histpT_p->SetYTitle("dN/dp_{T} ");

    TH1D *histpT_pbar = new TH1D("histpT_pbar", expression6,200,0,10); 
    histpT_pbar->SetXTitle("p_{T} (GeV/c)");
    histpT_pbar->SetYTitle("dN/dp_{T} ");

    TH1D *histpT_all = new TH1D("histpT_all", expression7,200,0,10); 
    histpT_all->SetXTitle("p_{T} (GeV/c)");
    histpT_all->SetYTitle("dN/dp_{T} ");

    TH1D *histeta_all = new TH1D("histeta_all", expression8,200,-1.5,1.5); 
    histeta_all->SetXTitle("eta");
    histeta_all->SetYTitle("dN/d#eta");

    TH1D *histphi_all = new TH1D("histphi_all", expression9,200,0,2*M_PI); 
    histphi_all->SetXTitle("phi");
    histphi_all->SetYTitle("dN/d#phi");

    TH2D *hist_pT_phi_eta_all = new TH2D("hist_pT_phi_eta_all",expression10,200,0,2*M_PI,200,-1.1,1.1); 
    hist_pT_phi_eta_all->SetXTitle("#phi (radians)");
    hist_pT_phi_eta_all->SetYTitle("#eta (pseudo-rapidity)");
    hist_pT_phi_eta_all->SetZTitle("dp_{T}^{ch.}/d#phid#eta");

    TH1D *histpsi_1 = new TH1D("histpsi_1", expression11,100,-0.5*M_PI,2.5*M_PI); 
    histpsi_1 ->SetXTitle("#Psi_{EP,1} (radians)");
    histpsi_1 ->SetYTitle("dN/d#Psi_{EP,1}");

    TH1D *histpsi_2 = new TH1D("histpsi_2", expression12,100,-0.5*M_PI,2.5*M_PI); 
    histpsi_2 ->SetXTitle("#Psi_{EP,2} (radians)");
    histpsi_2 ->SetYTitle("dN/d#Psi_{EP,2}");

    TH1D *histpsi_3 = new TH1D("histpsi_3", expression13,100,-0.5*M_PI,2.5*M_PI); 
    histpsi_3 ->SetXTitle("#Psi_{EP,3} (radians)");
    histpsi_3 ->SetYTitle("dN/d#Psi_{EP,3}");

    TH1D *histpsi_4 = new TH1D("histpsi_4", expression14,100,-0.5*M_PI,2.5*M_PI); 
    histpsi_4 ->SetXTitle("#Psi_{EP,4} (radians)");
    histpsi_4 ->SetYTitle("dN/d#Psi_{EP,4}");

    
 
    for(int i =0; i< size(); i++){
    

        histpsi_1->Fill(events[i].recoPsiN(1));
        histpsi_2->Fill(events[i].recoPsiN(2));
        histpsi_3->Fill(events[i].recoPsiN(3));
        histpsi_4->Fill(events[i].recoPsiN(4));

        for(int j =0; j< (events[i]).size(); j++ ){

            if((events[i])[j].Kf() == 211) histpT_piPlus->Fill((events[i])[j].Pt());
            if((events[i])[j].Kf() == -211) histpT_piMinus->Fill((events[i])[j].Pt());
            if((events[i])[j].Kf() == 321) histpT_kPlus->Fill((events[i])[j].Pt());
            if((events[i])[j].Kf() == -321) histpT_kMinus->Fill((events[i])[j].Pt());
            if((events[i])[j].Kf() == 2212) histpT_p->Fill((events[i])[j].Pt());
            if((events[i])[j].Kf() == -2212) histpT_pbar->Fill((events[i])[j].Pt());


            histphi_all->Fill((events[i])[j].Phi()); 
            histeta_all->Fill((events[i])[j].Eta()); 
            hist_pT_phi_eta_all->Fill((events[i])[j].Phi(),(events[i])[j].Eta(),(events[i])[j].Pt());
            histpT_all->Fill( (events[i])[j].Pt() );

        }

    }


    
    histpT_piPlus->Scale(1.0/histpT_piPlus->Integral());
    histpT_piMinus->Scale(1.0/histpT_piMinus->Integral());
    histpT_kPlus->Scale(1.0/histpT_kPlus->Integral());
    histpT_kMinus->Scale(1.0/histpT_kMinus->Integral());
    histpT_p->Scale(1.0/histpT_p->Integral());
    histpT_pbar->Scale(1.0/histpT_pbar->Integral());

    histphi_all->Scale(1.0/histphi_all->Integral());
    histeta_all->Scale(1.0/histeta_all->Integral());
    hist_pT_phi_eta_all->Scale(1.0/hist_pT_phi_eta_all->Integral());
    histphi_all->Scale(1.0/histphi_all->Integral());

    histpT_piPlus->Write();
    histpT_piMinus->Write();
    histpT_kPlus->Write();
    histpT_kMinus->Write();
    histpT_p->Write();
    histpT_pbar->Write();
    histphi_all->Write();
    histeta_all->Write();
    hist_pT_phi_eta_all->Write();
    histphi_all->Write();

    histpsi_1->Write();
    histpsi_2->Write();
    histpsi_3->Write();
    histpsi_4->Write();

    qa->Close();

}


TENNGEN_END_NAMESPACE
