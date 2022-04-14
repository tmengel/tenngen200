#ifndef TANNGENPYTHIA_H
#define TANNGENPYTHIA_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>



#include "Pythia8/Pythia.h"


using namespace Pythia8;
using namespace std;

const int MAXTracks = 2000;

void TannGenPythia(int nEvent,float pthardmin, float pthardmax, int batchnum, float etarange){



    string PtHatMinString = "PhaseSpace:pTHatMin = " + std::to_string(pthardmin);
    string PtHatMaxString = "PhaseSpace:pTHatMax = " + std::to_string(pthardmax);
    string NumEventString = "Main::numberOfEvents =" +std::to_string(nEvent);

    Pythia pythia;
    Event& event = pythia.event;

    pythia.readString(NumEventString);
    pythia.readString("Main:timesAllowErrors = 3");
    pythia.readString("Beams:idA =  2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 200");
    pythia.readString("HardQCD:all = on");
    pythia.readString("PromptPhoton:all = on");
    pythia.readString(PtHatMinString);
    pythia.readString(PtHatMaxString);
    pythia.readString("HadronLevel:Hadronize = on");
    pythia.readString("Tune:pp = 14 ");

    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("ParticleDecays:tau0Max = 10");

    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Next:numberCount = 100000");
    pythia.readString("Next:numberShowInfo = 0");
    
    pythia.init();

    int nAbort = 3;
    int iAbort = 0;

    for (int iEvent = 0; iEvent < nEvent; ++iEvent) { /*Creates event and fills the TTree */
       
        
        if (!pythia.next()) { /*Check to make sure pythia didn't crash or output junk*/
            if (++iAbort < nAbort) continue;
            cout << " Event generation aborted prematurely, owing to error!\n";
            break;
        }

        for (int i = 0; i < event.size(); ++i){
             if (event[i].isFinal() && event[i].isCharged()) {
                if(abs(event[i].eta())<=etarange){
                    event[i].px(),event[i].py(),event[i].pz(),event[i].e()
                    allparticles.push_back(PseudoJet(event[i].px(),event[i].py(),event[i].pz(),event[i].e()));
                    allparticles[numpart].set_user_index(truthind);
                    pythiaparticles.push_back(PseudoJet(event[i].px(),event[i].py(),event[i].pz(),event[i].e()));
                    pythiaparticles[truthpart].set_user_index(truthind);
                    numpart++;
                    truthpart++;
                }
        }
        

        for(int i =0;i<allparticles.size();i++){
            if(abs(allparticles[i].rap())>maxRap){maxRap = abs(allparticles[i].rap());}
        }
        
        ghost_maxrap = maxRap + (2.0*R);
        selectorRap = maxRap - R;

        fastjet::GhostedAreaSpec area_spec(ghost_maxrap);
        fastjet::AreaDefinition area_def(active_area_explicit_ghosts, area_spec);
        fastjet::AreaDefinition area_def_bkgd(active_area_explicit_ghosts, GhostedAreaSpec(ghost_maxrap));

        fastjet::Selector selector = SelectorAbsRapMax(selectorRap) * (!SelectorNHardest(2));
        fastjet::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
        fastjet::JetMedianBackgroundEstimator bkgd_estimator_pythia(selector, jet_def_bkgd, area_def_bkgd);

        fastjet::Subtractor subtractor(&bkgd_estimator);
        fastjet::Subtractor subtractor_pythia(&bkgd_estimator_pythia);


        #if FASTJET_VERSION_NUMBER >= 30100
        //subtractor.set_use_rho_m(true);
        //subtractor.set_safe_mass(true);
        #endif

        fastjet::ClusterSequenceArea cs_all(allparticles, jet_def, area_def);
        fastjet::ClusterSequenceArea cs_py(pythiaparticles, jet_def, area_def);
        fastjet::ClusterSequenceArea cs_tg(tgparticles, jet_def, area_def);

        vector<PseudoJet> all_jets = sorted_by_pt(cs_all.inclusive_jets(ptjetmin));
        vector<PseudoJet> pythia_jets = sorted_by_pt(cs_py.inclusive_jets(ptjetmin));
        vector<PseudoJet> fake_jets = sorted_by_pt(cs_tg.inclusive_jets(ptjetmin));

        bkgd_estimator.set_particles(allparticles);
        bkgd_estimator_pythia.set_particles(pythiaparticles);

        fastjet::BackgroundEstimate bkgd_estimate= bkgd_estimator.estimate();
        fastjet::BackgroundEstimate bkgd_estimate_pythia= bkgd_estimator_pythia.estimate();

    
        vector<PseudoJet> subtracted_all_jets = subtractor(all_jets);
        vector<PseudoJet> subtracted_pythia_jets = subtractor_pythia(pythia_jets);
        vector<PseudoJet> subtracted_fake_jets = subtractor(fake_jets);

        for (unsigned int i=0; i<all_jets.size(); i++){ /*fill TTrees with jet information*/
           
            if (all_jets[i].pt2() >= ptjetmin*ptjetmin && abs(all_jets[i].eta())<=(etarange-R)){

                vector<PseudoJet> constituents =  sorted_by_pt(all_jets[i].constituents());

                Event_BackgroundDensity = bkgd_estimate.rho();
                PtHardBin = pthardmin;
                JetPt = subtracted_all_jets[i].pt();
                JetPt_Raw = all_jets[i].pt();
                JetEta = all_jets[i].eta();
                JetPhi = all_jets[i].phi();
                JetArea = all_jets[i].area();
                NumTracks = constituents.size();

                for (unsigned int j = 0; j < constituents.size(); j++) {

                    TracksPt[j] = constituents[j].pt();
                    TracksEta[j] = constituents[j].eta();
                    TracksPhi[j] = constituents[j].phi();

                    totalPT+= constituents[j].pt();
                    JetAngularity+= TracksPt[j]*TMath::Sqrt((JetEta-TracksEta[j])*(JetEta-TracksEta[j])+(JetPhi-TracksPhi[j])*(JetPhi-TracksPhi[j]));

                    momentPT = (constituents[j].pt()-(JetPt_Raw/NumTracks));
                    
                    JetConstituent_1Moment+= TMath::Power(momentPT,1.0);
                    JetConstituent_2Moment+= TMath::Power(momentPT,2.0);
                    JetConstituent_3Moment+= TMath::Power(momentPT,3.0);
                    JetConstituent_4Moment+= TMath::Power(momentPT,4.0);

                    if(constituents[j].user_index()==1){
                        realPT+=constituents[j].pt();
                    }
                }


                JetConstituent_1Moment = JetConstituent_1Moment/NumTracks;
                JetConstituent_2Moment = JetConstituent_2Moment/NumTracks;
                JetConstituent_3Moment = JetConstituent_3Moment/NumTracks;
                JetConstituent_4Moment = JetConstituent_4Moment/NumTracks;

                JetAngularity= JetAngularity/JetPt_Raw;              
                TruePtFraction = realPT/totalPT;
                AverageConstituentPt = JetPt_Raw/NumTracks;


                outTree->Fill();
               
                for (unsigned int j = 0; j < NumTracks; j++) {
                    TracksPt[j] =0;
                    TracksEta[j] =0;
                    TracksPhi[j] = 0;                   
                }

                Event_BackgroundDensity =0;
                TruePtFraction = 0;    
                PtHardBin= 0;      
                JetPt = 0;
                JetPt_Raw = 0;
                JetEta =0;
                JetPhi = 0;
                JetArea = 0;
                JetAngularity = 0;
                NumTracks =0;
                AverageConstituentPt=0;
                JetConstituent_1Moment=0;
                JetConstituent_2Moment=0;
                JetConstituent_3Moment=0;
                JetConstituent_4Moment=0;

                realPT=0;
                totalPT =0;
                momentPT=0;
                
            }
       

        }


        for (unsigned int i=0; i<pythia_jets.size(); i++){ /*fill TTrees with jet information*/
           
            if (pythia_jets[i].pt2() >= ptjetmin*ptjetmin && abs(subtracted_pythia_jets[i].eta())<=(etarange-R)){

                vector<PseudoJet> constituents_pythia =  sorted_by_pt(pythia_jets[i].constituents());

                Event_BackgroundDensity_pythia = bkgd_estimate_pythia.rho();
                PtHardBin_pythia = pthardmin;
                JetPt_pythia = subtracted_pythia_jets[i].pt();
                JetPt_Raw_pythia = pythia_jets[i].pt();
                JetEta_pythia = pythia_jets[i].eta();
                JetPhi_pythia = pythia_jets[i].phi();
                JetArea_pythia = pythia_jets[i].area();
                NumTracks_pythia = constituents_pythia.size();

                for (unsigned int j = 0; j < constituents_pythia.size(); j++) {

                    TracksPt_pythia[j] = constituents_pythia[j].pt();
                    TracksEta_pythia[j] = constituents_pythia[j].eta();
                    TracksPhi_pythia[j] = constituents_pythia[j].phi();

                    totalPT+= constituents_pythia[j].pt();
                    JetAngularity_pythia+= TracksPt_pythia[j]*TMath::Sqrt((JetEta_pythia-TracksEta_pythia[j])*(JetEta_pythia-TracksEta_pythia[j])+(JetPhi_pythia-TracksPhi_pythia[j])*(JetPhi_pythia-TracksPhi_pythia[j]));

                    momentPT = (constituents_pythia[j].pt()-(JetPt_pythia/NumTracks_pythia));
                    
                    JetConstituent_1Moment_pythia+= TMath::Power(momentPT,1.0);
                    JetConstituent_2Moment_pythia+= TMath::Power(momentPT,2.0);
                    JetConstituent_3Moment_pythia+= TMath::Power(momentPT,3.0);
                    JetConstituent_4Moment_pythia+= TMath::Power(momentPT,4.0);

                    if(constituents_pythia[j].user_index()==1){
                        realPT+=constituents_pythia[j].pt();
                    }
                }

                JetConstituent_1Moment_pythia = JetConstituent_1Moment_pythia/NumTracks_pythia;
                JetConstituent_2Moment_pythia = JetConstituent_2Moment_pythia/NumTracks_pythia;
                JetConstituent_3Moment_pythia = JetConstituent_3Moment_pythia/NumTracks_pythia;
                JetConstituent_4Moment_pythia = JetConstituent_4Moment_pythia/NumTracks_pythia;

                JetAngularity_pythia= JetAngularity_pythia/JetPt_Raw_pythia;              
                TruePtFraction_pythia = realPT/totalPT;
                AverageConstituentPt_pythia = JetPt_Raw_pythia/NumTracks_pythia;


                pythiaTree->Fill();
               
                for (unsigned int j = 0; j < NumTracks_pythia; j++) {
                    TracksPt_pythia[j] =0;
                    TracksEta_pythia[j] =0;
                    TracksPhi_pythia[j] = 0;
                   
                }

                Event_BackgroundDensity_pythia =0;
                TruePtFraction_pythia = 0;    
                PtHardBin_pythia= 0;      
                JetPt_pythia = 0;
                JetPt_Raw_pythia = 0;
                JetEta_pythia =0;
                JetPhi_pythia = 0;
                JetArea_pythia = 0;
                JetAngularity_pythia = 0;
                NumTracks_pythia =0;
                AverageConstituentPt_pythia=0;
                JetConstituent_1Moment_pythia=0;
                JetConstituent_2Moment_pythia=0;
                JetConstituent_3Moment_pythia=0;
                JetConstituent_4Moment_pythia=0;

                realPT=0;
                totalPT =0;
                momentPT=0;
                
            }
       

        }


        for (unsigned int i=0; i<fake_jets.size(); i++){ /*fill TTrees with jet information*/
           
            if (fake_jets[i].pt2() >= ptjetmin*ptjetmin && abs(subtracted_fake_jets[i].eta())<=(etarange-R)){

                vector<PseudoJet> constituents_background =  sorted_by_pt(fake_jets[i].constituents());

                Event_BackgroundDensity_background = bkgd_estimate.rho();
                PtHardBin_background = pthardmin;
                JetPt_background = subtracted_fake_jets[i].pt();
                JetPt_Raw_background = fake_jets[i].pt();
                JetEta_background = fake_jets[i].eta();
                JetPhi_background = fake_jets[i].phi();
                JetArea_background = fake_jets[i].area();
                NumTracks_background = constituents_background.size();

                for (unsigned int j = 0; j < constituents_background.size(); j++) {

                    TracksPt_background[j] = constituents_background[j].pt();
                    TracksEta_background[j] = constituents_background[j].eta();
                    TracksPhi_background[j] = constituents_background[j].phi();

                    totalPT+= constituents_background[j].pt();
                    JetAngularity_background+= TracksPt_background[j]*TMath::Sqrt((JetEta_background-TracksEta_background[j])*(JetEta_background-TracksEta_background[j])+(JetPhi_background-TracksPhi_background[j])*(JetPhi_background-TracksPhi_background[j]));

                    momentPT = (constituents_background[j].pt()-(JetPt_background/NumTracks_background));
                    JetConstituent_1Moment_background+= TMath::Power(momentPT,1.0);
                    JetConstituent_2Moment_background+= TMath::Power(momentPT,2.0);
                    JetConstituent_3Moment_background+= TMath::Power(momentPT,3.0);
                    JetConstituent_4Moment_background+= TMath::Power(momentPT,4.0);

                    if(constituents_background[j].user_index()==1){
                        realPT+=constituents_background[j].pt();
                    }
                }

                JetConstituent_1Moment_background = JetConstituent_1Moment_background/NumTracks_background;
                JetConstituent_2Moment_background = JetConstituent_2Moment_background/NumTracks_background;
                JetConstituent_3Moment_background = JetConstituent_3Moment_background/NumTracks_background;
                JetConstituent_4Moment_background = JetConstituent_4Moment_background/NumTracks_background;

                JetAngularity_background= JetAngularity_background/JetPt_Raw_background;              
                TruePtFraction_background = realPT/totalPT;
                AverageConstituentPt_background = JetPt_Raw_background/NumTracks_background;


                backgroundTree->Fill();
               
                for (unsigned int j = 0; j < NumTracks_background; j++) {
                    TracksPt_background[j] =0;
                    TracksEta_background[j] =0;
                    TracksPhi_background[j] = 0;
                   
                }

                Event_BackgroundDensity_background =0;
                TruePtFraction_background = 0;    
                PtHardBin_background= 0;      
                JetPt_background = 0;
                JetPt_Raw_background = 0;
                JetEta_background =0;
                JetPhi_background = 0;
                JetArea_background = 0;
                JetAngularity_background = 0;
                NumTracks_background =0;
                AverageConstituentPt_background=0;
                JetConstituent_1Moment_background=0;
                JetConstituent_2Moment_background=0;
                JetConstituent_3Moment_background=0;
                JetConstituent_4Moment_background=0;

                realPT=0;
                totalPT =0;
                momentPT=0;
                
            }
       

        }


        if(iEvent%(nEvent-1)==0||iEvent%(nEvent/10)==0){
        printf("-----------------------------Event Summary-------------------------------------\n");
        printf("-------------------------------------------------------------------------------\n");
        printf("%10s %10s %15s  %15s\n","evnt #", "jets", "tracks",  "ptmin");
        printf("%10i %10i %15i  %15.8f\n",iEvent, all_jets.size(), numpart,  pthardmin);
        printf("-------------------------------------------------------------------------------\n");
        printf("-------------------------------------------------------------------------------\n");
        }

        
    } 
    gSystem->ProcessEvents();
    // gROOT->Reset();
    outFile->Write();
    delete outTree;
    delete pythiaTree;
    delete backgroundTree;
    



}



#endif