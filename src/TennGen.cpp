#include "TennGen.h"
#include <bits/stdc++.h>
#include <iostream>
#include <sys/stat.h>
#include "TGraphErrors.h"
#include <sys/types.h>

TENNGEN_BEGIN_NAMESPACE

void TennGen::init(){
  dir = settings.getOutputDir(); 

  if (mkdir(dir.c_str(), 0777) == -1) std::cerr << "Output " << strerror(errno) << endl;
        
  settings.setOutputDir(settings.getOutputDir()+"/");   

  if(settings.batchmode()) runBatch();
  
  if(settings.qamode()) runQA();

  if(settings.streammode()) runStream();

  std::cout << "======================================="<<std::endl;
  std::cout << "======================================="<<std::endl;
                        
}
void TennGen::runBatch(){
        if(settings.batchmode()){

                std::cout << "======================================="<<std::endl;
                std::cout << "TennGen Settings: BATCH MODE "  << std::endl ;
                std::cout << "======================================="<<std::endl;
                std::cout << "collison energy: " << settings.getCollEn() <<std::endl;
                std::cout << "nEvents: " << settings.getNevents() << std::endl;
                std::cout << "Cent: " << settings.getCentBin() << std::endl;
                std::cout << "Eta: " << settings.getEtaRange() << std::endl;

                if(settings.getVN(1) == 1) std::cout << "v1: On " <<  std::endl;
                if(settings.getVN(1) == 0) std::cout << "v1: Off " <<  std::endl;

                if(settings.getVN(2) == 1) std::cout << "v2: On " <<  std::endl;
                if(settings.getVN(2) == 0) std::cout << "v2: Off " <<  std::endl;

                if(settings.getVN(3) == 1) std::cout << "v3: On " <<  std::endl;
                if(settings.getVN(3) == 0) std::cout << "v3: Off " <<  std::endl;

                if(settings.getVN(4) == 1) std::cout << "v4: On " <<  std::endl;
                if(settings.getVN(4) == 0) std::cout << "v4: Off " <<  std::endl;

                if(settings.getCollEn()==5020){
                        if(settings.getVN(5) == 1) std::cout << "v5: On " <<  std::endl;
                        if(settings.getVN(5) == 0) std::cout << "v5: Off " <<  std::endl;
                }
                
                if(settings.getPsiN(1) == -1.0)  std::cout << "psi1: Random" << std::endl;
                if(settings.getPsiN(1) != -1.0)  std::cout << "psi1: " << settings.getPsiN(1) << std::endl;

                if(settings.getPsiN(2) == -1.0)  std::cout << "psi2: Random" << std::endl;
                if(settings.getPsiN(2) != -1.0)  std::cout << "psi2: " << settings.getPsiN(2) << std::endl;

                if(settings.getPsiN(3) == -1.0)  std::cout << "psi3: Random" << std::endl;
                if(settings.getPsiN(3) != -1.0)  std::cout << "psi3: " << settings.getPsiN(3) << std::endl;

                if(settings.getPsiN(4) == -1.0)  std::cout << "psi4: Random" << std::endl;
                if(settings.getPsiN(4) != -1.0)  std::cout << "psi4: " << settings.getPsiN(4) << std::endl;

                if(settings.getCollEn()==5020){
                        if(settings.getPsiN(5) == -1.0)  std::cout << "psi5: Random" << std::endl;
                        if(settings.getPsiN(5) != -1.0)  std::cout << "psi5: " << settings.getPsiN(5) << std::endl;
                }

                std::cout << "OutDir: " << settings.getOutputDir() << std::endl;

                if(settings.getHistos()) std::cout << "Histos: true "<< std::endl;
                if(!settings.getHistos()) std::cout << "Histos: false "<< std::endl;


                if(settings.getTTree()) std::cout << "TTree: true "<< std::endl;
                if(!settings.getTTree()) std::cout << "TTree: false "<< std::endl;

                std::cout<< "======================================="<<std::endl;
                std::cout<< "======================================="<<std::endl;

                if(settings.getNevents()>50000){

                        int lastbatch = settings.getNevents()%50000;
                        int numbatches = (settings.getNevents()-lastbatch)/50000;
                        int filelabel;
                        if(lastbatch !=0) filelabel = numbatches+1;
                        else filelabel = numbatches;
                        TGSettings TempParams = settings;
                        TempParams.setNevents(50000);
                        if(settings.getCollEn()==200){
                                tg.config(TempParams);
                                tg.seed(fRandom);
                        }

                        for(int k =0; k< (numbatches);k++){

                                TempParams.setOutputDir(settings.getOutputDir()+std::to_string(k)+"_of_"+std::to_string(filelabel)+"_");
                               // tFileHistos=TempParams.getOutputDir()+"TG_"+std::to_string(TempParams.getNevents())+"_eta"+std::to_string(settings.getEtaRange())+"_"+"Cent"+std::to_string(settings.getCentBin())+"_"+std::to_string(settings.getVN(1))+std::to_string(settings.getVN(2))+std::to_string(settings.getVN(3))+std::to_string(settings.getVN(4))+"_Histos.root";
                                //tFileTree=TempParams.getOutputDir()+"TG_"+std::to_string(TempParams.getNevents())+"_eta"+std::to_string(settings.getEtaRange())+"_"+"Cent"+std::to_string(settings.getCentBin())+"_"+std::to_string(settings.getVN(1))+std::to_string(settings.getVN(2))+std::to_string(settings.getVN(3))+std::to_string(settings.getVN(4))+"_TTree.root";

                                std::cout << "OutDir(s): " << TempParams.getOutputDir() << std::endl;
                                if(settings.getCollEn()==200) TG200 tgBatch(TempParams, fRandom);
                               
                                //if(settings.getCollEn()==5020) TG5020 tg(TempParams,fRandom);
                              
                                
                                std::cout << "Event " << (k+1)*50000 << " of " << settings.getNevents() << " completed" <<std::endl;
                                //events.clear();

                        }

                        if(lastbatch!=0){
                                TempParams.setNevents(lastbatch);
                                TempParams.setOutputDir(settings.getOutputDir()+std::to_string(numbatches+1)+"_of_"+std::to_string(filelabel)+"_");
                                tFileHistos=TempParams.getOutputDir()+"TG_"+std::to_string(TempParams.getNevents())+"_eta"+std::to_string(settings.getEtaRange())+"_"+"Cent"+std::to_string(settings.getCentBin())+"_"+std::to_string(settings.getVN(1))+std::to_string(settings.getVN(2))+std::to_string(settings.getVN(3))+std::to_string(settings.getVN(4))+"_Histos.root";
                                tFileTree=TempParams.getOutputDir()+"TG_"+std::to_string(TempParams.getNevents())+"_eta"+std::to_string(settings.getEtaRange())+"_"+"Cent"+std::to_string(settings.getCentBin())+"_"+std::to_string(settings.getVN(1))+std::to_string(settings.getVN(2))+std::to_string(settings.getVN(3))+std::to_string(settings.getVN(4))+"_TTree.root";
                                if(settings.getCollEn()==200) TG200 tgBatch(TempParams, fRandom);
                                
                                //if(settings.getCollEn()==5020) TG5020 tg(TempParams,fRandom);
                                //std::cout << "OutDir(s): " << TempParams.getOutputDir() << std::endl;
                
                                std::cout << "Event " << events.size() << " of " << settings.getNevents() << " completed" <<std::endl;
                        }
               

                }
                else{
                
                        if(settings.getCollEn()==200){
                                TG200 tgBatch(settings, fRandom);
                                events = tgBatch.events();
                        } 
                        // if(settings.getCollEn()==5020) TG5020 tg(TempParams,fRandom);
                        std::cout << "Event "<< events.size() << " completed." << std::endl;
                        std::cout << "======================================="<<std::endl;
                        tFileHistos=settings.getOutputDir()+"TG_"+std::to_string(settings.getNevents())+"_eta"+std::to_string(settings.getEtaRange())+"_"+"Cent"+std::to_string(settings.getCentBin())+"_"+std::to_string(settings.getVN(1))+std::to_string(settings.getVN(2))+std::to_string(settings.getVN(3))+std::to_string(settings.getVN(4))+"_Histos.root";
                        tFileTree=settings.getOutputDir()+"TG_"+std::to_string(settings.getNevents())+"_eta"+std::to_string(settings.getEtaRange())+"_"+"Cent"+std::to_string(settings.getCentBin())+"_"+std::to_string(settings.getVN(1))+std::to_string(settings.getVN(2))+std::to_string(settings.getVN(3))+std::to_string(settings.getVN(4))+"_TTree.root";
                        //events.clear();
  
                }

        }
        else {cout<<"ERROR TENNGEN NOT IN STREAM MODE"<<endl;}
}
void TennGen::runQA(){
        if(settings.qamode()){
                std::cout << "======================================="<<std::endl;
                std::cout << "Running TennGen QA "  << std::endl ;
                std::cout << "======================================="<<std::endl;
                std::cout << "collison energy: " << settings.getCollEn() <<std::endl;
                std::cout << "nEvents: 10000000"  << std::endl;
                std::cout << "Cent: ALL" << std::endl;
                std::cout << "Eta: " << settings.getEtaRange() << std::endl;

                TG200 tgQA(settings, fRandom);

        }
        else {cout<<"ERROR TENNGEN NOT IN QA MODE"<<endl;}
}
void TennGen::runStream(){

        if(settings.streammode()){
                std::cout << "======================================="<<std::endl;
                std::cout << "TennGen Settings: STREAMING "  << std::endl ;
                std::cout << "======================================="<<std::endl;
                std::cout << "collison energy: " << settings.getCollEn() <<std::endl;
                //std::cout << "nEvents: " << settings.getNevents() << std::endl;
                std::cout << "Cent: " << settings.getCentBin() << std::endl;
                std::cout << "Eta: " << settings.getEtaRange() << std::endl;

                if(settings.getVN(1) == 1) std::cout << "v1: On " <<  std::endl;
                if(settings.getVN(1) == 0) std::cout << "v1: Off " <<  std::endl;

                if(settings.getVN(2) == 1) std::cout << "v2: On " <<  std::endl;
                if(settings.getVN(2) == 0) std::cout << "v2: Off " <<  std::endl;

                if(settings.getVN(3) == 1) std::cout << "v3: On " <<  std::endl;
                if(settings.getVN(3) == 0) std::cout << "v3: Off " <<  std::endl;

                if(settings.getVN(4) == 1) std::cout << "v4: On " <<  std::endl;
                if(settings.getVN(4) == 0) std::cout << "v4: Off " <<  std::endl;

                if(settings.getCollEn()==5020){
                        if(settings.getVN(5) == 1) std::cout << "v5: On " <<  std::endl;
                        if(settings.getVN(5) == 0) std::cout << "v5: Off " <<  std::endl;
                }
                
                if(settings.getPsiN(1) == -1.0)  std::cout << "psi1: Random" << std::endl;
                if(settings.getPsiN(1) != -1.0)  std::cout << "psi1: " << settings.getPsiN(1) << std::endl;

                if(settings.getPsiN(2) == -1.0)  std::cout << "psi2: Random" << std::endl;
                if(settings.getPsiN(2) != -1.0)  std::cout << "psi2: " << settings.getPsiN(2) << std::endl;

                if(settings.getPsiN(3) == -1.0)  std::cout << "psi3: Random" << std::endl;
                if(settings.getPsiN(3) != -1.0)  std::cout << "psi3: " << settings.getPsiN(3) << std::endl;

                if(settings.getPsiN(4) == -1.0)  std::cout << "psi4: Random" << std::endl;
                if(settings.getPsiN(4) != -1.0)  std::cout << "psi4: " << settings.getPsiN(4) << std::endl;

                if(settings.getCollEn()==5020){
                        if(settings.getPsiN(5) == -1.0)  std::cout << "psi5: Random" << std::endl;
                        if(settings.getPsiN(5) != -1.0)  std::cout << "psi5: " << settings.getPsiN(5) << std::endl;
                }
                TGSettings TempParams = settings;
                TempParams.setNevents(50000);
                if(settings.getCollEn()==200){
                                tg.config(TempParams);
                                tg.seed(fRandom);
                                tg.genEvents();
                } 
                STREAMING = true;
                StreamInt =0;
                

        }
        else {cout<<"ERROR TENNGEN NOT IN STREAM MODE"<<endl;}

}
TGEvent& TennGen::next(){
        if(STREAMING){
                if(StreamInt == 49999){
                        tg.genEvents();
                        StreamInt = 0;
                }
                StreamInt++;
                return tg[StreamInt-1];
        }
        else{
                settings.setStream(true);
                runStream();
                if(StreamInt == 49999){
                        tg.genEvents();
                        StreamInt = 0;
                }
                StreamInt++;
                return tg[StreamInt-1];
        }



}
const int TG200::raw_high[]={608,430,311,145};
const int TG200::raw_low[]={431,312,146,56};
const int TG200::raw_bin[]= {177,118,165,89};

const int TG200::correction_high[]={778,499,344,153};
const int TG200::correction_low[]={500,345,154,58};


const float TG200::piplusRatios[] = {0.396738386,0.400412797,0.400958581,0.40230616};
const float TG200::piminusRatios[] = {0.396738386,0.400412797,0.400958581,0.40230616};
const float TG200::kaplusRatios[] = {0.063108127,0.061919505,0.060984334,0.059731231};
const float TG200::kaminusRatios[] = {0.06118953,0.059236326, 0.058838257,0.057484421};
const float TG200::proRatios[] = {0.043099904,0.041486068,0.042385006,0.042604604};
const float TG200::pbarRatios[] = {0.03295875,0.032404541,0.033371486,0.034295646};

float TG200::HarmonicFunction(const int harmonicN){
                       
                float pt;
                if(ptPart > 4.0){  pt = 4.0;  }
                else if(ptPart < 0.5){ return 0.0; }
                else{ pt = ptPart; }
           
                
                if(harmonicCent==0){
                        if( KF==211|| KF==-211){
                                if(harmonicN==0){
                                return -0.0143392+0.0759773*pt-0.0154024*pt*pt-0.000915513*pt*pt*pt+0.000264875*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                return -0.00558048+0.0136122*pt+0.0356842*pt*pt-0.0164884*pt*pt*pt+0.00195171*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return -0.00256332-0.00483424*pt+0.0390473*pt*pt-0.016438*pt*pt*pt+0.00205353*pt*pt*pt*pt;
                                }
                        }
                        else if(KF==321|| KF==-321){
                                if(harmonicN==0){
                                return -0.0287711+0.0777838*pt-0.0110908*pt*pt-0.00251265*pt*pt*pt+0.000435634*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                return -0.00843333+0.00413273*pt+0.03708*pt*pt-0.0140296*pt*pt*pt+0.00138853*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return -0.0147847+0.0238129*pt+0.00149144*pt*pt+0.00122299*pt*pt*pt-0.000583605*pt*pt*pt*pt;
                                }
                        }
                        else if(KF==2212|| KF==-2212){
                                if(harmonicN==0){
                                return 0.0139405-0.0778337*pt+0.122313*pt*pt-0.0441081*pt*pt*pt+0.00502149*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                return 0.0204243-0.0952528*pt+0.118049*pt*pt-0.0375978*pt*pt*pt+0.00388916*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return 0.0368641-0.132059*pt+0.140577*pt*pt-0.0465704*pt*pt*pt+0.00527664*pt*pt*pt*pt;
                                }
                        }
                }

                else if(harmonicCent==1){
                        if( KF==211|| KF==-211){
                                if(harmonicN==0){
                                return -0.0164604+0.119236*pt-0.0140501*pt*pt-0.00715798*pt*pt*pt+0.00137041*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                return -0.00615225+0.0205719*pt+0.0373663*pt*pt-0.0176272*pt*pt*pt+0.00201875*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return 0.00343526-0.0156758*pt+0.058631*pt*pt-0.0234438*pt*pt*pt+0.00258536*pt*pt*pt*pt;
                                }
                        }
                        else if(KF==321|| KF==-321){
                                if(harmonicN==0){
                                return -0.0335949+0.108473*pt+0.00189956*pt*pt-0.0119077*pt*pt*pt+0.00177362*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                return -0.0141866+0.0239085*pt+0.0233996*pt*pt-0.0081926*pt*pt*pt+0.000430152*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return 0.0178734-0.0725294*pt+0.105302*pt*pt-0.0391775*pt*pt*pt+0.00466704*pt*pt*pt*pt;
                                }
                        }
                        else if(KF==2212|| KF==-2212){
                                if(harmonicN==0){
                                return 0.0147481-0.0885341*pt+0.16892*pt*pt-0.0604128*pt*pt*pt+0.00661366*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                return 0.020801-0.0910493*pt+0.118184*pt*pt-0.0365487*pt*pt*pt+0.0035618*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return 0.0300511-0.108966*pt+0.122315*pt*pt-0.0365423*pt*pt*pt+0.00350489*pt*pt*pt*pt;
                                }
                        }
                }

                else if(harmonicCent==2){
                        if( KF==211|| KF==-211){
                                if(harmonicN==0){
                                return -0.0220529+0.172125*pt-0.0353618*pt*pt-0.003559*pt*pt*pt+0.00113968*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                return -0.0066372+0.0262161*pt+0.0372216*pt*pt-0.0187145*pt*pt*pt+0.00228567*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return 0.000152642+0.00135534*pt+0.0523496*pt*pt-0.0225954*pt*pt*pt+0.0025451*pt*pt*pt*pt;
                                }
                        }
                        else if(KF==321|| KF==-321){
                                if(harmonicN==0){
                                return -0.0424241+0.152629*pt-0.00506494*pt*pt-0.0151633*pt*pt*pt+0.00254353*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                return -0.0149325+0.0253627*pt+0.0329371*pt*pt-0.0153877*pt*pt*pt+0.00170996*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return -0.0171898+0.0261749*pt+0.032913*pt*pt-0.0180592*pt*pt*pt+0.00240376*pt*pt*pt*pt;
                                }
                        }
                        else if(KF==2212|| KF==-2212){
                                if(harmonicN==0){
                                return 0.0128407-0.0812974*pt+0.196424*pt*pt-0.0729275*pt*pt*pt+0.0081403*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                return 0.0216277-0.0905268*pt+0.125852*pt*pt-0.0410326*pt*pt*pt+0.00433817*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return 0.0296393-0.113592*pt+0.137947*pt*pt-0.0424535*pt*pt*pt+0.00422479*pt*pt*pt*pt;
                                }
                        }
                }

                else if(harmonicCent==3){
                        if( KF==211|| KF==-211){
                                if(harmonicN==0){
                                return -0.0273469+0.215291*pt-0.0580156*pt*pt+0.0015503*pt*pt*pt+0.00068957*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                return -0.00634738+0.0244379*pt+0.0472794*pt*pt-0.0265474*pt*pt*pt+0.00383202*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return 0.00529299-0.0155944*pt+0.0851034*pt*pt-0.0399046*pt*pt*pt+0.00537977*pt*pt*pt*pt;
                                }
                        }
                        else if(KF==321|| KF==-321){
                                if(harmonicN==0){
                                return -0.0457415+0.184931*pt-0.0184578*pt*pt-0.0130774*pt*pt*pt+0.00241422*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                return 0.00059818-0.0174573*pt+0.0752039*pt*pt-0.0318181*pt*pt*pt+0.00386052*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return 0.00319935-0.0357498*pt+0.0956003*pt*pt-0.0389201*pt*pt*pt+0.0046787*pt*pt*pt*pt;
                                }
                        }
                        else if(KF==2212|| KF==-2212){
                                if(harmonicN==0){
                                return 0.00914554-0.0597874*pt+0.203465*pt*pt-0.0797661*pt*pt*pt+0.00929514*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                return 0.0304227-0.111558*pt+0.150866*pt*pt-0.0511995*pt*pt*pt+0.00556649*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return -0.0025491-0.0227755*pt+0.0628781*pt*pt-0.0165041*pt*pt*pt+0.00111185*pt*pt*pt*pt;
                                }
                        }
                }

                else if(harmonicCent==4){
                        if( KF==211|| KF==-211){
                                if(harmonicN==0){
                                return -0.0300557+0.23959*pt-0.0712208*pt*pt+0.004233*pt*pt*pt+0.000504197*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                return -0.0047109+0.0195728*pt+0.0522525*pt*pt-0.0282469*pt*pt*pt+0.00377098*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return -0.0132215+0.0468*pt+0.0341852*pt*pt-0.0206421*pt*pt*pt+0.00294137*pt*pt*pt*pt;
                                }
                        }
                        else if(KF==321|| KF==-321){
                        if(harmonicN==0){
                                return -0.0425067+0.191418*pt-0.0147714*pt*pt-0.0177701*pt*pt*pt+0.00341417*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                return -0.000136675-0.0175618*pt+0.0863983*pt*pt-0.0430817*pt*pt*pt+0.00620464*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return -0.0570229+0.17675*pt-0.123802*pt*pt+0.0478088*pt*pt*pt-0.00638515*pt*pt*pt*pt;
                                }
                        }
                        else if(KF==2212|| KF==-2212){
                        if(harmonicN==0){
                                return 0.0054852-0.0327023*pt+0.19693*pt*pt-0.0815048*pt*pt*pt+0.0098101*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                return 0.0109575-0.0600514*pt+0.115052*pt*pt-0.0418587*pt*pt*pt+0.00470501*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return 0.04566-0.148386*pt+0.193706*pt*pt-0.0675996*pt*pt*pt+0.00792379*pt*pt*pt*pt;
                                }
                        }
                }

                else if(harmonicCent==5){
                        if( KF==211|| KF==-211){
                                if(harmonicN==0){
                                return -0.0327924+0.250176*pt-0.0765101*pt*pt+0.00390845*pt*pt*pt+0.000819225*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                return -0.003819+0.0112537*pt+0.0588299*pt*pt-0.0333377*pt*pt*pt+0.0046983*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return -0.00266408+0.00640717*pt+0.023585*pt*pt-0.0121294*pt*pt*pt+0.00172664*pt*pt*pt*pt;
                                }
                        }
                        else if(KF==321|| KF==-321){
                                if(harmonicN==0){
                                return -0.0131631+0.0325158*pt-0.00707803*pt*pt+0.000728541*pt*pt*pt-8.91768e-05*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                return -0.0571195+0.249406*pt-0.074045*pt*pt+0.00463722*pt*pt*pt+0.000540562*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return -0.0143528+0.0402737*pt+0.0106742*pt*pt-0.00873702*pt*pt*pt+0.000978765*pt*pt*pt*pt;
                                }
                        }
                        else if(KF==2212|| KF==-2212){
                                if(harmonicN==0){
                                                return -0.00897794+0.0202506*pt+0.159824*pt*pt-0.0719297*pt*pt*pt+0.00894275*pt*pt*pt*pt;
                                }
                                else if(harmonicN==1){
                                                return -0.00854808+0.0237419*pt+0.00737678*pt*pt+0.00711372*pt*pt*pt-0.00275382*pt*pt*pt*pt;
                                }
                                else if(harmonicN==2){
                                return 0.0;//0+1.764e-315*pt+2.122e-314*pt*pt+4.68595e-310*pt*pt*pt+1.94786e+118*pt*pt*pt*pt;
                                }
                        }
                }
                
            return 0.0;
}
float TG200::dNdPhi(){
        return (1.0 + 2.0*(v1*cos(phiPart-psi1)+v2*cos(2.0*(phiPart-psi2))+v3*cos(3.0*(phiPart-psi3))+v4*cos(4.0*(phiPart-psi4))))/(2.0*M_PI);
}
int TG200::SpectraCentBin(){
        if(multiplicity>=500){return 0;}
        else if(multiplicity>=345){return 1;}
        else if(multiplicity>=154){return 2;}
        else if(multiplicity>=58){return 3;}
        else{return 6;}
}
int TG200::HarmonicCentBin(){
        if(multiplicity>=500){return 0;}
        else if(multiplicity>=345){return 1;}
        else if(multiplicity>=234 ){return 2;}
        else if(multiplicity>=154){return 3;}
        else if(multiplicity>=99){return 4;}
        else if(multiplicity>=58){return 5;}
        else{return 6;}
}
void TG200::clearEventBuffer(){
        multiplicity = 0;
        myparticles.clear();
        tmpEvent.clear();
        


}
void TG200::clearDistroBuffer(){
        delete ptDistroPip;
        delete ptDistroPim;
        delete ptDistroKap;
        delete ptDistroKam;
        delete ptDistroPro;
        delete ptDistroPbar;
        delete MultiDistro;
}
void TG200::genEvents(){

        tgEvents.clear();
        while(tgEvents.size()!=settings.getNevents()){ 
            
            
            if(doPsi1 ==-1.0) psi1 =  fRandom->Uniform(0.0,2.0*M_PI);
            if(doPsi2 ==-1.0) psi2 =  fRandom->Uniform(0.0,2.0*M_PI);
            if(doPsi3 ==-1.0) psi3 =  fRandom->Uniform(0.0,2.0*M_PI);
            if(doPsi4 ==-1.0) psi4 =  fRandom->Uniform(0.0,2.0*M_PI);

            while(SpectraCentBin() != setCent){   multiplicity = (int)MultiDistro->GetRandom(fRandom);  }

            spectraCent = SpectraCentBin();
            harmonicCent = HarmonicCentBin();
            //multiplicity = int( 2.0*etaRange*( multiplicity+int( fRandom->Uniform(correction_low[spectraCent], correction_high[spectraCent]) ) ) );
            multiplicity = int((etaRange/0.5)*multiplicity);

            partNumbers[0] = int(piplusRatios[spectraCent]*multiplicity);
            partNumbers[1] = int(piminusRatios[spectraCent]*multiplicity);
            partNumbers[2] = int(kaplusRatios[spectraCent]*multiplicity);
            partNumbers[3] = int(kaminusRatios[spectraCent]*multiplicity);
            partNumbers[4] = int(proRatios[spectraCent]*multiplicity);
            partNumbers[5] = int(pbarRatios[spectraCent]*multiplicity);

        
            for(int i=0; i < 6; i++){

                for(int j =0; j<partNumbers[i]; j++){

                    if(i==0){
                        ptPart = ptDistroPip->GetRandom(fRandom);
                        KF = 211;
                    }
                    if(i==1){
                        ptPart = ptDistroPim->GetRandom(fRandom);
                        KF = -211;
                    }


                    if(i==2){
                        ptPart = ptDistroKap->GetRandom(fRandom);
                        KF = 321;                        
                    }
                    if(i==3){
                        ptPart = ptDistroKam->GetRandom(fRandom);
                        KF = -321;
                    }


                    if(i==4){
                        ptPart = ptDistroPro->GetRandom(fRandom);
                        KF = 2212;
                    }
                    if(i==5){
                        ptPart = ptDistroPbar->GetRandom(fRandom);
                        KF = -2212;
                    }

                    etaPart = fRandom->Uniform(-etaRange, etaRange);

                    v2 = doV2*HarmonicFunction(0);
                    v3 = doV3*HarmonicFunction(1);
                    v4 = doV4*HarmonicFunction(2);
                    if(v2<0.0) v2=0.0;
                    if(v3<0.0) v3=0.0;
                    if(v4<0.0) v4=0.0;
                    v1 = doV1*0.02*v2;
                    
                    maxPhi = 2.0*(1.02*abs(v2)+abs(v3)+abs(v4))/M_PI;
                    
                    CHECK =0;
                    while(CHECK!=1){
                            dndphi = fRandom->Uniform(0.0,maxPhi);
                            phiPart = fRandom->Uniform(0,2.0*M_PI);
                            if(dndphi < dNdPhi()) CHECK = 1;
                    }
                    myparticles+= TGParticle(ptPart,etaPart,phiPart,KF);
                }

            }
            tmpEvent = TGEvent(myparticles, spectraCent, multiplicity);
            tmpEvent.Eta(etaRange);
           
           
            tgEvents+= tmpEvent;
            clearEventBuffer();
        

        }
        //myFile->Close();
}
void TG200::genEventsQA(){
        
        string qaString = settings.getOutputDir()+"TG_200GeVQA.root";
        TFile *qa = new TFile(qaString.c_str(),"RECREATE");
        TDirectory *cent0 = qa->mkdir("0-10%");
        TDirectory *cent1 = qa->mkdir("10-20%");
        TDirectory *cent2 = qa->mkdir("20-30%");
        TDirectory *cent3 = qa->mkdir("30-40%");
        TDirectory *cent4 = qa->mkdir("40-50%");
        TDirectory *cent5 = qa->mkdir("50-60%");

       
        TH1D *histpsi_1 = new TH1D("histpsi_1", "#Psi_{EP,1} Thrown distribution for all particles",100,-0.5*M_PI,2.5*M_PI); 
        histpsi_1 ->SetXTitle("#Psi_{EP,1} (radians)");
        histpsi_1 ->SetYTitle("dN/d#Psi_{EP,1}");

        TH1D *histpsi_2 = new TH1D("histpsi_2", "#Psi_{EP,1} Thrown distribution for all particles",100,-0.5*M_PI,2.5*M_PI); 
        histpsi_2 ->SetXTitle("#Psi_{EP,2} (radians)");
        histpsi_2 ->SetYTitle("dN/d#Psi_{EP,2}");

        TH1D *histpsi_3 = new TH1D("histpsi_3", "#Psi_{EP,1} Thrown distribution for all particles",100,-0.5*M_PI,2.5*M_PI); 
        histpsi_3 ->SetXTitle("#Psi_{EP,3} (radians)");
        histpsi_3 ->SetYTitle("dN/d#Psi_{EP,3}");

        TH1D *histpsi_4 = new TH1D("histpsi_4", "#Psi_{EP,1} Thrown distribution for all particles",100,-0.5*M_PI,2.5*M_PI); 
        histpsi_4 ->SetXTitle("#Psi_{EP,4} (radians)");
        histpsi_4 ->SetYTitle("dN/d#Psi_{EP,4}");


        TH2D *v1_pi_h = new TH2D("v1_pi_h" , "v_{1} vs p_{T} for #pi^{+}, #pi^{-}"  , 200 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_pi_h = new TH2D("v2_pi_h" , "v_{2} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v3_pi_h = new TH2D("v3_pi_h" , "v_{3} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_pi_h = new TH2D("v4_pi_h" , "v_{4} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v1_K_h = new TH2D("v1_K_h" , "v_{1} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_K_h = new TH2D("v2_K_h" , "v_{2} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v3_K_h = new TH2D("v3_K_h" , "v_{3} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_K_h = new TH2D("v4_K_h" , "v_{4} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v1_P_h = new TH2D("v1_P_h" , "v_{1} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_P_h = new TH2D("v2_P_h" , "v_{2} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3); 
        TH2D *v3_P_h = new TH2D("v3_P_h" , "v_{3} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_P_h = new TH2D("v4_P_h" , "v_{4} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);

        TH2D *v1_pi_h1 = new TH2D("v1_pi_h1" , "v_{1} vs p_{T} for #pi^{+}, #pi^{-}"  , 200 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_pi_h1 = new TH2D("v2_pi_h1" , "v_{2} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v3_pi_h1 = new TH2D("v3_pi_h1" , "v_{3} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_pi_h1 = new TH2D("v4_pi_h1" , "v_{4} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v1_K_h1 = new TH2D("v1_K_h1" , "v_{1} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_K_h1 = new TH2D("v2_K_h1" , "v_{2} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v3_K_h1 = new TH2D("v3_K_h1" , "v_{3} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_K_h1 = new TH2D("v4_K_h1" , "v_{4} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v1_P_h1 = new TH2D("v1_P_h1" , "v_{1} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_P_h1 = new TH2D("v2_P_h1" , "v_{2} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3); 
        TH2D *v3_P_h1 = new TH2D("v3_P_h1" , "v_{3} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_P_h1 = new TH2D("v4_P_h1" , "v_{4} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);

        TH2D *v1_pi_h2 = new TH2D("v1_pi_h2" , "v_{1} vs p_{T} for #pi^{+}, #pi^{-}"  , 200 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_pi_h2 = new TH2D("v2_pi_h2" , "v_{2} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v3_pi_h2 = new TH2D("v3_pi_h2" , "v_{3} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_pi_h2 = new TH2D("v4_pi_h2" , "v_{4} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v1_K_h2 = new TH2D("v1_K_h2" , "v_{1} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_K_h2 = new TH2D("v2_K_h2" , "v_{2} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v3_K_h2 = new TH2D("v3_K_h2" , "v_{3} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_K_h2 = new TH2D("v4_K_h2" , "v_{4} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v1_P_h2 = new TH2D("v1_P_h2" , "v_{1} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_P_h2 = new TH2D("v2_P_h2" , "v_{2} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3); 
        TH2D *v3_P_h2 = new TH2D("v3_P_h2" , "v_{3} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_P_h2 = new TH2D("v4_P_h2" , "v_{4} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);

        TH2D *v1_pi_h3 = new TH2D("v1_pi_h3" , "v_{1} vs p_{T} for #pi^{+}, #pi^{-}"  , 200 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_pi_h3 = new TH2D("v2_pi_h3" , "v_{2} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v3_pi_h3 = new TH2D("v3_pi_h3" , "v_{3} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_pi_h3 = new TH2D("v4_pi_h3" , "v_{4} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v1_K_h3 = new TH2D("v1_K_h3" , "v_{1} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_K_h3 = new TH2D("v2_K_h3" , "v_{2} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v3_K_h3 = new TH2D("v3_K_h3" , "v_{3} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_K_h3 = new TH2D("v4_K_h3" , "v_{4} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v1_P_h3 = new TH2D("v1_P_h3" , "v_{1} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_P_h3 = new TH2D("v2_P_h3" , "v_{2} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3); 
        TH2D *v3_P_h3 = new TH2D("v3_P_h3" , "v_{3} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_P_h3 = new TH2D("v4_P_h3" , "v_{4} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);

        TH2D *v1_pi_h4 = new TH2D("v1_pi_h4" , "v_{1} vs p_{T} for #pi^{+}, #pi^{-}"  , 200 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_pi_h4 = new TH2D("v2_pi_h4" , "v_{2} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v3_pi_h4 = new TH2D("v3_pi_h4" , "v_{3} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_pi_h4 = new TH2D("v4_pi_h4" , "v_{4} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v1_K_h4 = new TH2D("v1_K_h4" , "v_{1} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_K_h4 = new TH2D("v2_K_h4" , "v_{2} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v3_K_h4 = new TH2D("v3_K_h4" , "v_{3} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_K_h4 = new TH2D("v4_K_h4" , "v_{4} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v1_P_h4 = new TH2D("v1_P_h4" , "v_{1} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_P_h4 = new TH2D("v2_P_h4" , "v_{2} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3); 
        TH2D *v3_P_h4 = new TH2D("v3_P_h4" , "v_{3} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_P_h4 = new TH2D("v4_P_h4" , "v_{4} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);

        TH2D *v1_pi_h5 = new TH2D("v1_pi_h5" , "v_{1} vs p_{T} for #pi^{+}, #pi^{-}"  , 200 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_pi_h5 = new TH2D("v2_pi_h5" , "v_{2} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v3_pi_h5 = new TH2D("v3_pi_h5" , "v_{3} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_pi_h5 = new TH2D("v4_pi_h5" , "v_{4} vs p_{T} for #pi^{+}, #pi^{-}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v1_K_h5 = new TH2D("v1_K_h5" , "v_{1} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_K_h5 = new TH2D("v2_K_h5" , "v_{2} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v3_K_h5 = new TH2D("v3_K_h5" , "v_{3} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_K_h5 = new TH2D("v4_K_h5" , "v_{4} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v1_P_h5 = new TH2D("v1_P_h5" , "v_{1} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v2_P_h5 = new TH2D("v2_P_h5" , "v_{2} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3); 
        TH2D *v3_P_h5 = new TH2D("v3_P_h5" , "v_{3} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
        TH2D *v4_P_h5 = new TH2D("v4_P_h5" , "v_{4} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);

        nEvent = 25000;
       

        TH1D* multiplicity_all_h = new TH1D("multiplicity_all_h","multiplicity_all_h",1800,0,1800);
        TH1D *multiplicity_cent0_h = new TH1D("multiplicity_cent0_h" , "multiplicity distribution for cent 0" ,1800 , 0 , 1800);
        TH1D *multiplicity_cent1_h = new TH1D("multiplicity_cent1_h" , "multiplicity distribution for cent 1" ,1800 , 0 , 1800);
        TH1D *multiplicity_cent2_h = new TH1D("multiplicity_cent2_h" , "multiplicity distribution for cent 2" ,1800 , 0 , 1800);
        TH1D *multiplicity_cent3_h = new TH1D("multiplicity_cent3_h" , "multiplicity distribution for cent 3" ,1800 , 0 , 1800);


        for(int i=0; i< 4; i++){
            setCent = i;
            getRootDistros();
            for(int j =0;j<nEvent; j++){

                if(doPsi1 ==-1.0) psi1 =  fRandom->Uniform(0.0,2.0*M_PI);
                if(doPsi2 ==-1.0) psi2 =  fRandom->Uniform(0.0,2.0*M_PI);
                if(doPsi3 ==-1.0) psi3 =  fRandom->Uniform(0.0,2.0*M_PI);
                if(doPsi4 ==-1.0) psi4 =  fRandom->Uniform(0.0,2.0*M_PI);

                histpsi_1->Fill(psi1);
                histpsi_2->Fill(psi2);
                histpsi_3->Fill(psi3);
                histpsi_4->Fill(psi4);
                
                while(SpectraCentBin() != setCent){   multiplicity = (int)MultiDistro->GetRandom(fRandom);  }
                spectraCent = SpectraCentBin();
                harmonicCent = HarmonicCentBin();
                 //multiplicity = int( 2.0*etaRange*( multiplicity+int( fRandom->Uniform(correction_low[spectraCent], correction_high[spectraCent]) ) ) );
                multiplicity = int((etaRange/0.5)*multiplicity);
                
                multiplicity_all_h->Fill(multiplicity);
                if(setCent == 0) multiplicity_cent0_h->Fill(multiplicity);
                if(setCent == 1) multiplicity_cent1_h->Fill(multiplicity);
                if(setCent == 2) multiplicity_cent2_h->Fill(multiplicity);
                if(setCent == 3) multiplicity_cent3_h->Fill(multiplicity);

                
                partNumbers[0] = int(piplusRatios[spectraCent]*multiplicity);
                partNumbers[1] = int(piminusRatios[spectraCent]*multiplicity);
                partNumbers[2] = int(kaplusRatios[spectraCent]*multiplicity);
                partNumbers[3] = int(kaminusRatios[spectraCent]*multiplicity);
                partNumbers[4] = int(proRatios[spectraCent]*multiplicity);
                partNumbers[5] = int(pbarRatios[spectraCent]*multiplicity);

                for(int k=0; k < 6; k++){

                        for(int l =0; l<partNumbers[k]; l++){

                                if(k==0){
                                        ptPart = ptDistroPip->GetRandom(fRandom);
                                        KF = 211;
        
                                }
                                if(k==1){
                                        ptPart = ptDistroPim->GetRandom(fRandom);
                                        KF = -211;
                                }


                                if(k==2){
                                        ptPart = ptDistroKap->GetRandom(fRandom);
                                        KF = 321;                        
                                }
                                if(k==3){
                                        ptPart = ptDistroKam->GetRandom(fRandom);
                                        KF = -321;
                                }


                                if(k==4){
                                        ptPart = ptDistroPro->GetRandom(fRandom);
                                        KF = 2212;
                                }
                                if(k==5){
                                        ptPart = ptDistroPbar->GetRandom(fRandom);
                                        KF = -2212;
                                }
                                
                                v2 = doV2*HarmonicFunction(0);
                                v3 = doV3*HarmonicFunction(1);
                                v4 = doV4*HarmonicFunction(2);
                                if(v2<0.0) v2=0.0;
                                if(v3<0.0) v3=0.0;
                                if(v4<0.0) v4=0.0;
                                v1 = doV1*0.02*v2;

                                if(KF==211||KF==-211){

                                        if(harmonicCent ==0) v1_pi_h->Fill(ptPart,v1);
                                        else if(harmonicCent ==1) v1_pi_h1->Fill(ptPart,v1);
                                        else if(harmonicCent ==2) v1_pi_h2->Fill(ptPart,v1);
                                        else if(harmonicCent ==3) v1_pi_h3->Fill(ptPart,v1);
                                        else if(harmonicCent ==4) v1_pi_h4->Fill(ptPart,v1);
                                        else if(harmonicCent ==5) v1_pi_h5->Fill(ptPart,v1);

                                        if(harmonicCent ==0) v2_pi_h->Fill(ptPart,v2);
                                        else if(harmonicCent ==1) v2_pi_h1->Fill(ptPart,v2);
                                        else if(harmonicCent ==2) v2_pi_h2->Fill(ptPart,v2);
                                        else if(harmonicCent ==3) v2_pi_h3->Fill(ptPart,v2);
                                        else if(harmonicCent ==4) v2_pi_h4->Fill(ptPart,v2);
                                        else if(harmonicCent ==5) v2_pi_h5->Fill(ptPart,v2);

                                        if(harmonicCent ==0) v3_pi_h->Fill(ptPart,v3);
                                        else if(harmonicCent ==1) v3_pi_h1->Fill(ptPart,v3);
                                        else if(harmonicCent ==2) v3_pi_h2->Fill(ptPart,v3);
                                        else if(harmonicCent ==3) v3_pi_h3->Fill(ptPart,v3);
                                        else if(harmonicCent ==4) v3_pi_h4->Fill(ptPart,v3);
                                        else if(harmonicCent ==5) v3_pi_h5->Fill(ptPart,v3);

                                        if(harmonicCent ==0) v4_pi_h->Fill(ptPart,v4);
                                        else if(harmonicCent ==1) v4_pi_h1->Fill(ptPart,v4);
                                        else if(harmonicCent ==2) v4_pi_h2->Fill(ptPart,v4);
                                        else if(harmonicCent ==3) v4_pi_h3->Fill(ptPart,v4);
                                        else if(harmonicCent ==4) v4_pi_h4->Fill(ptPart,v4);
                                        else if(harmonicCent ==5) v4_pi_h5->Fill(ptPart,v4);

                                }
                                if(KF==321||KF==-321){
                                        if(harmonicCent ==0) v1_K_h->Fill(ptPart,v1);
                                        else if(harmonicCent ==1) v1_K_h1->Fill(ptPart,v1);
                                        else if(harmonicCent ==2) v1_K_h2->Fill(ptPart,v1);
                                        else if(harmonicCent ==3) v1_K_h3->Fill(ptPart,v1);
                                        else if(harmonicCent ==4) v1_K_h4->Fill(ptPart,v1);
                                        else if(harmonicCent ==5) v1_K_h5->Fill(ptPart,v1);

                                        if(harmonicCent ==0) v2_K_h->Fill(ptPart,v2);
                                        else if(harmonicCent ==1) v2_K_h1->Fill(ptPart,v2);
                                        else if(harmonicCent ==2) v2_K_h2->Fill(ptPart,v2);
                                        else if(harmonicCent ==3) v2_K_h3->Fill(ptPart,v2);
                                        else if(harmonicCent ==4) v2_K_h4->Fill(ptPart,v2);
                                        else if(harmonicCent ==5) v2_K_h5->Fill(ptPart,v2);

                                        if(harmonicCent ==0) v3_K_h->Fill(ptPart,v3);
                                        else if(harmonicCent ==1) v3_K_h1->Fill(ptPart,v3);
                                        else if(harmonicCent ==2) v3_K_h2->Fill(ptPart,v3);
                                        else if(harmonicCent ==3) v3_K_h3->Fill(ptPart,v3);
                                        else if(harmonicCent ==4) v3_K_h4->Fill(ptPart,v3);
                                        else if(harmonicCent ==5) v3_K_h5->Fill(ptPart,v3);

                                        if(harmonicCent ==0) v4_K_h->Fill(ptPart,v4);
                                        else if(harmonicCent ==1) v4_K_h1->Fill(ptPart,v4);
                                        else if(harmonicCent ==2) v4_K_h2->Fill(ptPart,v4);
                                        else if(harmonicCent ==3) v4_K_h3->Fill(ptPart,v4);
                                        else if(harmonicCent ==4) v4_K_h4->Fill(ptPart,v4);
                                        else if(harmonicCent ==5) v4_K_h5->Fill(ptPart,v4);

                                }
                                if(KF==2212||KF==-2212){
                                        
                                        if(harmonicCent ==0) v1_P_h->Fill(ptPart,v1);
                                        else if(harmonicCent ==1) v1_P_h1->Fill(ptPart,v1);
                                        else if(harmonicCent ==2) v1_P_h2->Fill(ptPart,v1);
                                        else if(harmonicCent ==3) v1_P_h3->Fill(ptPart,v1);
                                        else if(harmonicCent ==4) v1_P_h4->Fill(ptPart,v1);
                                        else if(harmonicCent ==5) v1_P_h5->Fill(ptPart,v1);

                                        if(harmonicCent ==0) v2_P_h->Fill(ptPart,v2);
                                        else if(harmonicCent ==1) v2_P_h1->Fill(ptPart,v2);
                                        else if(harmonicCent ==2) v2_P_h2->Fill(ptPart,v2);
                                        else if(harmonicCent ==3) v2_P_h3->Fill(ptPart,v2);
                                        else if(harmonicCent ==4) v2_P_h4->Fill(ptPart,v2);
                                        else if(harmonicCent ==5) v2_P_h5->Fill(ptPart,v2);

                                        if(harmonicCent ==0) v3_P_h->Fill(ptPart,v3);
                                        else if(harmonicCent ==1) v3_P_h1->Fill(ptPart,v3);
                                        else if(harmonicCent ==2) v3_P_h2->Fill(ptPart,v3);
                                        else if(harmonicCent ==3) v3_P_h3->Fill(ptPart,v3);
                                        else if(harmonicCent ==4) v3_P_h4->Fill(ptPart,v3);
                                        else if(harmonicCent ==5) v3_P_h5->Fill(ptPart,v3);

                                        if(harmonicCent ==0) v4_P_h->Fill(ptPart,v4);
                                        else if(harmonicCent ==1) v4_P_h1->Fill(ptPart,v4);
                                        else if(harmonicCent ==2) v4_P_h2->Fill(ptPart,v4);
                                        else if(harmonicCent ==3) v4_P_h3->Fill(ptPart,v4);
                                        else if(harmonicCent ==4) v4_P_h4->Fill(ptPart,v4);
                                        else if(harmonicCent ==5) v4_P_h5->Fill(ptPart,v4);
                                
                                
                                }


                        
                        }

                }

                clearEventBuffer();



            }
           clearDistroBuffer();
        }

        // multiplicity_cent0_h->Scale(1./multiplicity_cent0_h->Integral());
        // multiplicity_cent1_h->Scale(1./multiplicity_cent1_h->Integral());
        // multiplicity_cent2_h->Scale(1./multiplicity_cent2_h->Integral());
        // multiplicity_cent3_h->Scale(1./multiplicity_cent3_h->Integral());
        multiplicity_all_h->Scale(1./multiplicity_all_h->Integral());

        multiplicity_cent0_h->SetLineColor(kRed);
        multiplicity_cent1_h->SetLineColor(kBlue);
        multiplicity_cent2_h->SetLineColor(kGreen);
        multiplicity_cent3_h->SetLineColor(kCyan);


        TCanvas *c2 = new TCanvas("c2","c2",800,600);
        c2->cd();
        gPad->SetLogy();
        gPad->SetGrid();
        multiplicity_all_h->Draw("HIST");
        //multiplicity_all_h->Rebin(1);
        multiplicity_all_h->SetLineWidth(2);
        multiplicity_all_h->SetLineColor(kBlack);
        multiplicity_all_h->SetFillColorAlpha(kBlack,0.25);
        multiplicity_all_h->SetFillStyle(3001);

        multiplicity_all_h->GetXaxis()->SetTitle("N_{ch}");
        multiplicity_all_h->GetYaxis()->SetTitle("#frac{1}{N} dN/dN_{ch}^{raw}");
        multiplicity_all_h->GetYaxis()->SetTitleOffset(1.5);
        multiplicity_all_h->GetYaxis()->SetRangeUser(0.0001,1);
        multiplicity_all_h->GetXaxis()->SetRangeUser(0,1800);
        multiplicity_all_h->GetXaxis()->SetLabelSize(0.03);
        multiplicity_all_h->GetYaxis()->SetLabelSize(0.03);
        multiplicity_all_h->GetXaxis()->SetTitleSize(0.03);
        multiplicity_all_h->GetYaxis()->SetTitleSize(0.03);
        multiplicity_all_h->GetXaxis()->SetTitleOffset(1.1);
        multiplicity_all_h->GetYaxis()->SetTitleOffset(1.1);
        multiplicity_all_h->SetTitle("");
        multiplicity_all_h->GetXaxis()->CenterTitle();

        multiplicity_cent0_h->Draw("HIST same");
        multiplicity_cent0_h->SetFillColorAlpha(kRed,0.5);
        multiplicity_cent0_h->SetFillStyle(3001);
        multiplicity_cent1_h->Draw("HIST same");
        multiplicity_cent1_h->SetFillColorAlpha(kBlue,0.5);
        multiplicity_cent1_h->SetFillStyle(3001);
        multiplicity_cent2_h->Draw("HIST same");
        multiplicity_cent2_h->SetFillColorAlpha(kGreen,0.5);
        multiplicity_cent2_h->SetFillStyle(3001);
        multiplicity_cent3_h->Draw("HIST same");
        multiplicity_cent3_h->SetFillColorAlpha(kCyan,0.5);
        multiplicity_cent3_h->SetFillStyle(3001);

        TLegend *legend = new TLegend(0.55,0.65,0.85,0.85);
        legend->SetBorderSize(0);
        legend->SetFillColor(0);

        legend->AddEntry(multiplicity_all_h,"STAR 200 GeV Au+Au N_{ch}^{raw}","lf");
        legend->AddEntry(multiplicity_cent0_h,"0-10%","lf");
        legend->AddEntry(multiplicity_cent1_h,"10-20%","lf");
        legend->AddEntry(multiplicity_cent2_h,"20-40%","lf");
        legend->AddEntry(multiplicity_cent3_h,"40-60%","lf");
        legend->Draw();
        c2->SetTitle("Multiplicity");
        c2->Write();

        histpsi_1->Write();
        histpsi_2->Write();
        histpsi_3->Write();
        histpsi_4->Write();

        cent0->cd();
        v1_pi_h->Write(); 
        v2_pi_h->Write(); 
        v3_pi_h->Write(); 
        v4_pi_h->Write(); 
        v1_K_h->Write(); 
        v2_K_h->Write(); 
        v3_K_h->Write(); 
        v4_K_h->Write(); 
        v1_P_h->Write(); 
        v2_P_h->Write();  
        v3_P_h->Write(); 
        v4_P_h->Write(); 
        cent1->cd();
        v1_pi_h1->Write(); 
        v2_pi_h1->Write(); 
        v3_pi_h1->Write(); 
        v4_pi_h1->Write(); 
        v1_K_h1->Write(); 
        v2_K_h1->Write(); 
        v3_K_h1->Write(); 
        v4_K_h1->Write(); 
        v1_P_h1->Write(); 
        v2_P_h1->Write();  
        v3_P_h1->Write(); 
        v4_P_h1->Write(); 
        cent2->cd();
        v1_pi_h2->Write(); 
        v2_pi_h2->Write(); 
        v3_pi_h2->Write(); 
        v4_pi_h2->Write(); 
        v1_K_h2->Write(); 
        v2_K_h2->Write(); 
        v3_K_h2->Write(); 
        v4_K_h2->Write(); 
        v1_P_h2->Write(); 
        v2_P_h2->Write();  
        v3_P_h2->Write(); 
        v4_P_h2->Write(); 
        cent3->cd();
        v1_pi_h3->Write(); 
        v2_pi_h3->Write(); 
        v3_pi_h3->Write(); 
        v4_pi_h3->Write(); 
        v1_K_h3->Write(); 
        v2_K_h3->Write(); 
        v3_K_h3->Write(); 
        v4_K_h3->Write(); 
        v1_P_h3->Write(); 
        v2_P_h3->Write();  
        v3_P_h3->Write(); 
        v4_P_h3->Write(); 
        cent4->cd();
        v1_pi_h4->Write(); 
        v2_pi_h4->Write(); 
        v3_pi_h4->Write(); 
        v4_pi_h4->Write(); 
        v1_K_h4->Write(); 
        v2_K_h4->Write(); 
        v3_K_h4->Write(); 
        v4_K_h4->Write(); 
        v1_P_h4->Write(); 
        v2_P_h4->Write();  
        v3_P_h4->Write(); 
        v4_P_h4->Write(); 
        cent5->cd();
        v1_pi_h5->Write(); 
        v2_pi_h5->Write(); 
        v3_pi_h5->Write(); 
        v4_pi_h5->Write(); 
        v1_K_h5->Write(); 
        v2_K_h5->Write(); 
        v3_K_h5->Write(); 
        v4_K_h5->Write(); 
        v1_P_h5->Write(); 
        v2_P_h5->Write();  
        v3_P_h5->Write(); 
        v4_P_h5->Write(); 
        
        myFile->Close();
        qa->Close();    
               


}
void TG200::getRootDistros(){

        ptHistoPip = "pi+_pt_cent"+std::to_string(setCent);
        ptHistoPim = "pi-_pt_cent"+std::to_string(setCent);
        ptHistoKap = "ka+_pt_cent"+std::to_string(setCent);
        ptHistoKam = "ka-_pt_cent"+std::to_string(setCent);
        ptHistoPro = "pro_pt_cent"+std::to_string(setCent);
        ptHistoPbar = "pba_pt_cent"+std::to_string(setCent);

        ptDistroPip = (TH1F*)myFile->Get(ptHistoPip.c_str());
        ptDistroPim = (TH1F*)myFile->Get(ptHistoPim.c_str());
        ptDistroKap = (TH1F*)myFile->Get(ptHistoKap.c_str());
        ptDistroKam = (TH1F*)myFile->Get(ptHistoKam.c_str());
        ptDistroPro = (TH1F*)myFile->Get(ptHistoPro.c_str());
        ptDistroPbar = (TH1F*)myFile->Get(ptHistoPbar.c_str());


        TH1F* All_Multiplicity_Distro = (TH1F*)myFile->Get("MultiplicityDistro");
        TH1F *multiplicity_cent_h = new TH1F("multiplicity_cent_h" , "multiplicity distribution for cent 0" ,800 , 0 , 800);

        int raw_temp_low[9] = {14,30,56,94,145,217,312,431,510};
        int raw_temp_high[9]={29,55,93,144,216,311,430,509,608};
        float raw_temp_mean[9]= {22.5,43.1,74.8,120,181,264,370,470,559};
        int corr_mean[9]= {22,45,78,126,195,287,421,558,691};
        float corr_err[9] = {2,3,6,9,14,20,30,40,49};

        float nbin_raw[9];
        float diff_mean[9];
        float raw_error[9];

        for(int i = 0;i<9;i++){
                nbin_raw[i] = raw_temp_high[i]-raw_temp_low[i];
                diff_mean[i] = corr_mean[i]-raw_temp_mean[i];
                raw_error[i] = nbin_raw[i]/2;
        }

        TGraphErrors *mygraph = new TGraphErrors(9,raw_temp_mean,diff_mean,raw_error,corr_err);
        TF1* fit_pol = new TF1("fit_pol","pol3",0,600);
        mygraph->Fit(fit_pol,"BRQ");
        float p0 = fit_pol->GetParameter(0);
        float p1 = fit_pol->GetParameter(1);
        float p2 = fit_pol->GetParameter(2);
        float p3 = fit_pol->GetParameter(3);

        TF1* fit_final = new TF1("fit_final","[0]+[1]*x+[2]*x**2+[3]*x**3",0,800);
        fit_final->SetParameters(p0,p1,p2,p3);
        fit_final->FixParameter(0,p0);
        fit_final->FixParameter(1,p1);
        fit_final->FixParameter(2,p2);
        fit_final->FixParameter(3,p3);

        for(int i=0;i<1000000;i++){
                float multtemp = All_Multiplicity_Distro->GetRandom();
                multtemp+=fit_final->Eval(multtemp);
                multiplicity_cent_h->Fill(multtemp);
        }
        multiplicity_cent_h->Scale(1/multiplicity_cent_h->Integral());
        MultiDistro= (TH1F*)multiplicity_cent_h->Clone("MultiplicityDistro");
       // MultiDistro = (TH1F*)myFile->Get("MultiplicityDistro");


}
void TG200::config(TGSettings inputSettings) {

        settings = inputSettings;
        tgEvents.clear();
        clearEventBuffer();
       // clearDistroBuffer();

        doPsi1 = 0;
        doPsi2 = 0;
        doPsi3 = 0;
        doPsi4 = 0;

        nEvent = settings.getNevents();
        etaRange = settings.getEtaRange();
        setCent = settings.getCentBin();
        doV1 = settings.getVN(1);
        doV2 = settings.getVN(2);
        doV3 = settings.getVN(3);
        doV4 = settings.getVN(4);

        psi1 = settings.getPsiN(1);
        if(psi1 ==-1.0 ) doPsi1 = -1.0;
        psi2 = settings.getPsiN(2);
        if(psi2 ==-1.0 ) doPsi2 = -1.0;
        psi3 = settings.getPsiN(3);
        if(psi3 ==-1.0 ) doPsi3 = -1.0;
        psi4 = settings.getPsiN(4);
        if(psi4 ==-1.0 ) doPsi4 = -1.0;

        getRootDistros();


}
void TG200::seed(TRandom3* FRandom){
        fRandom = FRandom;

}
TENNGEN_END_NAMESPACE