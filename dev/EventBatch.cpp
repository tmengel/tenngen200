#include <iostream>
#include <string>
#include <chrono>
#include <vector>
#include <iomanip>


// #include "TannGenEmbed.h"
//#include "PythiaGen.h"
//#include "ptPythiaTest.h"
//#include "TannGenPythia.h"
#include "TannGen.h"

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[]){



        TRandom3* fRandom= new TRandom3(3);


        int nEvent = 1000;
        int cent = 0;
        float etaRange = 1.1;

        int v1 =1;
        int v2 =1;
        int v3 =1;
        int v4 =1;

        int numBatches = 22;
        double pTlimitLow[23] = {5.,7.,8.,10.,11.,12.,15.,16.,17.0,20.,21.,23.,25.,27.,30.,35.,40.,45.,50.,55.,60.,65.,-1.};


        double R = 0.4;
        double ptLead = 5;

        if(argc >=2) nEvent = atoi(argv[1]);
        if(argc >=3) cent = atoi(argv[2]);
        if(argc >=4) etaRange = 1.0*atoi(argv[3]);

        if(argc >=5) v1 = atoi(argv[4]);
        if(argc >=6) v2 = atoi(argv[5]);
        if(argc >=7) v3 = atoi(argv[6]);
        if(argc >=8) v4 = atoi(argv[7]);

        string filedir = "pythiaData";

        auto start = high_resolution_clock::now();
        std::vector<TGEvent> myEvents =  TannGen(nEvent,fRandom,cent,etaRange);
        

        
        //ptPythiaTest(testfile);
        // // std::vector<TGEvent> PythiaGen;
        for(int i = 0; i<15; i++){
            
        }
         
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        cout << "=======================================================================\n";
        cout <<setprecision(6) <<fixed;
        cout <<"Events:  " << numBatches*nEvent << " Total time: " << duration.count()/1000000 << " seconds " << "Average: " << duration.count()/(numBatches*nEvent) << " microseconds/event" << endl;
        cout << "=======================================================================\n";
        return 0;

}



