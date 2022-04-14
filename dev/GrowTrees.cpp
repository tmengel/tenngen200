#include <iostream>
#include <string>
#include <chrono>
#include <vector>
#include <iomanip>

#include "PythiaGen.h"

int main(int argc, char *argv[]){


    int nEvent = 1000000;
    int cent = 0;
    float etaRange = 1.1;
    int v1 =1;
    int v2 =1;
    int v3 =1;
    int v4 =1;


    if(argc >=2) nEvent = atoi(argv[1]);
    if(argc >=3) cent = atoi(argv[2]);
    if(argc >=4) etaRange = 1.0*atoi(argv[3]);
    if(argc >=5) v1 = atoi(argv[4]);
    if(argc >=6) v2 = atoi(argv[5]);
    if(argc >=7) v3 = atoi(argv[6]);
    if(argc >=8) v4 = atoi(argv[7]);

    int numBatches = 22;
    double pTlimitLow[23] = {5.,7.,8.,10.,11.,12.,15.,16.,17.0,20.,21.,23.,25.,27.,30.,35.,40.,45.,50.,55.,60.,65.,-1.};

    string OutDir = "/home/tmengel/tenngen200/src/output/200/pythia";

    for(int i =0; i<(numBatches-1);i++) PythiaGen(nEvent, pTlimitLow[i], pTlimitLow[i+1],i,OutDir);
  

    return 0;

}