#include <iostream>
#include <string>

#include "TennGen.h"

using namespace tenngen;


int main(){

    TennGen tg;

    tg.setcollen(200);
    tg.setnevent(1000);
    tg.setcent(0); // 0-10%
    tg.seteta(1.1); // abs value between 0 and 1.1

    tg.setvN(1,true);
    tg.setvN(2,true);
    tg.setvN(3,true);
    tg.setvN(4,true);
    
    tg.setpsiN(1,-1.0); // -1.0 means it will be randomly choosen between 0-2 pi
    tg.setpsiN(2,0.0);
    tg.setpsiN(3,-1.0);
    tg.setpsiN(4,0.0);
    tg.setOutputDir("output-200");

    tg.do_Histos(true); // will make debug root plots
    tg.do_TTree(true); // will make TTree 
    tg.set_Batch(false);
    tg.set_QA(false);
    tg.set_Stream(true);
    tg.init();


    // all of above tg can be achieve via
    //tg.defaultSettings200(true);

     
    cout << "Multiplicity |Eta| centralitiy \n";
    for(int j =0; j < 100000; j++){ cout << j << "\t" << tg.next() << endl;  }


    return 0;

}