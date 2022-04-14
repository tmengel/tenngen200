#include <iostream>
#include <string>

#include "TennGen.h"

using namespace tenngen;


int main(){

    TennGen tg;

    tg.setcollen(200);
    tg.setnevent(100);
    tg.setcent(1); // 0-10%
    tg.seteta(1.1); // abs value between 0 and 1.1

    tg.setvN(1,true);
    tg.setvN(2,true);
    tg.setvN(3,true);
    tg.setvN(4,true);
    
    tg.setpsiN(1,-1.0); // -1.0 means it will be randomly choosen between 0-2 pi
    tg.setpsiN(2,0.0);
    tg.setpsiN(3,-1.0);
    tg.setpsiN(4,0.0);
    tg.do_QA(true);
    

   // tg.do_Histos(false); // will make debug root plots
   // tg.do_TTree(false); // will make TTree 

    // all of above tg can be achieve via
    //tg.defaultSettings200(true);


    tg.setOutputDir("output-200");
    tg.init();

    // cout << "Multiplicity |Eta| centralitiy \n";
    // cout << tg[999] << endl; 
    // cout << "particles: pt eta phi kf" <<endl;
    // for(int j =0; j < tg[999].size(); j++){ cout << "\t" << tg[999][j] << endl;  }


    return 0;

}