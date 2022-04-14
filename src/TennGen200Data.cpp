#include "TennGenFunctions.h"
using namespace tenngen;
int main(){

    std::string rootfilename = "AuAu200Data.root";
    tenngen::makeDistros200(rootfilename.c_str());
    return 0;

}