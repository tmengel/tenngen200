#include "TennGenFunctions.h"

int main(int argc, char *argv[]){

    std::string rootfilename = "AuAu200Data.root";
    rootfilename = argv[1];
    tenngen::makeDistros200(rootfilename.c_str());
    return 0;

}