#include "optics.h"
#include "io_utils.h"
#include <iostream>
#include <fstream>

using namespace std;

int ReadOpticalParameter(const std::string& filename,
                         bool&       simpleOptic,
                         int&        numMed,
                         int&        uniformBoundary,
                         double&     envRefIdx,
                         TMedOptic*& medOptic){
  ifstream fin(filename);
  if(!fin.good()){
    cerr << "\tCould not open file: " << filename << endl;
    return -1;
  }

  int type;
  skip_comments(fin); fin >> type;
  if     (type == 1){ cerr << "\tOptical parameter per region\n";  simpleOptic = true;  }
  else if(type == 2){ cerr << "\tOptical parameter per element\n"; simpleOptic = false; }

  skip_comments(fin); fin >> numMed;
  medOptic = new TMedOptic[numMed+1];

  for(int i = 1; i <= numMed; i++){
    skip_comments(fin);
    fin >> medOptic[i].mua >> medOptic[i].mus >> medOptic[i].g >> medOptic[i].RefIdx;

    if(medOptic[i].mua < 0 || medOptic[i].mus < 0){
      cerr << "\tmua and mus must be >= 0.\n"; return -1;
    }
    if(medOptic[i].g < 0 || medOptic[i].g > 1){
      cerr << "\tg must be in [0, 1].\n"; return -1;
    }
    if(medOptic[i].RefIdx < 1){
      cerr << "\tRefractive index must be >= 1.\n"; return -1;
    }
    if(medOptic[i].mua < 1e-10 && medOptic[i].mus < 1e-10){
      cerr << "\tmua and mus cannot both be zero. "
              "For glass, set g=1 and mua/mus to small non-zero values.\n";
      return -1;
    }

    double g = medOptic[i].g;
    medOptic[i].OneMinsGG = 1.0 - g*g;
    medOptic[i].OneMinsG  = 1.0 - g;
    medOptic[i].OneAddGG  = 1.0 + g*g;
    medOptic[i].OneAddG   = 1.0 + g;
    medOptic[i].TwoG      = 2.0 * g;
    medOptic[i].MUAMUS    = medOptic[i].mua + medOptic[i].mus;
    medOptic[i].IMUAMUS   = 1.0 / medOptic[i].MUAMUS;
    medOptic[i].pdwa      = medOptic[i].mua / medOptic[i].MUAMUS;
  }

  skip_comments(fin); fin >> type;
  if(type == 1){
    cerr << "\tUniform environment refractive index.\n";
    uniformBoundary = 0;
    skip_comments(fin); fin >> envRefIdx;
  }else if(type == 2){
    cerr << "\tMatched boundary (no reflection).\n";
    uniformBoundary = 1;
  }
  return 0;
}
