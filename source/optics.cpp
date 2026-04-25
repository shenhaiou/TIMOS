#include "optics.h"
#include "io_utils.h"
#include <iostream>
#include <fstream>

using namespace std;

TiResult ReadOpticalParameter(const std::string& filename, SimContext& ctx){
  ifstream fin(filename);
  if(!fin.good()){ cout << "\tCould not open file: " << filename << "\n"; return std::unexpected("error"); }

  int type;
  skip_comments(fin); fin >> type;
  if     (type==1){ cout<<"\tOptical parameter per region\n";  ctx.simpleOptic=true;  }
  else if(type==2){ cout<<"\tOptical parameter per element\n"; ctx.simpleOptic=false; }

  skip_comments(fin); fin >> ctx.numMed;
  ctx.medOptic.resize(ctx.numMed+1);

  for(int i = 1; i <= ctx.numMed; i++){
    skip_comments(fin);
    fin >> ctx.medOptic[i].mua >> ctx.medOptic[i].mus
        >> ctx.medOptic[i].g   >> ctx.medOptic[i].RefIdx;

    if(ctx.medOptic[i].mua<0 || ctx.medOptic[i].mus<0){
      cout<<"\tmua and mus must be >= 0.\n"; return std::unexpected("error");
    }
    if(ctx.medOptic[i].g<0 || ctx.medOptic[i].g>1){
      cout<<"\tg must be in [0,1].\n"; return std::unexpected("error");
    }
    if(ctx.medOptic[i].RefIdx<1){
      cout<<"\tRefractive index must be >= 1.\n"; return std::unexpected("error");
    }
    if(ctx.medOptic[i].mua<1e-10 && ctx.medOptic[i].mus<1e-10){
      cout<<"\tmua and mus cannot both be zero.\n"; return std::unexpected("error");
    }
    double g = ctx.medOptic[i].g;
    ctx.medOptic[i].OneMinsGG = 1.0-g*g;
    ctx.medOptic[i].OneMinsG  = 1.0-g;
    ctx.medOptic[i].OneAddGG  = 1.0+g*g;
    ctx.medOptic[i].OneAddG   = 1.0+g;
    ctx.medOptic[i].TwoG      = 2.0*g;
    ctx.medOptic[i].MUAMUS    = ctx.medOptic[i].mua + ctx.medOptic[i].mus;
    ctx.medOptic[i].IMUAMUS   = 1.0/ctx.medOptic[i].MUAMUS;
    ctx.medOptic[i].pdwa      = ctx.medOptic[i].mua/ctx.medOptic[i].MUAMUS;
  }

  skip_comments(fin); fin >> type;
  if(type==1){
    cout<<"\tUniform environment refractive index.\n";
    ctx.uniformBoundary = 0;
    skip_comments(fin); fin >> ctx.envRefIdx;
  }else if(type==2){
    cout<<"\tMatched boundary (no reflection).\n";
    ctx.uniformBoundary = 1;
  }
  return {};
}
