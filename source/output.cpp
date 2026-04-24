#include "output.h"
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

void AbsorptionToFluence(SimContext& ctx,
                         double& sufThreshhold, double& intThreshhold){
  double* CAbs=new double[ctx.numElem+1]();
  double total=0, maxSuf=0, maxInt=0;
  for(int i=1;i<=ctx.numBoundaryTrig;i++){
    ctx.surfMeas[i]/=ctx.triNodes[ctx.boundaryTrigs[i]].Area;
    if(ctx.surfMeas[i]>maxSuf) maxSuf=ctx.surfMeas[i];
  }
  for(int i=1;i<=ctx.numElem;i++){
    if(ctx.medOptic[ctx.elems[i].MedIdx].mua==0.0){ ctx.absorption[i]=-1.0; }
    else{
#ifdef NOPUREGLASS
      CAbs[i]=ctx.absorption[i];
#else
      if(ctx.medOptic[ctx.elems[i].MedIdx].g<1.0) CAbs[i]=ctx.absorption[i];
#endif
      ctx.absorption[i]=(ctx.absorption[i]/ctx.elemNodes[i].Vol)/ctx.medOptic[ctx.elems[i].MedIdx].mua;
      if(ctx.absorption[i]>maxInt) maxInt=ctx.absorption[i];
    }
  }
  sort(CAbs,CAbs+ctx.numElem+1);
  for(int i=0;i<=ctx.numElem;i++) total+=CAbs[i];
  cerr<<"Absorbed Fraction: "<<total/double(ctx.totalPhoton)<<"\n";
  sufThreshhold=maxSuf/1e20; intThreshhold=maxInt/1e20;
  delete[] CAbs;
}

void TimeAbsorptionToFluence(SimContext& ctx,
                             double& sufThreshhold, double& intThreshhold){
  double* CAbs=new double[ctx.numElem+1]();
  double total=0, maxSuf=0, maxInt=0;
  for(int i=1;i<=ctx.numBoundaryTrig;i++)
    for(int k=0;k<ctx.numTimeStep;k++){
      ctx.timeSurfMeas[i][k]/=ctx.triNodes[ctx.boundaryTrigs[i]].Area;
      if(ctx.timeSurfMeas[i][k]>maxSuf) maxSuf=ctx.timeSurfMeas[i][k];
    }
  for(int i=1;i<=ctx.numElem;i++){
    if(ctx.medOptic[ctx.elems[i].MedIdx].mua==0.0)
      for(int k=0;k<ctx.numTimeStep;k++) ctx.timeAbsorption[i][k]=-1.0;
    else
      for(int k=0;k<ctx.numTimeStep;k++){
        if(ctx.medOptic[ctx.elems[i].MedIdx].g<1.0) CAbs[i]+=ctx.timeAbsorption[i][k];
        ctx.timeAbsorption[i][k]=(ctx.timeAbsorption[i][k]/ctx.elemNodes[i].Vol)
                                  /ctx.medOptic[ctx.elems[i].MedIdx].mua;
        if(ctx.timeAbsorption[i][k]>maxInt) maxInt=ctx.timeAbsorption[i][k];
      }
  }
  sort(CAbs,CAbs+ctx.numElem+1);
  for(int i=0;i<=ctx.numElem;i++) total+=CAbs[i];
  cerr<<"Absorbed Fraction: "<<total/double(ctx.totalPhoton)<<"\n";
  sufThreshhold=maxSuf/1e20; intThreshhold=maxInt/1e20;
  delete[] CAbs;
}

TiResult WriteResultASCII(const std::string& opt_f,const std::string& fem_f,
                      const std::string& src_f,const std::string& out_f,
                      SimContext& ctx, double sufTh, double intTh, int fmt){
  const std::string fb="/tmp/timos_tmp_result.dat";
  ofstream fout(out_f);
  if(!fout.good()){ cerr<<"Cannot write to "<<out_f<<"; using "<<fb<<"\n"; fout.open(fb); if(!fout.good()) return std::unexpected("cannot open output file"); }
  fout<<"% Optical filename:    "<<opt_f<<"\n"
      <<"% Fem mesh filename:   "<<fem_f<<"\n"
      <<"% Source filename:     "<<src_f<<"\n"
      <<"%      Num Photon:     "<<ctx.totalPhoton<<"\n"
      <<"% Num of threads:      "<<ctx.numThread<<"\n"
      <<"% Start random stream: "<<ctx.startRandIdx<<"\n";
  if(fmt==1||fmt==3){
    fout<<"1 "<<ctx.numBoundaryTrig<<" 1\n";
    for(int i=1;i<=ctx.numBoundaryTrig;i++){
      fout<<ctx.triNodes[ctx.boundaryTrigs[i]].N[0]<<" \t"
          <<ctx.triNodes[ctx.boundaryTrigs[i]].N[1]<<" \t"
          <<ctx.triNodes[ctx.boundaryTrigs[i]].N[2]<<" \t"
          <<ctx.triNodes[ctx.boundaryTrigs[i]].Area<<" \t"
          <<(ctx.surfMeas[i]<sufTh?0.0:ctx.surfMeas[i])<<"\n";
    }
  }
  if(fmt==2||fmt==3){
    fout<<"2 "<<ctx.numElem<<" 1\n";
    for(int i=1;i<=ctx.numElem;i++){
      fout<<ctx.elemNodes[i].N[0]<<" \t"<<ctx.elemNodes[i].N[1]<<" \t"
          <<ctx.elemNodes[i].N[2]<<" \t"<<ctx.elemNodes[i].N[3]<<" \t"
          <<ctx.elemNodes[i].Vol<<" \t"
          <<(ctx.absorption[i]<intTh?0.0:ctx.absorption[i])<<"\n";
    }
  }
  return {};
}

TiResult TimeWriteResultASCII(const std::string& opt_f,const std::string& fem_f,
                          const std::string& src_f,const std::string& out_f,
                          SimContext& ctx, double sufTh, double intTh, int fmt){
  const std::string fb="/tmp/timos_tmp_result.dat";
  ofstream fout(out_f);
  if(!fout.good()){ cerr<<"Cannot write to "<<out_f<<"; using "<<fb<<"\n"; fout.open(fb); if(!fout.good()) return std::unexpected("cannot open output file"); }
  fout<<"% Optical filename:    "<<opt_f<<"\n"
      <<"% Fem mesh filename:   "<<fem_f<<"\n"
      <<"% Source filename:     "<<src_f<<"\n"
      <<"%      Num Photon:     "<<ctx.totalPhoton<<"\n"
      <<"% Time domain. Step: "<<ctx.timeStep<<" ns, Steps: "<<ctx.numTimeStep<<"\n"
      <<"% Num of threads:      "<<ctx.numThread<<"\n"
      <<"% Start random stream: "<<ctx.startRandIdx<<"\n";
  if(fmt==1||fmt==3){
    fout<<"1 "<<ctx.numBoundaryTrig<<" "<<ctx.numTimeStep<<"\n";
    for(int i=1;i<=ctx.numBoundaryTrig;i++){
      fout<<ctx.triNodes[ctx.boundaryTrigs[i]].N[0]<<" \t"
          <<ctx.triNodes[ctx.boundaryTrigs[i]].N[1]<<" \t"
          <<ctx.triNodes[ctx.boundaryTrigs[i]].N[2]<<" \t"
          <<ctx.triNodes[ctx.boundaryTrigs[i]].Area<<" \t";
      for(int j=0;j<ctx.numTimeStep;j++) fout<<(ctx.timeSurfMeas[i][j]<sufTh?0.0:ctx.timeSurfMeas[i][j])<<" \t";
      fout<<"\n";
    }
  }
  if(fmt==2||fmt==3){
    fout<<"2 "<<ctx.numElem<<" "<<ctx.numTimeStep<<"\n";
    for(int i=1;i<=ctx.numElem;i++){
      fout<<ctx.elemNodes[i].N[0]<<" \t"<<ctx.elemNodes[i].N[1]<<" \t"
          <<ctx.elemNodes[i].N[2]<<" \t"<<ctx.elemNodes[i].N[3]<<" \t"
          <<ctx.elemNodes[i].Vol<<" \t";
      for(int j=0;j<ctx.numTimeStep;j++) fout<<(ctx.timeAbsorption[i][j]<intTh?0.0:ctx.timeAbsorption[i][j])<<" \t";
      fout<<"\n";
    }
  }
  return {};
}
