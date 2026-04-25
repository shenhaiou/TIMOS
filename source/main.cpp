#include "simcontext.h"
#include "mesh.h"
#include "optics.h"
#include "source_io.h"
#include "simulation.h"
#include "output.h"

#include <iostream>
#include <string>
#include <string_view>
#include <filesystem>
#include <thread>
#include <vector>
#include <sys/time.h>

using namespace std;

static void print_usage(){
  cerr << "\nTIM-OS version 3.0\n"
       << "Usage: timos -p opt_file -f mesh_file -s src_file -m [s|i|si] "
          "-o out_file [-T dt steps] [-t threads] [-r start_rand]\n\n"
       << "  -m s|i|si   output format: surface / internal / both\n"
       << "  -T dt N     time-domain simulation: step dt (ns), N steps\n"
       << "  -t N        thread count (default 1)\n"
       << "  -r N        starting RNG stream index (default 1)\n\n"
       << "All lengths in mm; optical coefficients in mm⁻¹.\n";
}

static bool check_file(const std::string& p){
  std::error_code ec;
  if(!std::filesystem::exists(p,ec)||ec){ cerr<<"No such file: "<<p<<"\n"; return false; }
  if(!std::filesystem::is_regular_file(p,ec)||ec){ cerr<<"Not a regular file: "<<p<<"\n"; return false; }
  return true;
}

static int flag_code(std::string_view s){
  if(s.size()>=2&&s[0]=='-'){
    switch(s[1]){ case 'p':return 1;case 'f':return 2;case 's':return 3;
                  case 'm':return 4;case 'o':return 5;case 't':return 6;
                  case 'r':return 7;case 'T':return 8; }
  }
  return -1;
}

static bool parse_argu(int argc, char* argv[],
                        string& opt_f,string& fem_f,string& src_f,string& out_f,
                        int& fmt, SimContext& ctx){
  int i=1;
  bool hasOpt=false,hasFem=false,hasSrc=false,hasOut=false;
  bool hasFmt=false,hasTh=false,hasRd=false;
  while(i<argc){
    switch(flag_code(argv[i])){
    case 1: opt_f=argv[++i]; if(!check_file(opt_f))return false; hasOpt=true; break;
    case 2: fem_f=argv[++i]; if(!check_file(fem_f))return false; hasFem=true; break;
    case 3: src_f=argv[++i]; if(!check_file(src_f))return false; hasSrc=true; break;
    case 4:{
      std::string_view m=argv[++i];
      if(m=="s")fmt=1; else if(m=="i")fmt=2; else if(m=="si"||m=="is")fmt=3;
      else if(m=="t")fmt=4; else{ cerr<<"Unknown format '"<<m<<"'; using 's'.\n";fmt=1; }
      hasFmt=true; break;
    }
    case 5: out_f=argv[++i]; hasOut=true; break;
    case 6: ctx.numThread=std::stoi(argv[++i]); if(ctx.numThread<=0||ctx.numThread>=1024)ctx.numThread=64; hasTh=true; break;
    case 7: ctx.startRandIdx=std::stoi(argv[++i]); hasRd=true; break;
    case 8:
      ctx.timeDomain=true;
      ctx.timeStep=std::stod(argv[++i]);
      ctx.numTimeStep=std::stoi(argv[++i]);
      if(ctx.numTimeStep<=1){ cerr<<"numTimeStep must be > 1.\n"; return false; }
      ctx.invTimeStep=1.0/ctx.timeStep;
      ctx.invLightSpeedMutTimeStep=INV_LIGHT_SPEED*ctx.invTimeStep;
      break;
    default: print_usage(); return false;
    }
    i++;
  }
  if(!hasOpt||!hasFem||!hasSrc||!hasOut){ print_usage(); return false; }
  if(!hasTh)  ctx.numThread=1;
  if(!hasFmt) fmt=1;
  if(!hasRd)  ctx.startRandIdx=1;
  else if(ctx.startRandIdx<=0||ctx.startRandIdx+ctx.numThread>1024){
    cerr<<"startRandIdx out of range.\n"; return false;
  }
  return true;
}

static double wall_time(){
  struct timeval tv; gettimeofday(&tv,nullptr);
  return tv.tv_sec+tv.tv_usec/1e6;
}

int main(int argc, char* argv[]){
  cerr<<"\n------------------------------------------------\n";
  if(argc<9){ print_usage(); return -1; }

  SimContext ctx;
  string opt_f,fem_f,src_f,out_f;
  int fmt=1;

  if(!parse_argu(argc,argv,opt_f,fem_f,src_f,out_f,fmt,ctx)) return -1;

  cerr<<"Optical  filename: "<<opt_f<<"\n"
      <<"Fem data filename: "<<fem_f<<"\n"
      <<"Source   filename: "<<src_f<<"\n";

  // Lambda: report error and return -1 if TiResult holds an error
  auto check = [](TiResult r, const char* step) -> bool {
    if(!r){ cerr << "Error in " << step << ": " << r.error() << "\n"; return false; }
    return true;
  };

  cerr<<"Begin read optical parameters.\n";
  if(!check(ReadOpticalParameter(opt_f,ctx),"ReadOpticalParameter")) return -1;
  cerr<<"End read optical parameters.\n";

  cerr<<"Begin read finite element file.\n";
  if(!check(fem_read(fem_f,ctx),"fem_read")) return -1;
  cerr<<"End read finite element file.\n";

  cerr<<"Begin Preprocessing the finite element mesh.\n";
  if(!check(PreProcessor(ctx),"PreProcessor")) return -1;
  cerr<<"End Preprocessing the finite element mesh.\n";

  cerr<<"Begin read source file.\n";
  {
    auto src_result = ReadSource(src_f, ctx);
    if(!src_result){ cerr<<"Error reading source: "<<src_result.error()<<"\n"; return -1; }
    ctx.totalPhoton = *src_result;
  }
  cerr<<"End read source file.\n";
  if(ctx.totalPhoton<=0){ cerr<<"No photons in source file.\n"; return -1; }

  if(!check(Prepare_Source(ctx),"Prepare_Source")) return -1;

  cerr<<"       Num Photon: "<<ctx.totalPhoton<<"\n";
  if(ctx.timeDomain)
    cerr<<"Time domain: dt="<<ctx.timeStep<<" ns, steps="<<ctx.numTimeStep<<"\n";
  cerr<<"Output: "<<out_f<<"  format="<<fmt
      <<"  threads="<<ctx.numThread<<"  startRand="<<ctx.startRandIdx<<"\n";

  // Allocate result arrays
  if(!ctx.timeDomain){
    ctx.surfMeas  .assign(ctx.numBoundaryTrig+1, 0.0);
    ctx.absorption.assign(ctx.numElem+1,         0.0);
  }else{
    ctx.timeSurfMeas  .assign(ctx.numBoundaryTrig+1,vector<double>(ctx.numTimeStep,0.0));
    ctx.timeAbsorption.assign(ctx.numElem+1,         vector<double>(ctx.numTimeStep,0.0));
  }

  ctx.numIntersections=0; ctx.numSteps=0; ctx.simedPhoton=0;

  // Launch threads via std::jthread — auto-joins on scope exit
  double t0=wall_time();
  {
    std::vector<ThreadArg> args(ctx.numThread);
    std::vector<std::jthread> threads;
    threads.reserve(ctx.numThread);
    for(int i=0;i<ctx.numThread;i++){
      args[i]={i,&ctx};
      threads.emplace_back([&args,i]{ ThreadPhotonPropagation(&args[i]); });
    }
  }
  double t1=wall_time();

  cerr<<"Total simulation time: "<<t1-t0<<" sec\n"
      <<"Num of Intersection:   "<<ctx.numIntersections<<"\n"
      <<"Num of Step:   "<<ctx.numSteps<<"\n"
      <<"Prepare output data:\n";

  double sufTh=0,intTh=0;
  if(ctx.timeDomain) TimeAbsorptionToFluence(ctx,sufTh,intTh);
  else               AbsorptionToFluence    (ctx,sufTh,intTh);

  cerr<<"Write output file\n";
  TiResult wr = ctx.timeDomain
    ? TimeWriteResultASCII(opt_f,fem_f,src_f,out_f,ctx,sufTh,intTh,fmt)
    :      WriteResultASCII(opt_f,fem_f,src_f,out_f,ctx,sufTh,intTh,fmt);
  if(!wr){ cerr<<"Write failed: "<<wr.error()<<"\n"; return -1; }

  cerr<<"Done\n------------------------------------------------\n";
  // ctx goes out of scope here; all std::vector members free automatically.
  return 0;
}
