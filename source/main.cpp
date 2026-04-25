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
#include <getopt.h>    // POSIX getopt_long
#include <sys/time.h>

using namespace std;

// --------------------------------------------------------------------- usage
static void print_usage(){
  cout <<
    "\nTIM-OS v3.3 — Tetrahedral-mesh Inhomogeneous Monte-Carlo Optical Simulator\n\n"
    "Usage:\n"
    "  timos -p OPT -f MESH -s SRC -m MODE [options]\n\n"
    "Required:\n"
    "  -p, --optical FILE    optical parameter file (.opt)\n"
    "  -f, --mesh    FILE    finite-element mesh file (.mesh)\n"
    "  -s, --source  FILE    light source definition (.source)\n"
    "  -m, --mode    MODE    output mode: s=surface, i=internal, si=both\n\n"
    "Optional:\n"
    "  -o, --output  FILE    output filename (.dat) [omit for benchmark mode]\n"
    "  -t, --threads N       number of threads (default: 1)\n"
    "  -r, --rand-start N    starting RNG stream index for distributed runs\n"
    "                        (default: 1; max streams = 1024)\n"
    "  -T, --time-step DT    enable time-domain mode; DT = step size in ns\n"
    "  -N, --num-steps N     number of time steps (required when -T is set)\n"
    "  -h, --help            show this message\n\n"
    "Notes:\n"
    "  All lengths in mm; optical coefficients in mm^-1.\n"
    "  Omitting -o enables Benchmark Mode (skips writing results to disk).\n"
    "  For distributed runs across machines use non-overlapping -r values,\n"
    "  e.g. machine 1: -r 1 -t 16 ; machine 2: -r 17 -t 16.\n\n";
}

// ---------------------------------------------------------------- file check
static bool check_file(const std::string& p){
  std::error_code ec;
  if(!std::filesystem::exists(p,ec)||ec){ cout<<"No such file: "<<p<<"\n"; return false; }
  if(!std::filesystem::is_regular_file(p,ec)||ec){ cout<<"Not a regular file: "<<p<<"\n"; return false; }
  return true;
}

// ------------------------------------------------------------ argument parser
// Uses POSIX getopt_long so both short flags (-p) and long options (--optical)
// are accepted. The non-standard -T flag takes a second positional argument
// (num_steps) which is consumed from optind after getopt sees the DT value.
static bool parse_argu(int argc, char* argv[],
                        string& opt_f, string& fem_f,
                        string& src_f, string& out_f,
                        int& fmt, SimContext& ctx){
  // clang-format off
  static const option long_opts[] = {
    {"optical",    required_argument, nullptr, 'p'},
    {"mesh",       required_argument, nullptr, 'f'},
    {"source",     required_argument, nullptr, 's'},
    {"mode",       required_argument, nullptr, 'm'},
    {"output",     required_argument, nullptr, 'o'},
    {"threads",    required_argument, nullptr, 't'},
    {"rand-start", required_argument, nullptr, 'r'},
    {"time-step",  required_argument, nullptr, 'T'},
    {"num-steps",  required_argument, nullptr, 'N'},
    {"help",       no_argument,       nullptr, 'h'},
    {nullptr,      0,                 nullptr,  0 }
  };
  // clang-format on
  static const char* short_opts = "p:f:s:m:o:t:r:T:N:h";

  bool hasOpt=false, hasFem=false, hasSrc=false, hasOut=false;
  bool hasFmt=false, hasTh=false,  hasRd=false;

  int c;
  while((c = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1){
    switch(c){
    case 'p':
      opt_f = optarg; if(!check_file(opt_f)) return false; hasOpt=true; break;
    case 'f':
      fem_f = optarg; if(!check_file(fem_f)) return false; hasFem=true; break;
    case 's':
      src_f = optarg; if(!check_file(src_f)) return false; hasSrc=true; break;
    case 'm': {
      std::string_view m = optarg;
      if     (m=="s")            fmt=1;
      else if(m=="i")            fmt=2;
      else if(m=="si"||m=="is")  fmt=3;
      else{ cout<<"Unknown mode '"<<m<<"'; using 's'.\n"; fmt=1; }
      hasFmt=true; break;
    }
    case 'o':
      out_f = optarg; hasOut=true; break;
    case 't':
      try {
        ctx.numThread = std::stoi(optarg);
      } catch (...) { cout << "Invalid thread count: " << optarg << "\n"; return false; }
      if(ctx.numThread<=0||ctx.numThread>=1024) ctx.numThread=64;
      hasTh=true; break;
    case 'r':
      try {
        ctx.startRandIdx = std::stoi(optarg);
      } catch (...) { cout << "Invalid random start index: " << optarg << "\n"; return false; }
      hasRd=true; break;
    case 'T':
      // -T DT : enables time-domain, sets step size; -N must also be provided.
      ctx.timeDomain = true;
      try {
        ctx.timeStep = std::stod(optarg);
      } catch (...) { cout << "Invalid time step: " << optarg << "\n"; return false; }
      ctx.invTimeStep = 1.0 / ctx.timeStep;
      ctx.invLightSpeedMutTimeStep = INV_LIGHT_SPEED * ctx.invTimeStep;
      // Legacy form: if next arg is a bare number (not a flag), treat it as
      // num_steps to preserve backwards-compatibility with "-T dt N" usage.
      if(optind < argc && argv[optind][0] != '-'){
        try {
          ctx.numTimeStep = std::stoi(argv[optind++]);
        } catch (...) { optind--; } // backtrack if next arg wasn't a valid int
      }
      break;
    case 'N':
      try {
        ctx.numTimeStep = std::stoi(optarg);
      } catch (...) { cout << "Invalid number of steps: " << optarg << "\n"; return false; }
      break;
    case 'h':
      print_usage(); return false;
    default:
      print_usage(); return false;
    }
  }

  if(!hasOpt||!hasFem||!hasSrc){ print_usage(); return false; }
  if(!hasTh)  ctx.numThread    = 1;
  if(!hasFmt) fmt              = 1;
  if(!hasRd)  ctx.startRandIdx = 1;
  else if(ctx.startRandIdx<=0 || ctx.startRandIdx+ctx.numThread>1024){
    cout<<"startRandIdx out of range.\n"; return false;
  }
  if(ctx.timeDomain && ctx.numTimeStep<=1){
    cout<<"--num-steps N must be > 1 for time-domain mode.\n"; return false;
  }
  return true;
}

static double wall_time(){
  struct timeval tv; gettimeofday(&tv,nullptr);
  return tv.tv_sec+tv.tv_usec/1e6;
}

int main(int argc, char* argv[]){
  cout<<"\n------------------------------------------------\n";
  if(argc<2){ print_usage(); return -1; }

  SimContext ctx;
  string opt_f,fem_f,src_f,out_f;
  int fmt=1;

  if(!parse_argu(argc,argv,opt_f,fem_f,src_f,out_f,fmt,ctx)) return -1;

  cout<<"Optical  filename: "<<opt_f<<"\n"
      <<"Fem data filename: "<<fem_f<<"\n"
      <<"Source   filename: "<<src_f<<"\n";

  // Lambda: report error and return -1 if TiResult holds an error
  auto check = [](TiResult r, const char* step) -> bool {
    if(!r){ cout << "Error in " << step << ": " << r.error() << "\n"; return false; }
    return true;
  };

  cout<<"Begin read optical parameters.\n";
  if(!check(ReadOpticalParameter(opt_f,ctx),"ReadOpticalParameter")) return -1;
  cout<<"End read optical parameters.\n";

  double mt0 = wall_time();
  cout<<"Begin read finite element file.\n";
  if(!check(fem_read(fem_f,ctx),"fem_read")) return -1;
  cout<<"End read finite element file (" << wall_time()-mt0 << " sec).\n";

  double pt0 = wall_time();
  cout<<"Begin Preprocessing the finite element mesh.\n";
  if(!check(PreProcessor(ctx),"PreProcessor")) return -1;
  cout<<"End Preprocessing the finite element mesh (" << wall_time()-pt0 << " sec).\n";

  cout<<"Begin read source file.\n";
  {
    auto src_result = ReadSource(src_f, ctx);
    if(!src_result){ cout<<"Error reading source: "<<src_result.error()<<"\n"; return -1; }
    ctx.totalPhoton = *src_result;
  }
  cout<<"End read source file.\n";
  if(ctx.totalPhoton<=0){ cout<<"No photons in source file.\n"; return -1; }

  if(!check(Prepare_Source(ctx),"Prepare_Source")) return -1;

  cout<<"       Num Photon: "<<ctx.totalPhoton<<"\n";
  if(ctx.timeDomain)
    cout<<"Time domain: dt="<<ctx.timeStep<<" ns, steps="<<ctx.numTimeStep<<"\n";
  cout<<"Output: "<<out_f<<"  format="<<fmt
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
  double st0=wall_time();
  {
    std::vector<ThreadArg> args(ctx.numThread);
    std::vector<std::jthread> threads;
    threads.reserve(ctx.numThread);
    for(int i=0;i<ctx.numThread;i++){
      args[i]={i,&ctx};
      threads.emplace_back([&args,i]{ ThreadPhotonPropagation(&args[i]); });
    }
  }
  double st1=wall_time();

  cout<<"Total simulation time: "<<st1-st0<<" sec\n"
      <<"Num of Intersection:   "<<ctx.numIntersections<<"\n"
      <<"Num of Step:   "<<ctx.numSteps<<"\n"
      <<"Prepare output data:\n";

  double sufTh=0,intTh=0;
  if(ctx.timeDomain) TimeAbsorptionToFluence(ctx,sufTh,intTh);
  else               AbsorptionToFluence    (ctx,sufTh,intTh);

  if(!out_f.empty()){
    double wt0 = wall_time();
    cout<<"Write output file\n";
    TiResult wr = ctx.timeDomain
      ? TimeWriteResultASCII(opt_f,fem_f,src_f,out_f,ctx,sufTh,intTh,fmt)
      :      WriteResultASCII(opt_f,fem_f,src_f,out_f,ctx,sufTh,intTh,fmt);
    if(!wr){ cout<<"Write failed: "<<wr.error()<<"\n"; return -1; }
    cout<<"End write output file (" << wall_time()-wt0 << " sec).\n";
  }else{
    cout<<"Benchmark mode: Skipping output file write.\n";
  }

  cout<<"Done\n------------------------------------------------\n";
  // ctx goes out of scope here; all std::vector members free automatically.
  return 0;
}
