#include "globals.h"
#include "constants.h"
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
       << "  -m s     surface fluence\n"
       << "  -m i     internal fluence\n"
       << "  -m si    both\n"
       << "  -T dt N  time-domain: step dt (ns), N steps\n"
       << "  -t N     thread count (default 1)\n"
       << "  -r N     starting random-stream index (default 1)\n\n"
       << "All lengths in mm; optical coefficients in mm⁻¹.\n"
       << "Up to 1024 independent RNG streams; distribute across machines with -r.\n";
}

static bool check_input_file(const std::string& path){
  std::error_code ec;
  if(!std::filesystem::exists(path, ec) || ec){ cerr << "No such file: " << path << "\n"; return false; }
  if(!std::filesystem::is_regular_file(path, ec) || ec){ cerr << "Not a regular file: " << path << "\n"; return false; }
  return true;
}

static int flag_code(std::string_view s){
  if(s.size() >= 2 && s[0] == '-'){
    switch(s[1]){
    case 'p': return 1; case 'f': return 2; case 's': return 3;
    case 'm': return 4; case 'o': return 5; case 't': return 6;
    case 'r': return 7; case 'T': return 8;
    }
  }
  return -1;
}

static bool parse_argu(int argc, char* argv[],
                        string& opt_f, string& fem_f, string& src_f, string& out_f,
                        int& fmt, int& nThread, int& startRand,
                        bool& td, double& dt, double& invDt, double& invLsDt, int& nSteps){
  int i = 1;
  bool hasOpt=false, hasFem=false, hasSrc=false, hasOut=false;
  bool hasFmt=false, hasTh=false,  hasRd=false;

  while(i < argc){
    switch(flag_code(argv[i])){
    case 1: opt_f = argv[++i]; if(!check_input_file(opt_f)) return false; hasOpt=true; break;
    case 2: fem_f = argv[++i]; if(!check_input_file(fem_f)) return false; hasFem=true; break;
    case 3: src_f = argv[++i]; if(!check_input_file(src_f)) return false; hasSrc=true; break;
    case 4: {
      std::string_view m = argv[++i];
      if     (m=="s")           fmt=1;
      else if(m=="i")           fmt=2;
      else if(m=="si"||m=="is") fmt=3;
      else{ cerr<<"Unknown format '"<<m<<"'; using 's'.\n"; fmt=1; }
      hasFmt=true; break;
    }
    case 5: out_f = argv[++i]; hasOut=true; break;
    case 6: nThread = std::stoi(argv[++i]); if(nThread<=0||nThread>=1024) nThread=64; hasTh=true; break;
    case 7: startRand = std::stoi(argv[++i]); hasRd=true; break;
    case 8:
      td     = true;
      dt     = std::stod(argv[++i]);
      nSteps = std::stoi(argv[++i]);
      if(nSteps <= 1){ cerr<<"NumTimeStep must be > 1.\n"; return false; }
      invDt  = 1.0 / dt;
      invLsDt = INV_LIGHT_SPEED * invDt;
      break;
    default: print_usage(); return false;
    }
    i++;
  }

  if(!hasOpt||!hasFem||!hasSrc||!hasOut){ print_usage(); return false; }
  if(!hasTh)  nThread  = 1;
  if(!hasFmt) fmt      = 1;
  if(!hasRd)  startRand = 1;
  else if(startRand<=0 || startRand+nThread>1024){ cerr<<"StartRandIdx out of range.\n"; return false; }
  return true;
}

static double wall_time(){
  struct timeval tv;
  gettimeofday(&tv, nullptr);
  return tv.tv_sec + tv.tv_usec / 1e6;
}

int main(int argc, char* argv[]){
  cerr << "\n------------------------------------------------\n";
  if(argc < 9){ print_usage(); return -1; }

  string opt_f, fem_f, src_f, out_f;
  int fmt = 1;

  if(!parse_argu(argc, argv, opt_f, fem_f, src_f, out_f, fmt,
                 g_NumThread, g_StartRandIdx,
                 g_TimeDomain, g_TimeStep, g_InvTimeStep,
                 g_InvLightSpeedMutTimeStep, g_NumTimeStep))
    return -1;

  cerr << "Optical  filename: " << opt_f << "\n"
       << "Fem data filename: " << fem_f << "\n"
       << "Source   filename: " << src_f << "\n";

  cerr << "Begin read optical parameters.\n";
  if(ReadOpticalParameter(opt_f, g_SimpleOptic, g_NumMed, g_UniformBoundary, g_EnvRefIdx, g_MedOptic) < 0)
    return -1;
  cerr << "End read optical parameters.\n";

  cerr << "Begin read finite element file.\n";
  if(fem_read(fem_f, g_SimpleOptic, g_NumNode, g_NumElem, g_NumMed, g_Nodes, g_ElemNodes, g_Elems) < 0)
    return -1;
  cerr << "End read finite element file.\n";

  cerr << "Begin Preprocessing the finite element mesh.\n";
  if(PreProcessor(g_NumElem, g_NumNode, g_NumMed, g_NumBoundaryTrig, g_NumTrig,
                  g_Nodes, g_ElemNodes, g_Elems, g_TriNodes, g_Triangles, g_BoundaryTrigs) < 0)
    return -1;
  cerr << "End Preprocessing the finite element mesh.\n";

  cerr << "Begin read source file.\n";
  g_TotalPhoton = ReadSource(src_f, g_NumSource, g_Sources);
  cerr << "End read source file.\n";
  if(g_TotalPhoton <= 0) return -1;

  if(Prepare_Source(g_NumNode, g_NumElem, g_NumTrig, g_NumSource,
                    g_Sources, g_Elems, g_Nodes, g_ElemNodes, g_TriNodes, g_Triangles) < 0)
    return -1;

  cerr << "       Num Photon: " << g_TotalPhoton << "\n";
  if(g_TimeDomain)
    cerr << "Time domain: dt=" << g_TimeStep << " ns, steps=" << g_NumTimeStep << "\n";
  cerr << "Output filename:   " << out_f   << "\n"
       << "Output format:     " << fmt     << "\n"
       << "Threads:           " << g_NumThread   << "\n"
       << "Start RNG stream:  " << g_StartRandIdx << "\n";

  // Allocate result arrays
  if(!g_TimeDomain){
    g_SurfMeas   = new double[g_NumBoundaryTrig+1]();
    g_Absorption = new double[g_NumElem+1]();
  }else{
    g_TimeSurfMeas  .assign(g_NumBoundaryTrig+1, std::vector<double>(g_NumTimeStep, 0.0));
    g_TimeAbsorption.assign(g_NumElem+1,          std::vector<double>(g_NumTimeStep, 0.0));
  }

  g_NumIntersections = 0; g_NumSteps = 0; g_SimedPhoton = 0;

  // Launch simulation threads
  double t0 = wall_time();
  {
    std::vector<std::jthread> threads;
    threads.reserve(g_NumThread);
    for(int i = 0; i < g_NumThread; i++)
      threads.emplace_back([i]{ ThreadPhotonPropagation((void*)(uintptr_t)i); });
  }
  double t1 = wall_time();

  cerr << "Total simulation time: " << t1-t0 << " sec\n"
       << "Num of Intersection:   " << g_NumIntersections  << "\n"
       << "Num of Step:   "         << g_NumSteps  << "\n"
       << "Prepare output data:\n";

  double sufTh = 0, intTh = 0;
  if(g_TimeDomain)
    TimeAbsorptionToFluence(g_NumElem, g_NumBoundaryTrig, g_NumTimeStep, g_TotalPhoton,
                            sufTh, intTh, g_TriNodes, g_MedOptic, g_BoundaryTrigs,
                            g_Elems, g_ElemNodes, g_TimeSurfMeas, g_TimeAbsorption);
  else
    AbsorptionToFluence(g_NumElem, g_NumBoundaryTrig, g_TotalPhoton,
                        sufTh, intTh, g_TriNodes, g_MedOptic, g_BoundaryTrigs,
                        g_Elems, g_ElemNodes, g_SurfMeas, g_Absorption);

  cerr << "Write output file\n";
  if(g_TimeDomain)
    TimeWriteResultASCII(opt_f, fem_f, src_f, out_f, g_NumTimeStep, g_TimeStep,
                         g_TotalPhoton, g_NumThread, g_StartRandIdx, sufTh, intTh,
                         g_NumElem, g_NumBoundaryTrig, g_TriNodes, g_BoundaryTrigs,
                         g_Elems, g_ElemNodes, g_TimeSurfMeas, g_TimeAbsorption, fmt);
  else
    WriteResultASCII(opt_f, fem_f, src_f, out_f, g_TotalPhoton, g_NumThread, g_StartRandIdx,
                     sufTh, intTh, g_NumElem, g_NumBoundaryTrig, g_TriNodes, g_BoundaryTrigs,
                     g_Elems, g_ElemNodes, g_SurfMeas, g_Absorption, fmt);

  cerr << "Done\n------------------------------------------------\n";

  delete[] g_Sources; delete[] g_MedOptic; delete[] g_BoundaryTrigs;
  delete[] g_ElemNodes; delete[] g_Elems; delete[] g_Nodes;
  delete[] g_TriNodes; delete[] g_Triangles;
  if(!g_TimeDomain){ delete[] g_SurfMeas; delete[] g_Absorption; }
  else{ g_TimeSurfMeas.clear(); g_TimeAbsorption.clear(); }

  return 0;
}
