#include "output.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <format>
#include <vector>
#include <charconv>

using namespace std;

void AbsorptionToFluence(SimContext& ctx,
                         double& sufThreshhold, double& intThreshhold){
  vector<double> CAbs(ctx.numElem + 1, 0.0);
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
  sort(CAbs.begin(), CAbs.end());
  for(double v : CAbs) total += v;
  cout<<"Absorbed Fraction: "<<total/double(ctx.totalPhoton)<<"\n";
  sufThreshhold=maxSuf/1e20; intThreshhold=maxInt/1e20;
}

void TimeAbsorptionToFluence(SimContext& ctx,
                             double& sufThreshhold, double& intThreshhold){
  vector<double> CAbs(ctx.numElem + 1, 0.0);
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
  sort(CAbs.begin(), CAbs.end());
  for(double v : CAbs) total += v;
  cout<<"Absorbed Fraction: "<<total/double(ctx.totalPhoton)<<"\n";
  sufThreshhold=maxSuf/1e20; intThreshhold=maxInt/1e20;
}

static inline void fast_append_double(std::string& s, double val){
  char buf[64];
  auto [ptr, ec] = std::to_chars(buf, buf + sizeof(buf), val, std::chars_format::general, 6);
  if (ec == std::errc()) {
    s.append(buf, ptr - buf);
  } else {
    s += "0";
  }
  s += ' ';
}

TiResult WriteResultASCII(const std::string& opt_f,const std::string& fem_f,
                      const std::string& src_f,const std::string& out_f,
                      SimContext& ctx, double sufTh, double intTh, int fmt){
  ofstream fout(out_f);
  if(!fout.good()) return std::unexpected("cannot open output file");

  std::vector<char> buffer(1024 * 1024);
  fout.rdbuf()->pubsetbuf(buffer.data(), buffer.size());

  fout << std::format("% Optical filename:    {}\n", opt_f)
       << std::format("% Fem mesh filename:   {}\n", fem_f)
       << std::format("% Source filename:     {}\n", src_f)
       << std::format("%      Num Photon:     {}\n", ctx.totalPhoton)
       << std::format("% Num of threads:      {}\n", ctx.numThread)
       << std::format("% Start random stream: {}\n", ctx.startRandIdx);

  std::string line;
  line.reserve(1024);

  if(fmt==1||fmt==3){
    fout << std::format("1 {} 1\n", ctx.numBoundaryTrig);
    for(int i=1;i<=ctx.numBoundaryTrig;i++){
      const auto& t = ctx.triNodes[ctx.boundaryTrigs[i]];
      double val = (ctx.surfMeas[i] < sufTh) ? 0.0 : ctx.surfMeas[i];
      line.clear();
      std::format_to(std::back_inserter(line), "{} {} {} {} ", t.N[0], t.N[1], t.N[2], t.Area);
      fast_append_double(line, val);
      line += '\n';
      fout << line;
    }
  }
  if(fmt==2||fmt==3){
    fout << std::format("2 {} 1\n", ctx.numElem);
    for(int i=1;i<=ctx.numElem;i++){
      const auto& e = ctx.elemNodes[i];
      double val = (ctx.absorption[i] < intTh) ? 0.0 : ctx.absorption[i];
      line.clear();
      std::format_to(std::back_inserter(line), "{} {} {} {} {} ", e.N[0], e.N[1], e.N[2], e.N[3], e.Vol);
      fast_append_double(line, val);
      line += '\n';
      fout << line;
    }
  }
  return {};
}

TiResult TimeWriteResultASCII(const std::string& opt_f,const std::string& fem_f,
                          const std::string& src_f,const std::string& out_f,
                          SimContext& ctx, double sufTh, double intTh, int fmt){
  ofstream fout(out_f);
  if(!fout.good()) return std::unexpected("cannot open output file");

  std::vector<char> buffer(1024 * 1024);
  fout.rdbuf()->pubsetbuf(buffer.data(), buffer.size());

  fout << std::format("% Optical filename:    {}\n", opt_f)
       << std::format("% Fem mesh filename:   {}\n", fem_f)
       << std::format("% Source filename:     {}\n", src_f)
       << std::format("%      Num Photon:     {}\n", ctx.totalPhoton)
       << std::format("% Time domain. Step: {} ns, Steps: {}\n", ctx.timeStep, ctx.numTimeStep)
       << std::format("% Num of threads:      {}\n", ctx.numThread)
       << std::format("% Start random stream: {}\n", ctx.startRandIdx);

  std::string line;
  line.reserve(512 + 24 * ctx.numTimeStep);

  if(fmt==1||fmt==3){
    fout << std::format("1 {} {}\n", ctx.numBoundaryTrig, ctx.numTimeStep);
    for(int i=1;i<=ctx.numBoundaryTrig;i++){
      const auto& t = ctx.triNodes[ctx.boundaryTrigs[i]];
      line.clear();
      std::format_to(std::back_inserter(line), "{} {} {} {} ", t.N[0], t.N[1], t.N[2], t.Area);
      for(int j=0;j<ctx.numTimeStep;j++){
        double val = (ctx.timeSurfMeas[i][j] < sufTh) ? 0.0 : ctx.timeSurfMeas[i][j];
        fast_append_double(line, val);
      }
      line += '\n';
      fout << line;
    }
  }
  if(fmt==2||fmt==3){
    fout << std::format("2 {} {}\n", ctx.numElem, ctx.numTimeStep);
    for(int i=1;i<=ctx.numElem;i++){
      const auto& e = ctx.elemNodes[i];
      line.clear();
      std::format_to(std::back_inserter(line), "{} {} {} {} {} ", e.N[0], e.N[1], e.N[2], e.N[3], e.Vol);
      for(int j=0;j<ctx.numTimeStep;j++){
        double val = (ctx.timeAbsorption[i][j] < intTh) ? 0.0 : ctx.timeAbsorption[i][j];
        fast_append_double(line, val);
      }
      line += '\n';
      fout << line;
    }
  }
  return {};
}

TiResult WriteGridASCII(const std::string& out_f, SimContext& ctx){
  std::string grid_out = out_f + ".grid";
  ofstream fout(grid_out);
  if(!fout.good()) return std::unexpected("cannot open grid output file");

  int nt = ctx.timeDomain ? ctx.numTimeStep : 1;
  fout << "% Cylindrical Grid: Nr=" << ctx.gridNr << " Ny=" << ctx.gridNy << " Nt=" << nt << "\n";
  fout << "% Rmax=" << ctx.gridRMax << " Ymax=" << ctx.gridYMax << "\n";
  fout << "3 " << ctx.gridNr << " " << ctx.gridNy << " " << nt << "\n";

  double dr = ctx.gridRMax / ctx.gridNr;
  double dy = ctx.gridYMax / ctx.gridNy;
  double PI = 3.14159265358979323846;

  std::string line;
  line.reserve(1024);
  for(int r=0; r<ctx.gridNr; r++){
    // Calculate volume of the cylindrical ring at radius r
    double r_in = r * dr;
    double r_out = (r + 1) * dr;
    double ring_vol = PI * (r_out * r_out - r_in * r_in) * dy;
    
    for(int y=0; y<ctx.gridNy; y++){
      line.clear();
      std::format_to(std::back_inserter(line), "{} {} ", r, y);
      for(int t=0; t<nt; t++){
        // Cylindrical Grid data already stores Fluence Weight (E / mua).
        // Complete the fluence calculation by dividing by ring volume.
        double val = ctx.cylindricalGrid[r][y][t] / ring_vol;
        fast_append_double(line, val);
      }
      line += '\n';
      fout << line;
    }
  }
  return {};
}


