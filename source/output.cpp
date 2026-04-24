#include "output.h"
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

void AbsorptionToFluence(int numElem, int numBoundaryTrig,
                         long long int TotalPhoton,
                         double sufThreshhold, double intThreshhold,
                         TTriNode* triNodes, TMedOptic* medOptic,
                         int* boundaryTrigs, TElem* Elems, TElemNode* ElemNodes,
                         double* SurfMeas, double* Absorption){
  double* CAbs = new double[numElem+1]();
  double TotalAbs = 0, maxSuf = 0, maxInt = 0;

  for(int i = 1; i <= numBoundaryTrig; i++){
    SurfMeas[i] /= triNodes[boundaryTrigs[i]].Area;
    if(SurfMeas[i] > maxSuf) maxSuf = SurfMeas[i];
  }
  for(int i = 1; i <= numElem; i++){
    if(medOptic[Elems[i].MedIdx].mua == 0.0){
      Absorption[i] = -1.0;
    }else{
#ifdef NOPUREGLASS
      CAbs[i] = Absorption[i];
#else
      if(medOptic[Elems[i].MedIdx].g < 1.0) CAbs[i] = Absorption[i];
#endif
      Absorption[i] = (Absorption[i] / ElemNodes[i].Vol) / medOptic[Elems[i].MedIdx].mua;
      if(Absorption[i] > maxInt) maxInt = Absorption[i];
    }
  }
  sort(CAbs, CAbs + numElem + 1);
  for(int i = 0; i <= numElem; i++) TotalAbs += CAbs[i];
  cerr << "Absorbed Fraction: " << TotalAbs / double(TotalPhoton) << endl;
  sufThreshhold = maxSuf / 1e20;
  intThreshhold = maxInt / 1e20;
  delete[] CAbs;
}

void TimeAbsorptionToFluence(int numElem, int numBoundaryTrig, int NumTimeStep,
                             long long int TotalPhoton,
                             double sufThreshhold, double intThreshhold,
                             TTriNode* triNodes, TMedOptic* medOptic,
                             int* boundaryTrigs, TElem* Elems, TElemNode* ElemNodes,
                             std::vector<std::vector<double>>& TimeSurfMeas,
                             std::vector<std::vector<double>>& TimeAbsorption){
  double* CAbs = new double[numElem+1]();
  double TotalAbs = 0, maxSuf = 0, maxInt = 0;

  for(int i = 1; i <= numBoundaryTrig; i++)
    for(int k = 0; k < NumTimeStep; k++){
      TimeSurfMeas[i][k] /= triNodes[boundaryTrigs[i]].Area;
      if(TimeSurfMeas[i][k] > maxSuf) maxSuf = TimeSurfMeas[i][k];
    }
  for(int i = 1; i <= numElem; i++){
    if(medOptic[Elems[i].MedIdx].mua == 0.0){
      for(int k = 0; k < NumTimeStep; k++) TimeAbsorption[i][k] = -1.0;
    }else{
      for(int k = 0; k < NumTimeStep; k++){
        if(medOptic[Elems[i].MedIdx].g < 1.0) CAbs[i] += TimeAbsorption[i][k];
        TimeAbsorption[i][k] = (TimeAbsorption[i][k] / ElemNodes[i].Vol)
                                / medOptic[Elems[i].MedIdx].mua;
        if(TimeAbsorption[i][k] > maxInt) maxInt = TimeAbsorption[i][k];
      }
    }
  }
  sort(CAbs, CAbs + numElem + 1);
  for(int i = 0; i <= numElem; i++) TotalAbs += CAbs[i];
  cerr << "Absorbed Fraction: " << TotalAbs / double(TotalPhoton) << endl;
  sufThreshhold = maxSuf / 1e20;
  intThreshhold = maxInt / 1e20;
  delete[] CAbs;
}

bool WriteResultASCII(const std::string& optical_filename,
                      const std::string& fem_filename,
                      const std::string& source_filename,
                      const std::string& output_filename,
                      long long int TotalPhoton, int NumThread, int StartRandIdx,
                      double sufThreshhold, double intThreshhold,
                      int numElem, int numBoundaryTrig,
                      TTriNode* triNodes, int* boundaryTrigs,
                      TElem* Elems, TElemNode* ElemNodes,
                      double* SurfMeas, double* Absorption,
                      int output_format){
  const std::string fallback = "/tmp/timos_tmp_result.dat";
  ofstream fout(output_filename);
  if(!fout.good()){
    cerr << "Could not write to " << output_filename
         << ". Falling back to " << fallback << endl;
    fout.open(fallback);
    if(!fout.good()) return false;
  }

  fout << "% Optical filename:    " << optical_filename << "\n"
       << "% Fem mesh filename:   " << fem_filename     << "\n"
       << "% Source filename:     " << source_filename  << "\n"
       << "%      Num Photon:     " << TotalPhoton      << "\n"
       << "% Num of threads:      " << NumThread        << "\n"
       << "% Start random stream: " << StartRandIdx     << "\n";

  int NumTimeStep = 1;
  if(output_format == 1 || output_format == 3){
    fout << "1 " << numBoundaryTrig << " " << NumTimeStep << "\n";
    for(int i = 1; i <= numBoundaryTrig; i++){
      fout << triNodes[boundaryTrigs[i]].N[0] << " \t"
           << triNodes[boundaryTrigs[i]].N[1] << " \t"
           << triNodes[boundaryTrigs[i]].N[2] << " \t"
           << triNodes[boundaryTrigs[i]].Area << " \t";
      fout << ((SurfMeas[i] < sufThreshhold) ? 0.0 : SurfMeas[i]) << "\n";
    }
  }
  if(output_format == 2 || output_format == 3){
    fout << "2 " << numElem << " " << NumTimeStep << "\n";
    for(int i = 1; i <= numElem; i++){
      fout << ElemNodes[i].N[0] << " \t" << ElemNodes[i].N[1] << " \t"
           << ElemNodes[i].N[2] << " \t" << ElemNodes[i].N[3] << " \t"
           << ElemNodes[i].Vol  << " \t";
      fout << ((Absorption[i] < intThreshhold) ? 0.0 : Absorption[i]) << "\n";
    }
  }
  return true;
}

bool TimeWriteResultASCII(const std::string& optical_filename,
                          const std::string& fem_filename,
                          const std::string& source_filename,
                          const std::string& output_filename,
                          int NumTimeStep, double TimeStep,
                          long long int TotalPhoton, int NumThread, int StartRandIdx,
                          double sufThreshhold, double intThreshhold,
                          int numElem, int numBoundaryTrig,
                          TTriNode* triNodes, int* boundaryTrigs,
                          TElem* Elems, TElemNode* ElemNodes,
                          std::vector<std::vector<double>>& TimeSurfMeas,
                          std::vector<std::vector<double>>& TimeAbsorption,
                          int output_format){
  const std::string fallback = "/tmp/timos_tmp_result.dat";
  ofstream fout(output_filename);
  if(!fout.good()){
    cerr << "Could not write to " << output_filename
         << ". Falling back to " << fallback << endl;
    fout.open(fallback);
    if(!fout.good()) return false;
  }

  fout << "% Optical filename:    " << optical_filename << "\n"
       << "% Fem mesh filename:   " << fem_filename     << "\n"
       << "% Source filename:     " << source_filename  << "\n"
       << "%      Num Photon:     " << TotalPhoton      << "\n"
       << "% Time domain setting. Step size: " << TimeStep
       << " ns, Num Step: " << NumTimeStep             << "\n"
       << "% Num of threads:      " << NumThread        << "\n"
       << "% Start random stream: " << StartRandIdx     << "\n";

  if(output_format == 1 || output_format == 3){
    fout << "1 " << numBoundaryTrig << " " << NumTimeStep << "\n";
    for(int i = 1; i <= numBoundaryTrig; i++){
      fout << triNodes[boundaryTrigs[i]].N[0] << " \t"
           << triNodes[boundaryTrigs[i]].N[1] << " \t"
           << triNodes[boundaryTrigs[i]].N[2] << " \t"
           << triNodes[boundaryTrigs[i]].Area << " \t";
      for(int j = 0; j < NumTimeStep; j++)
        fout << ((TimeSurfMeas[i][j] < sufThreshhold) ? 0.0 : TimeSurfMeas[i][j]) << " \t";
      fout << "\n";
    }
  }
  if(output_format == 2 || output_format == 3){
    fout << "2 " << numElem << " " << NumTimeStep << "\n";
    for(int i = 1; i <= numElem; i++){
      fout << ElemNodes[i].N[0] << " \t" << ElemNodes[i].N[1] << " \t"
           << ElemNodes[i].N[2] << " \t" << ElemNodes[i].N[3] << " \t"
           << ElemNodes[i].Vol  << " \t";
      for(int j = 0; j < NumTimeStep; j++)
        fout << ((TimeAbsorption[i][j] < intThreshhold) ? 0.0 : TimeAbsorption[i][j]) << " \t";
      fout << "\n";
    }
  }
  return true;
}
