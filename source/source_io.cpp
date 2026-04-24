#include "source_io.h"
#include "io_utils.h"
#include "constants.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

long long int ReadSource(const std::string& filename, SimContext& ctx){
  long long int total = 0;
  ifstream fin(filename);
  if(!fin.good()){ cerr<<"\tCould not open file: "<<filename<<"\n"; return -1; }

  int curLine = 0;
  skip_comments(fin); fin >> ctx.numSource; curLine++;
  if(fin.fail()) cerr<<"Error in source file at line "<<curLine<<"\n";

  ctx.sources = new TSource[ctx.numSource];
  for(int i = 0; i < ctx.numSource; i++){
    curLine++;
    skip_comments(fin); fin >> ctx.sources[i].SourceType;
    switch(ctx.sources[i].SourceType){
    case 1:
      fin >> ctx.sources[i].Position.X >> ctx.sources[i].Position.Y >> ctx.sources[i].Position.Z;
      if(fin.fail()) cerr<<"Error in source file at line "<<curLine<<"\n";
      break;
    case 2:
      fin >> ctx.sources[i].ElemIdx;
      break;
    case 11:
      fin >> ctx.sources[i].ElemIdx;
      fin >> ctx.sources[i].Position.X >> ctx.sources[i].Position.Y >> ctx.sources[i].Position.Z;
      fin >> ctx.sources[i].IncAngle.X >> ctx.sources[i].IncAngle.Y >> ctx.sources[i].IncAngle.Z;
      break;
    case 12:
      fin >> ctx.sources[i].SurfTriNodes[0] >> ctx.sources[i].SurfTriNodes[1] >> ctx.sources[i].SurfTriNodes[2];
      break;
    default:
      cerr<<"\tUnknown source type near line: "<<i+1<<"\n"; return -1;
    }
    skip_comments(fin); fin >> ctx.sources[i].NumPhoton;
    total += ctx.sources[i].NumPhoton;
  }
  return total;
}

int Prepare_Source(SimContext& ctx){
  cerr << "Begin check source\n";

  int* nodeTriStart = new int[ctx.numNode+1]();
  for(int i = 1; i <= ctx.numTrig; i++)
    if(nodeTriStart[ctx.triNodes[i].N[0]] == 0)
      nodeTriStart[ctx.triNodes[i].N[0]] = i;
  for(int i = ctx.numNode-1; i >= 1; i--)
    if(nodeTriStart[i] == 0) nodeTriStart[i] = nodeTriStart[i+1];

  double h[4];
  for(int i = 0; i < ctx.numSource; i++){
    if(ctx.sources[i].SourceType == 1){
      ctx.sources[i].ElemIdx = -1;
      for(int e = 1; e <= ctx.numElem; e++){
        double minH = NO_INTERSECTION;
        h[0]=ctx.elems[e].TriNorm[0]*ctx.sources[i].Position.X
            +ctx.elems[e].TriNorm[4]*ctx.sources[i].Position.Y
            +ctx.elems[e].TriNorm[8]*ctx.sources[i].Position.Z+ctx.elems[e].TriNorm[12];
        if(h[0]<0) continue;
        if(h[0]<minH) minH=h[0];
        h[1]=ctx.elems[e].TriNorm[1]*ctx.sources[i].Position.X
            +ctx.elems[e].TriNorm[5]*ctx.sources[i].Position.Y
            +ctx.elems[e].TriNorm[9]*ctx.sources[i].Position.Z+ctx.elems[e].TriNorm[13];
        if(h[1]<0) continue;
        if(h[1]<minH) minH=h[1];
        h[2]=ctx.elems[e].TriNorm[2]*ctx.sources[i].Position.X
            +ctx.elems[e].TriNorm[6]*ctx.sources[i].Position.Y
            +ctx.elems[e].TriNorm[10]*ctx.sources[i].Position.Z+ctx.elems[e].TriNorm[14];
        if(h[2]<0) continue;
        if(h[2]<minH) minH=h[2];
        h[3]=ctx.elems[e].TriNorm[3]*ctx.sources[i].Position.X
            +ctx.elems[e].TriNorm[7]*ctx.sources[i].Position.Y
            +ctx.elems[e].TriNorm[11]*ctx.sources[i].Position.Z+ctx.elems[e].TriNorm[15];
        if(h[3]<0) continue;
        if(h[3]<minH) minH=h[3];
        ctx.sources[i].ElemIdx = e;
        if(minH < 1e-6){
          cerr<<"\tSource "<<i+1<<" too close to boundary of element "<<e<<"; nudging.\n";
          int* N = ctx.elemNodes[e].N;
          double CX=(ctx.nodes[N[0]].X+ctx.nodes[N[1]].X+ctx.nodes[N[2]].X+ctx.nodes[N[3]].X)/4.0;
          double CY=(ctx.nodes[N[0]].Y+ctx.nodes[N[1]].Y+ctx.nodes[N[2]].Y+ctx.nodes[N[3]].Y)/4.0;
          double CZ=(ctx.nodes[N[0]].Z+ctx.nodes[N[1]].Z+ctx.nodes[N[2]].Z+ctx.nodes[N[3]].Z)/4.0;
          double DX=CX-ctx.sources[i].Position.X, DY=CY-ctx.sources[i].Position.Y, DZ=CZ-ctx.sources[i].Position.Z;
          double dist=sqrt(DX*DX+DY*DY+DZ*DZ);
          if(dist<1e-6){ ctx.sources[i].Position.X=CX; ctx.sources[i].Position.Y=CY; ctx.sources[i].Position.Z=CZ; }
          else{ ctx.sources[i].Position.X+=1e-6*DX; ctx.sources[i].Position.Y+=1e-6*DY; ctx.sources[i].Position.Z+=1e-6*DZ; }
        }
        break;
      }
      if(ctx.sources[i].ElemIdx <= 0){
        cerr<<"\tSource "<<i+1<<" not within phantom.\n";
        delete[] nodeTriStart; return -1;
      }

    }else if(ctx.sources[i].SourceType == 12){
      for(int j=0;j<=1;j++)
        for(int k=j+1;k<=2;k++)
          if(ctx.sources[i].SurfTriNodes[j]>ctx.sources[i].SurfTriNodes[k])
            std::swap(ctx.sources[i].SurfTriNodes[j],ctx.sources[i].SurfTriNodes[k]);

      int tIdx=-1;
      int s0=ctx.sources[i].SurfTriNodes[0];
      for(int j=nodeTriStart[s0]; j<nodeTriStart[s0+1]; j++)
        if(ctx.sources[i].SurfTriNodes[1]==ctx.triNodes[j].N[1] &&
           ctx.sources[i].SurfTriNodes[2]==ctx.triNodes[j].N[2]){ tIdx=j; break; }
      if(tIdx==-1){ cerr<<"\tWrong source 1\n"; delete[] nodeTriStart; return -1; }
      if(ctx.triangles[tIdx].Num_Elem==2){ cerr<<"\tWrong source 2\n"; delete[] nodeTriStart; return -1; }

      ctx.sources[i].SurfTriIdx = tIdx;
      int e = ctx.sources[i].ElemIdx = ctx.triangles[tIdx].ElemIdx[0];
      // SoA TriNorm: face j → indices [j],[j+4],[j+8]
      if(ctx.sources[i].SurfTriNodes[0]==ctx.elemNodes[e].N[0]){
        if(ctx.sources[i].SurfTriNodes[1]==ctx.elemNodes[e].N[1]){
          if(ctx.sources[i].SurfTriNodes[2]==ctx.elemNodes[e].N[2])
            { ctx.sources[i].IncAngle.X=ctx.elems[e].TriNorm[0]; ctx.sources[i].IncAngle.Y=ctx.elems[e].TriNorm[4]; ctx.sources[i].IncAngle.Z=ctx.elems[e].TriNorm[8]; }
          else
            { ctx.sources[i].IncAngle.X=ctx.elems[e].TriNorm[1]; ctx.sources[i].IncAngle.Y=ctx.elems[e].TriNorm[5]; ctx.sources[i].IncAngle.Z=ctx.elems[e].TriNorm[9]; }
        }else
          { ctx.sources[i].IncAngle.X=ctx.elems[e].TriNorm[2]; ctx.sources[i].IncAngle.Y=ctx.elems[e].TriNorm[6]; ctx.sources[i].IncAngle.Z=ctx.elems[e].TriNorm[10]; }
      }else
        { ctx.sources[i].IncAngle.X=ctx.elems[e].TriNorm[3]; ctx.sources[i].IncAngle.Y=ctx.elems[e].TriNorm[7]; ctx.sources[i].IncAngle.Z=ctx.elems[e].TriNorm[11]; }
    }
  }
  cerr << "End check source\n";
  delete[] nodeTriStart;
  return 0;
}
