#include "source_io.h"
#include "io_utils.h"
#include "constants.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

long long int ReadSource(const std::string& filename,
                         int&      numSource,
                         TSource*& sources){
  long long int l_totalPhoton = 0;
  ifstream fin(filename);
  if(!fin.good()){
    cerr << "\tCould not open file: " << filename << endl;
    return -1;
  }

  int curLine = 0;
  skip_comments(fin);
  fin >> numSource; curLine++;
  if(fin.fail()){ cerr << "Error in source file at line " << curLine << endl; }

  sources = new TSource[numSource];

  for(int i = 0; i < numSource; i++){
    curLine++;
    skip_comments(fin);
    fin >> sources[i].SourceType;
    switch(sources[i].SourceType){
    case 1: // isotropic point source
      fin >> sources[i].Position.X >> sources[i].Position.Y >> sources[i].Position.Z;
      if(fin.fail()){ cerr << "Error in source file at line " << curLine << endl; }
      break;
    case 2: // isotropic region source
      fin >> sources[i].ElemIdx;
      break;
    case 11: // pencil beam
      fin >> sources[i].ElemIdx;
      fin >> sources[i].Position.X >> sources[i].Position.Y >> sources[i].Position.Z;
      fin >> sources[i].IncAngle.X >> sources[i].IncAngle.Y >> sources[i].IncAngle.Z;
      break;
    case 12: // surface triangle region
      fin >> sources[i].SurfTriNodes[0] >> sources[i].SurfTriNodes[1] >> sources[i].SurfTriNodes[2];
      break;
    default:
      cerr << "\tUnknown source type near line: " << i+1 << endl;
      return -1;
    }

    skip_comments(fin);
    fin >> sources[i].NumPhoton;
    l_totalPhoton += sources[i].NumPhoton;
  }
  return l_totalPhoton;
}

int Prepare_Source(int&        numNode,
                   int&        numElem,
                   int&        numTrig,
                   int&        numSource,
                   TSource*&   sources,
                   TElem*&     elems,
                   TNode*&     nodes,
                   TElemNode*& elemNodes,
                   TTriNode*&  triNodes,
                   TTriangle*& triangles){
  cerr << "Begin check source" << endl;

  // Build index: nodeTriStartIdx[n] = first triangle whose N[0] == n
  int* nodeTriStartIdx = new int[numNode+1]();
  for(int i = 1; i <= numTrig; i++){
    if(nodeTriStartIdx[triNodes[i].N[0]] == 0)
      nodeTriStartIdx[triNodes[i].N[0]] = i;
  }
  for(int i = numNode-1; i >= 1; i--)
    if(nodeTriStartIdx[i] == 0)
      nodeTriStartIdx[i] = nodeTriStartIdx[i+1];

  double height[4];
  for(int i = 0; i < numSource; i++){
    if(sources[i].SourceType == 1){
      // Locate the element containing the point source
      sources[i].ElemIdx = -1;
      for(int e = 1; e <= numElem; e++){
        double Min_H = NO_INTERSECTION;
        // SoA TriNorm: face j at indices [j],[j+4],[j+8],[j+12]
        height[0] = elems[e].TriNorm[0]*sources[i].Position.X
                  + elems[e].TriNorm[4]*sources[i].Position.Y
                  + elems[e].TriNorm[8]*sources[i].Position.Z
                  + elems[e].TriNorm[12];
        if(height[0] < 0) goto L_ProcSource_Outside;
        if(height[0] < Min_H) Min_H = height[0];

        height[1] = elems[e].TriNorm[1]*sources[i].Position.X
                  + elems[e].TriNorm[5]*sources[i].Position.Y
                  + elems[e].TriNorm[9]*sources[i].Position.Z
                  + elems[e].TriNorm[13];
        if(height[1] < 0) goto L_ProcSource_Outside;
        if(height[1] < Min_H) Min_H = height[1];

        height[2] = elems[e].TriNorm[2]*sources[i].Position.X
                  + elems[e].TriNorm[6]*sources[i].Position.Y
                  + elems[e].TriNorm[10]*sources[i].Position.Z
                  + elems[e].TriNorm[14];
        if(height[2] < 0) goto L_ProcSource_Outside;
        if(height[2] < Min_H) Min_H = height[2];

        height[3] = elems[e].TriNorm[3]*sources[i].Position.X
                  + elems[e].TriNorm[7]*sources[i].Position.Y
                  + elems[e].TriNorm[11]*sources[i].Position.Z
                  + elems[e].TriNorm[15];
        if(height[3] < 0) goto L_ProcSource_Outside;
        if(height[3] < Min_H) Min_H = height[3];

        sources[i].ElemIdx = e;
        if(Min_H < 1e-6){
          cerr << "\tSource " << i+1 << " too close to boundary of element " << e
               << "; nudging toward centroid.\n";
          double CX = (nodes[elemNodes[e].N[0]].X + nodes[elemNodes[e].N[1]].X +
                       nodes[elemNodes[e].N[2]].X + nodes[elemNodes[e].N[3]].X) / 4.0;
          double CY = (nodes[elemNodes[e].N[0]].Y + nodes[elemNodes[e].N[1]].Y +
                       nodes[elemNodes[e].N[2]].Y + nodes[elemNodes[e].N[3]].Y) / 4.0;
          double CZ = (nodes[elemNodes[e].N[0]].Z + nodes[elemNodes[e].N[1]].Z +
                       nodes[elemNodes[e].N[2]].Z + nodes[elemNodes[e].N[3]].Z) / 4.0;
          double DX = CX - sources[i].Position.X;
          double DY = CY - sources[i].Position.Y;
          double DZ = CZ - sources[i].Position.Z;
          double dist = sqrt(DX*DX + DY*DY + DZ*DZ);
          if(dist < 1e-6){
            sources[i].Position.X = CX;
            sources[i].Position.Y = CY;
            sources[i].Position.Z = CZ;
          }else{
            sources[i].Position.X += 1e-6 * DX;
            sources[i].Position.Y += 1e-6 * DY;
            sources[i].Position.Z += 1e-6 * DZ;
          }
        }
        break;
      L_ProcSource_Outside:;
      }
      if(sources[i].ElemIdx <= 0){
        cerr << "\tSource " << i+1 << " is not within the phantom.\n";
        delete[] nodeTriStartIdx;
        return -1;
      }

    }else if(sources[i].SourceType == 12){
      // Surface triangle source: sort nodes, find triangle, assign element + normal
      for(int j = 0; j <= 1; j++)
        for(int k = j+1; k <= 2; k++)
          if(sources[i].SurfTriNodes[j] > sources[i].SurfTriNodes[k])
            std::swap(sources[i].SurfTriNodes[j], sources[i].SurfTriNodes[k]);

      int Source_TriIdx = -1;
      int start = nodeTriStartIdx[sources[i].SurfTriNodes[0]];
      int end   = nodeTriStartIdx[sources[i].SurfTriNodes[0]+1];
      for(int j = start; j < end; j++){
        if(sources[i].SurfTriNodes[1] == triNodes[j].N[1] &&
           sources[i].SurfTriNodes[2] == triNodes[j].N[2]){
          Source_TriIdx = j;
          break;
        }
      }
      if(Source_TriIdx == -1){ cerr << "\tWrong source 1\n"; delete[] nodeTriStartIdx; return -1; }
      if(triangles[Source_TriIdx].Num_Elem == 2){ cerr << "\tWrong source 2\n"; delete[] nodeTriStartIdx; return -1; }

      sources[i].SurfTriIdx = Source_TriIdx;
      sources[i].ElemIdx    = triangles[Source_TriIdx].ElemIdx[0];

      // Pick the matching face normal (SoA layout: face j at [j],[j+4],[j+8])
      int e = sources[i].ElemIdx;
      if(sources[i].SurfTriNodes[0] == elemNodes[e].N[0]){
        if(sources[i].SurfTriNodes[1] == elemNodes[e].N[1]){
          if(sources[i].SurfTriNodes[2] == elemNodes[e].N[2]){
            sources[i].IncAngle.X = elems[e].TriNorm[0];
            sources[i].IncAngle.Y = elems[e].TriNorm[4];
            sources[i].IncAngle.Z = elems[e].TriNorm[8];
          }else{
            sources[i].IncAngle.X = elems[e].TriNorm[1];
            sources[i].IncAngle.Y = elems[e].TriNorm[5];
            sources[i].IncAngle.Z = elems[e].TriNorm[9];
          }
        }else{
          sources[i].IncAngle.X = elems[e].TriNorm[2];
          sources[i].IncAngle.Y = elems[e].TriNorm[6];
          sources[i].IncAngle.Z = elems[e].TriNorm[10];
        }
      }else{
        sources[i].IncAngle.X = elems[e].TriNorm[3];
        sources[i].IncAngle.Y = elems[e].TriNorm[7];
        sources[i].IncAngle.Z = elems[e].TriNorm[11];
      }
    }
  }

  cerr << "End check source" << endl;
  delete[] nodeTriStartIdx;
  return 0;
}
