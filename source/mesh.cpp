#include "mesh.h"
#include "io_utils.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

using namespace std;

// Lexicographic order for deduplication during triangle enumeration.
inline bool operator<(TSimpleTri a, TSimpleTri b){
  return a.N[0] < b.N[0]
      || (a.N[0] == b.N[0] && a.N[1] < b.N[1])
      || (a.N[0] == b.N[0] && a.N[1] == b.N[1] && a.N[2] < b.N[2]);
}

int fem_read(const std::string& filename,
             bool         simpleOptic,
             int&         numNode,
             int&         numElem,
             int&         numMed,
             TNode*&      nodes,
             TElemNode*&  elemNodes,
             TElem*&      elems){

  ifstream fin(filename);
  if(!fin.good()){
    cerr << "\tCould not open file: " << filename << endl;
    return -1;
  }

  long long int curLine = 1;

  skip_comments(fin); fin >> numNode;
  if(fin.fail()){ cerr << "Error in fem file at line " << curLine << endl; }

  skip_comments(fin); fin >> numElem; curLine++;
  if(fin.fail()){ cerr << "Error in fem file at line " << curLine << endl; }

  if(!simpleOptic && numMed != numElem){
    cerr << "Error: Num Med should equal Num Elem for complex optic\n";
    return -1;
  }

  cerr << "\tNum Node: " << numNode << "\n\tNum Elem: " << numElem << endl;

  nodes     = new TNode    [numNode+1]();
  elemNodes = new TElemNode[numElem+1]();
  elems     = new TElem    [numElem+1];

  for(int i = 0; i <= numElem; i++)
    elemNodes[i].T[0] = elemNodes[i].T[1] = elemNodes[i].T[2] = elemNodes[i].T[3] = 0;

  for(int i = 1; i <= numNode; i++){
    skip_comments(fin);
    fin >> nodes[i].X >> nodes[i].Y >> nodes[i].Z;
    curLine++;
    if(fin.fail()){ cerr << "Error in fem file at line " << curLine << endl; }
  }

  for(int i = 1; i <= numElem; i++){
    skip_comments(fin);
    fin >> elemNodes[i].N[0] >> elemNodes[i].N[1]
        >> elemNodes[i].N[2] >> elemNodes[i].N[3]
        >> elems[i].MedIdx;

    for(int k = 0; k < 4; k++){
      if(elemNodes[i].N[k] < 1 || elemNodes[i].N[k] > numNode){
        cerr << "\tError in fem file at line: " << curLine << endl;
        return -1;
      }
    }
    curLine++;
    if(fin.fail()){ cerr << "Error at line " << curLine << endl; }

    if(simpleOptic){
      if(elems[i].MedIdx > numMed){
        cerr << "\tFEM file references medium " << elems[i].MedIdx
             << " but only " << numMed << " defined.\n";
        return -1;
      }
    }else{
      elems[i].MedIdx = i;
    }
    elemNodes[i].T[0] = elemNodes[i].T[1] = elemNodes[i].T[2] = elemNodes[i].T[3] = -1;

    // sort node indices of this element (small to large)
    for(int j = 0; j <= 2; j++)
      for(int k = j+1; k <= 3; k++)
        if(elemNodes[i].N[j] > elemNodes[i].N[k])
          std::swap(elemNodes[i].N[j], elemNodes[i].N[k]);
  }
  return 0;
}

int PreProcessor(int&        numElem,
                 int&        numNode,
                 int&        numMed,
                 int&        numBoundaryTrig,
                 int&        numTrig,
                 TNode*&     nodes,
                 TElemNode*& elemNodes,
                 TElem*&     elems,
                 TTriNode*&  triNodes,
                 TTriangle*& triangles,
                 int*&       boundaryTrigs){

  TSimpleTri* triangleList = new TSimpleTri[4*numElem+1];

  const int Ar[4][3] = {{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
  for(int i = 1; i <= numElem; i++){
    for(int j = 0; j < 4; j++){
      int k = (i-1)*4+j;
      triangleList[k].N[0] = elemNodes[i].N[Ar[j][0]];
      triangleList[k].N[1] = elemNodes[i].N[Ar[j][1]];
      triangleList[k].N[2] = elemNodes[i].N[Ar[j][2]];
      triangleList[k].NumConElem = 1;
      triangleList[k].ElemIdx[0] = i;
      triangleList[k].ElemIdx[1] = -1;
    }
  }
  triangleList[4*numElem].N[0] = triangleList[4*numElem].N[1] = triangleList[4*numElem].N[2] = numNode;

  sort(triangleList, triangleList + 4*numElem);

  int hole = 0, curTri = 0;
  numBoundaryTrig = 0;
  do{
    bool dup = (triangleList[curTri].N[0] == triangleList[curTri+1].N[0] &&
                triangleList[curTri].N[1] == triangleList[curTri+1].N[1] &&
                triangleList[curTri].N[2] == triangleList[curTri+1].N[2]);
    triangleList[hole].N[0] = triangleList[curTri].N[0];
    triangleList[hole].N[1] = triangleList[curTri].N[1];
    triangleList[hole].N[2] = triangleList[curTri].N[2];
    if(dup){
      triangleList[hole].NumConElem = 2;
      triangleList[hole].ElemIdx[0] = triangleList[curTri].ElemIdx[0];
      triangleList[hole].ElemIdx[1] = triangleList[curTri+1].ElemIdx[0];
      curTri += 2;
    }else{
      triangleList[hole].NumConElem = triangleList[curTri].NumConElem;
      triangleList[hole].ElemIdx[0] = triangleList[curTri].ElemIdx[0];
      triangleList[hole].ElemIdx[1] = triangleList[curTri].ElemIdx[1];
      curTri++;
      numBoundaryTrig++;
    }
    hole++;
  }while(curTri < 4*numElem);

  numTrig = hole;
  cerr << "\tNum_Trig: " << numTrig << endl;

  triNodes      = new TTriNode  [numTrig+1];
  triangles     = new TTriangle [numTrig+1];
  boundaryTrigs = new int       [numBoundaryTrig+1];

  int tempNumBoundaryTrig = 0;
  for(int p = 0; p < numTrig; p++){
    int L = p + 1;
    triNodes[L].N[0] = triangleList[p].N[0];
    triNodes[L].N[1] = triangleList[p].N[1];
    triNodes[L].N[2] = triangleList[p].N[2];

    if(triangleList[p].NumConElem == 1){
      triangles[L].Num_Elem    = 1;
      triangles[L].ElemIdx[0]  = triangleList[p].ElemIdx[0];
      triangles[L].ElemIdx[1]  = -1;
      triangles[L].BoundaryIdx = ++tempNumBoundaryTrig;
      boundaryTrigs[triangles[L].BoundaryIdx] = p+1;
    }else{
      triangles[L].Num_Elem   = 2;
      triangles[L].ElemIdx[0] = triangleList[p].ElemIdx[0];
      triangles[L].ElemIdx[1] = triangleList[p].ElemIdx[1];
    }

    for(int i = 0; i < triangles[L].Num_Elem; i++){
      if(triangles[L].ElemIdx[i] < 1 || triangles[L].ElemIdx[i] > numElem){
        cerr << "\tThe finite element mesh is not correct.\n";
        delete[] triangleList;
        return -1;
      }
      for(int j = 0; j <= 3; j++){
        if(elemNodes[triangles[L].ElemIdx[i]].T[j] == -1){
          elemNodes[triangles[L].ElemIdx[i]].T[j] = L;
          break;
        }
      }
    }
  }

  // Compute triangle normals, areas, and plane equations
  for(int i = 1; i <= numTrig; i++){
    int n1 = triNodes[i].N[0], n2 = triNodes[i].N[1], n3 = triNodes[i].N[2];
    double x1=nodes[n1].X, y1=nodes[n1].Y, z1=nodes[n1].Z;
    double x2=nodes[n2].X, y2=nodes[n2].Y, z2=nodes[n2].Z;
    double x3=nodes[n3].X, y3=nodes[n3].Y, z3=nodes[n3].Z;
    double ax=x2-x1, ay=y2-y1, az=z2-z1;
    double bx=x3-x1, by=y3-y1, bz=z3-z1;
    double nx = ay*bz-az*by, ny = az*bx-ax*bz, nz = ax*by-ay*bx;
    double len = sqrt(nx*nx+ny*ny+nz*nz);
    triNodes[i].Area = len / 2.0;
    triangles[i].X = nx/len;
    triangles[i].Y = ny/len;
    triangles[i].Z = nz/len;
    triangles[i].d = -(x1*triangles[i].X + y1*triangles[i].Y + z1*triangles[i].Z);
  }

  // Fill element normals in SoA TriNorm layout [nx0..3, ny0..3, nz0..3, d0..3]
  for(int i = 1; i <= numElem; i++){
    double x4 = (nodes[elemNodes[i].N[0]].X + nodes[elemNodes[i].N[1]].X +
                 nodes[elemNodes[i].N[2]].X + nodes[elemNodes[i].N[3]].X) / 4.0;
    double y4 = (nodes[elemNodes[i].N[0]].Y + nodes[elemNodes[i].N[1]].Y +
                 nodes[elemNodes[i].N[2]].Y + nodes[elemNodes[i].N[3]].Y) / 4.0;
    double z4 = (nodes[elemNodes[i].N[0]].Z + nodes[elemNodes[i].N[1]].Z +
                 nodes[elemNodes[i].N[2]].Z + nodes[elemNodes[i].N[3]].Z) / 4.0;

    for(int j = 0; j <= 3; j++){
      int tri = elemNodes[i].T[j];
      if(triangles[tri].Num_Elem == 1)
        elems[i].AdjElemIdx[j] = -triangles[tri].BoundaryIdx;
      else
        elems[i].AdjElemIdx[j] = (triangles[tri].ElemIdx[0] == i)
                                ? triangles[tri].ElemIdx[1]
                                : triangles[tri].ElemIdx[0];

      double nx = triangles[tri].X, ny = triangles[tri].Y;
      double nz = triangles[tri].Z, d  = triangles[tri].d;
      double t  = -(nx*x4 + ny*y4 + nz*z4 + d); // denominator is 1 (unit normal)

      if(t > 0){
        elems[i].TriNorm[j   ] = -nx;
        elems[i].TriNorm[j+ 4] = -ny;
        elems[i].TriNorm[j+ 8] = -nz;
        elems[i].TriNorm[j+12] = -d;
        elemNodes[i].Vol += triNodes[tri].Area * t / 3.0;
      }else if(t < 0){
        elems[i].TriNorm[j   ] = nx;
        elems[i].TriNorm[j+ 4] = ny;
        elems[i].TriNorm[j+ 8] = nz;
        elems[i].TriNorm[j+12] = d;
        elemNodes[i].Vol += -triNodes[tri].Area * t / 3.0;
      }
    }
  }

  cerr << "End Mesh prepare" << endl;
  delete[] triangleList;
  return 0;
}
