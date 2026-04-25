#include "mesh.h"
#include "io_utils.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;

// Lexicographic order for triangle deduplication.
static inline bool operator<(const TSimpleTri& a, const TSimpleTri& b){
  if (a.N[0] != b.N[0]) return a.N[0] < b.N[0];
  if (a.N[1] != b.N[1]) return a.N[1] < b.N[1];
  return a.N[2] < b.N[2];
}

// Faster skip_comments that avoids repeated std::ws/peek calls where possible
static inline void fast_skip_comments(std::istream& fin) {
  char c;
  while (fin >> std::ws && fin.peek() == '#') {
    std::string dummy;
    std::getline(fin, dummy);
  }
}

TiResult fem_read(const std::string& filename, SimContext& ctx){
  ifstream fin(filename);
  if(!fin.good()){ cout << "\tCould not open file: " << filename << "\n"; return std::unexpected("fem error"); }

  // Use a large buffer for reading
  std::vector<char> buffer(1024 * 1024);
  fin.rdbuf()->pubsetbuf(buffer.data(), buffer.size());

  fast_skip_comments(fin); fin >> ctx.numNode;
  fast_skip_comments(fin); fin >> ctx.numElem;

  if(!ctx.simpleOptic && ctx.numMed != ctx.numElem){
    cout << "Error: numMed must equal numElem for per-element optics.\n";
    return std::unexpected("fem error");
  }
  cout << "\tNum Node: " << ctx.numNode << "\n\tNum Elem: " << ctx.numElem << "\n";

  ctx.nodes.assign(ctx.numNode+1, TNode{});
  ctx.elemNodes.assign(ctx.numElem+1, TElemNode{});
  ctx.elems.resize(ctx.numElem+1);
  for(int i = 0; i <= ctx.numElem; i++)
    ctx.elemNodes[i].T[0] = ctx.elemNodes[i].T[1] = ctx.elemNodes[i].T[2] = ctx.elemNodes[i].T[3] = -1;

  for(int i = 1; i <= ctx.numNode; i++){
    fast_skip_comments(fin);
    fin >> ctx.nodes[i].X >> ctx.nodes[i].Y >> ctx.nodes[i].Z;
  }

  for(int i = 1; i <= ctx.numElem; i++){
    fast_skip_comments(fin);
    fin >> ctx.elemNodes[i].N[0] >> ctx.elemNodes[i].N[1]
        >> ctx.elemNodes[i].N[2] >> ctx.elemNodes[i].N[3]
        >> ctx.elems[i].MedIdx;

    if(ctx.simpleOptic){
      if(ctx.elems[i].MedIdx > ctx.numMed){
        cout << "\tFEM references medium " << ctx.elems[i].MedIdx
             << " but only " << ctx.numMed << " defined.\n";
        return std::unexpected("fem error");
      }
    }else{
      ctx.elems[i].MedIdx = i;
    }

    // Sort node indices small → large
    std::sort(ctx.elemNodes[i].N, ctx.elemNodes[i].N + 4);
  }
  return {};
}

TiResult PreProcessor(SimContext& ctx){
  int numElem = ctx.numElem, numNode = ctx.numNode;

  std::vector<TSimpleTri> tl(4 * numElem);
  const int Ar[4][3] = {{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
  for(int i = 1; i <= numElem; i++)
    for(int j = 0; j < 4; j++){
      int k = (i-1)*4+j;
      tl[k].N[0] = ctx.elemNodes[i].N[Ar[j][0]];
      tl[k].N[1] = ctx.elemNodes[i].N[Ar[j][1]];
      tl[k].N[2] = ctx.elemNodes[i].N[Ar[j][2]];
      tl[k].NumConElem = 1;
      tl[k].ElemIdx[0] = i; tl[k].ElemIdx[1] = -1;
    }

  // Use standard sort
  std::sort(tl.begin(), tl.end());

  int hole = 0, cur = 0, numBT = 0;
  int totalTl = 4 * numElem;
  while(cur < totalTl){
    bool hasNext = (cur + 1 < totalTl);
    bool dup = hasNext && (tl[cur].N[0]==tl[cur+1].N[0] && 
                           tl[cur].N[1]==tl[cur+1].N[1] && 
                           tl[cur].N[2]==tl[cur+1].N[2]);
    
    tl[hole].N[0] = tl[cur].N[0]; 
    tl[hole].N[1] = tl[cur].N[1]; 
    tl[hole].N[2] = tl[cur].N[2];
    
    if(dup){ 
      tl[hole].NumConElem=2; 
      tl[hole].ElemIdx[0]=tl[cur].ElemIdx[0];
      tl[hole].ElemIdx[1]=tl[cur+1].ElemIdx[0]; 
      cur+=2; 
    } else { 
      tl[hole].NumConElem=tl[cur].NumConElem; 
      tl[hole].ElemIdx[0]=tl[cur].ElemIdx[0];
      tl[hole].ElemIdx[1]=tl[cur].ElemIdx[1]; 
      cur++; 
      numBT++; 
    }
    hole++;
  }

  ctx.numTrig         = hole;
  ctx.numBoundaryTrig = numBT;
  cout << "\tNum_Trig: " << ctx.numTrig << "\n";

  ctx.triNodes.assign(ctx.numTrig+1, TTriNode{});
  ctx.triangles.assign(ctx.numTrig+1, TTriangle{});
  ctx.boundaryTrigs.assign(ctx.numBoundaryTrig+1, 0);

  int tempBT = 0;
  for(int p = 0; p < ctx.numTrig; p++){
    int L = p+1;
    ctx.triNodes[L].N[0] = tl[p].N[0];
    ctx.triNodes[L].N[1] = tl[p].N[1];
    ctx.triNodes[L].N[2] = tl[p].N[2];
    if(tl[p].NumConElem == 1){
      ctx.triangles[L].Num_Elem   = 1;
      ctx.triangles[L].ElemIdx[0] = tl[p].ElemIdx[0];
      ctx.triangles[L].ElemIdx[1] = -1;
      ctx.triangles[L].BoundaryIdx = ++tempBT;
      ctx.boundaryTrigs[tempBT] = p+1;
    }else{
      ctx.triangles[L].Num_Elem   = 2;
      ctx.triangles[L].ElemIdx[0] = tl[p].ElemIdx[0];
      ctx.triangles[L].ElemIdx[1] = tl[p].ElemIdx[1];
    }
    for(int i = 0; i < ctx.triangles[L].Num_Elem; i++){
      int eIdx = ctx.triangles[L].ElemIdx[i];
      for(int j = 0; j <= 3; j++)
        if(ctx.elemNodes[eIdx].T[j] == -1){
          ctx.elemNodes[eIdx].T[j] = L; break;
        }
    }
  }

  // Compute triangle normals / areas / plane equations
  for(int i = 1; i <= ctx.numTrig; i++){
    int n1=ctx.triNodes[i].N[0], n2=ctx.triNodes[i].N[1], n3=ctx.triNodes[i].N[2];
    double x1=ctx.nodes[n1].X,y1=ctx.nodes[n1].Y,z1=ctx.nodes[n1].Z;
    double x2=ctx.nodes[n2].X,y2=ctx.nodes[n2].Y,z2=ctx.nodes[n2].Z;
    double x3=ctx.nodes[n3].X,y3=ctx.nodes[n3].Y,z3=ctx.nodes[n3].Z;
    double ax=x2-x1,ay=y2-y1,az=z2-z1, bx=x3-x1,by=y3-y1,bz=z3-z1;
    double nx=ay*bz-az*by, ny=az*bx-ax*bz, nz=ax*by-ay*bx;
    double len = sqrt(nx*nx+ny*ny+nz*nz);
    ctx.triNodes[i].Area = len/2.0;
    ctx.triangles[i].X = nx/len; ctx.triangles[i].Y = ny/len; ctx.triangles[i].Z = nz/len;
    ctx.triangles[i].d = -(x1*ctx.triangles[i].X + y1*ctx.triangles[i].Y + z1*ctx.triangles[i].Z);
  }

  // Fill element SoA TriNorm [nx0..3, ny0..3, nz0..3, d0..3]
  for(int i = 1; i <= numElem; i++){
    int* N = ctx.elemNodes[i].N;
    double x4=(ctx.nodes[N[0]].X+ctx.nodes[N[1]].X+ctx.nodes[N[2]].X+ctx.nodes[N[3]].X)/4.0;
    double y4=(ctx.nodes[N[0]].Y+ctx.nodes[N[1]].Y+ctx.nodes[N[2]].Y+ctx.nodes[N[3]].Y)/4.0;
    double z4=(ctx.nodes[N[0]].Z+ctx.nodes[N[1]].Z+ctx.nodes[N[2]].Z+ctx.nodes[N[3]].Z)/4.0;

    ctx.elemNodes[i].Vol = 0;
    for(int j = 0; j <= 3; j++){
      int tri = ctx.elemNodes[i].T[j];
      ctx.elems[i].AdjElemIdx[j] = (ctx.triangles[tri].Num_Elem==1)
        ? -ctx.triangles[tri].BoundaryIdx
        : (ctx.triangles[tri].ElemIdx[0]==i ? ctx.triangles[tri].ElemIdx[1]
                                            : ctx.triangles[tri].ElemIdx[0]);
      double nx=ctx.triangles[tri].X, ny=ctx.triangles[tri].Y;
      double nz=ctx.triangles[tri].Z, d=ctx.triangles[tri].d;
      double t = -(nx*x4+ny*y4+nz*z4+d);
      if(t>0){ ctx.elems[i].TriNorm[j]=-nx; ctx.elems[i].TriNorm[j+4]=-ny;
               ctx.elems[i].TriNorm[j+8]=-nz; ctx.elems[i].TriNorm[j+12]=-d;
               ctx.elemNodes[i].Vol += ctx.triNodes[tri].Area*t/3.0; }
      else {   ctx.elems[i].TriNorm[j]=nx; ctx.elems[i].TriNorm[j+4]=ny;
               ctx.elems[i].TriNorm[j+8]=nz; ctx.elems[i].TriNorm[j+12]=d;
               ctx.elemNodes[i].Vol += -ctx.triNodes[tri].Area*t/3.0; }
    }
  }

  cout << "End Mesh prepare\n";
  return {};
}

