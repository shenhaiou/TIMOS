//-------------------------------------------------
// Copyright: 
// Biomedical Imaging Division 
// School of Biomedical Engineering and Sciences
// Virginia Tech
// 2011
//----------------------------------------------------

// To compile on macOS (M1):
// clang++ -O3 -std=c++17 -march=native -ffast-math -funroll-loops -flto \
//         -DHAS_ACCELERATE -framework Accelerate -o timos TIMOS.cpp

#define _REENTRANT
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#  include <arm_neon.h>
#  define TIMOS_NEON 1
#endif
#include <time.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <pthread.h>

#include <set>
#include <deque>
#include <algorithm>

#include <sys/stat.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <cmath>
#include <dispatch/dispatch.h>
#include "rng.h"

#include "timos.h"

using namespace std;

// Batched RNG buffer sizes. Buffers are refilled on exhaustion.
#define MAXRANDNUM  256
#define MAXRANDNUM3 768

// Fill buf[0..n) with uniform doubles in [lo, hi)
static inline void rng_fill_uniform(timos::Xoshiro256ss& rng, double* buf, int n, double lo, double hi) {
  const double scale = (hi - lo) * 0x1.0p-53;
  for (int i = 0; i < n; i++)
    buf[i] = (rng() >> 11) * scale + lo;
}

// Fill buf[0..n) with log(u) for u uniform in (0,1); stored value is <= 0.
// Step-size sampling uses -buf[i], yielding an exponential variate.
static inline void rng_fill_neglog(timos::Xoshiro256ss& rng, double* buf, int n) {
  for (int i = 0; i < n; i++) {
    double u = (rng() >> 11) * 0x1.0p-53;
    if (u <= 0.0) u = 1e-300;
    buf[i] = std::log(u);
  }
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Definition of global constant variables
// The name format for a global constant should be all:
// UP_CASE

// The speed of light in vacuum is in millimeter per nano-second.
const double LIGHT_SPEED = 2.99792458E2; // mm/ns
const double INV_LIGHT_SPEED = 1.0/LIGHT_SPEED;

// Similar to the definition in MCML but with higher threshold
const double G_COS_0_D   = 1.0 - 1.0E-14; 
// cos(0.0000001) > G_COS_0_D
const double G_COS_90_D  = 1.0E-7;
// cos(M_PI/2.0-0.00000001) < G_COS_90_D


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Variables for thread mutex — use os_unfair_lock (Apple-native, lower overhead than pthread_mutex)
#include <os/lock.h>
static os_unfair_lock Result_Lock = OS_UNFAIR_LOCK_INIT;
static os_unfair_lock Source_Lock = OS_UNFAIR_LOCK_INIT;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Definition of global variables
// The format for the name of a global variables is:
// g_VariableFormat 

// -----------------------------------------------------------
double g_TimeStep                 = 0.1; // 0.1 ns
double g_InvTimeStep              = 1/g_TimeStep;
double g_InvLightSpeedMutTimeStep = INV_LIGHT_SPEED*g_InvTimeStep;
int    g_NumTimeStep              = 50;
// -----------------------------------------------------------

// This is the global variable that stores the number of 
//      photon-tetrahedron intersection tests in the simulation.
long long int g_NumIntersections = 0;
long long int g_NumSteps = 0;


bool g_SimpleOptic     = true;
// true:  optical parameters by region.
//        In this setting, tetrahedra within the same region have  
//        the same optical parameter.
// 
// false: optical parameters by tetrahedron.
//        In this setting, each tetrahedron can have different 
//        optical parameter.
//        This feature is not fully tested.


int  g_UniformBoundary = 0; 
// 0: Uniform Boundary with one environment refractive index
// 1: Martched boundary (no reflection)
// The next two settings are not implemented yet.
// 2: Each region has different environment refractive index
// 3: Each surface triangle has different environment refractive index 


bool g_TimeDomain      = false;
// true  : time domain simulation
// false : normal mc simulation


int g_StartRandIdx = 1;
// The random number generator used in this program is Intel MT2203,     
//     which can have 1024 different independent streams. 
//     We can use this to distribute the program on a cluster 
//     such that each thread has different random stream.  
//     For example, we can run the program on 2 machines, 
//     each machine has 8 threads. In this case, the first machine
//     can have g_StartRandIdx = 1, and the second machine can have
//     g_StartRandIdx = 9.  


// This global variable controls how TIM-OS save the result
// This feature is not implimented yet.
int g_InternalResultFormat;
//  0 -> Ignore internal fluence quantification
//  1 -> Save fluence to tetrahedron
//  2 -> Save fluence to voxel grid (dx, dy, dz may not the same)
//       To use this, user also need to supply the ROI (minx, miny, minz, maxx, maxy, maxz)
//  3 -> Save fluence to cylindrical grid (dr, dz)
//  4 -> Save fluence to spherical grid (dr)


// can be set in optical parameter file
double g_EnvRefIdx = 1.0;
// Environment refractive index

int g_NumMed;    // Number of medium
int g_NumNode;   // Number of nodes
int g_NumElem;   // Number of tetrahedrons
int g_NumTrig;   // Number of triangles
int g_NumBoundaryTrig; // Number of boundary triangles
int g_NumSource; // Number of sources
int g_NumThread; // Number of threads


TMedOptic * g_MedOptic;  
// If g_SimpleOptic = true
//    The size of g_MedOptic is g_NumMed
// Else 
//    The size of g_MedOptic is g_NumMed = g_NumElem

TNode     * g_Nodes;     
// Node list. The index is start from 1.

TTriNode  * g_TriNodes;  
// Node index list for triangles

TElemNode * g_ElemNodes; 
// Node index list for g_Elems.

TElem     * g_Elems;     
// Tetrahedron element list. The Index is start from 1.

int       * g_BoundaryTrigs;     
// The boundary triangle index to triangle index. The Index is start form 1.

TTriangle * g_Triangles; 
// Triangle list. The index is start from 1.


double    * g_SurfMeas;  
// The surface fluence #photon/mm^2

double   ** g_TimeSurfMeas;
// Time domain suface fluence

double    * g_Absorption;   
// g_Absorption: Internal fluence 

double   ** g_TimeAbsorption;
// Time domain internal fluence

TSource   * g_Sources;
int         g_SourceIdx = 0;


long long int g_SimedPhoton = 0;
long long int g_TotalPhoton  = 0;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// operator to sort the TSimpleTri
// The reason is that we need to find out how many different triangles
//     in the finite element file. By define an order, we can then sort 
//     the triangles. 
// The algorithm is that we first add all triangles into a list.
//     Then, we sort the triangle list.
//     Then, we do an one-pass scan to find out the duplicated triangles,
//     which are the internal triangles.
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
inline bool operator < (TSimpleTri a, TSimpleTri b){
  if( a.N[0]<b.N[0] || 
      (a.N[0]==b.N[0] && a.N[1]<b.N[1]) || 
      (a.N[0]==b.N[0] && a.N[1]==b.N[1] && a.N[2]<b.N[2])){
    return true;
  }else{
    return false;
  }
}

// Skip whitespace then any number of '#'-to-end-of-line comment lines.
// Call before every logical token read in the three file parsers.
static inline void skip_comments(std::istream& fin) {
  fin >> std::ws;
  while(fin.good() && fin.peek() == '#') {
    std::string line;
    std::getline(fin, line);
    fin >> std::ws;
  }
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++
//
// Read finite element mesh from file
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++
int fem_read(char      *  filename,
	     bool         simpleOptic,
	     int       &  numNode, 
	     int       &  numElem,
	     int       &  numMed,
	     TNode     *& nodes, 
	     TElemNode *& elemNodes, 
	     TElem     *& elems){

  ifstream fin;
  fin.open(filename, ifstream::in);
  if(!fin.good()){
    cerr << "\tCould not open file: " << filename << endl;
    return -1;
  }

  // The format of the data file:
  // 1> numNode (int): index starts from 1
  // 2> numElem (int): index starts from 1
  // 3> List of nodes
  //    x, y, z (double)
  // 4> List of elems
  //    n1, n2, n3, n4, med (int)
  // The largest tested example contains more than 10,000,000 tetrahedra

  // Read the finite element mesh data.

  long long int curLine = 1;

  skip_comments(fin);
  fin >> numNode;
  if(fin.fail()){cerr <<"Error in fem file at line " << curLine <<endl;}

  skip_comments(fin);
  fin >> numElem;
  curLine++;
  if(fin.fail()){cerr <<"Error in fem file at line " << curLine <<endl;}
  
  if(!simpleOptic){
    if(numMed != numElem){
      cerr <<"Error in fem file at line " << curLine <<endl;
      cerr << "\tNum Med should equal Num Elem for complex optic";
      return -1;
    }
  }
  
  cerr << "\tNum Node: " << numNode<< endl;
  cerr << "\tNum Elem: " << numElem<< endl;
  
  nodes     = new TNode [numNode+1];
  elemNodes = new TElemNode[numElem+1];
  elems     = new TElem [numElem+1];
  

  for (int i=0; i<=numNode; i++){
    nodes[i].X = nodes[i].Y = nodes[i].Z = 0.;
  }

  for (int i=0; i<=numElem; i++){
    elemNodes[i].N[0] = elemNodes[i].N[1] = elemNodes[i].N[2] = elemNodes[i].N[3] = 0;
    elemNodes[i].T[0] = elemNodes[i].T[1] = elemNodes[i].T[2] = elemNodes[i].T[3] = 0;
    elemNodes[i].Vol = 0.;
  }


  for (int i=1; i<=numNode; i++){
    skip_comments(fin);
    fin >> nodes[i].X;
    fin >> nodes[i].Y;
    fin >> nodes[i].Z;
    curLine++;
    if(fin.fail()){cerr <<"Error in fem file at line " << curLine <<endl;}
  }

  for (int i=1; i<=numElem; i++){
    skip_comments(fin);
    fin >> elemNodes[i].N[0];
    fin >> elemNodes[i].N[1];
    fin >> elemNodes[i].N[2];
    fin >> elemNodes[i].N[3];
    fin >> elems[i].MedIdx;

    if(elemNodes[i].N[0]<1 || elemNodes[i].N[0]>numNode || 
       elemNodes[i].N[1]<1 || elemNodes[i].N[1]>numNode || 
       elemNodes[i].N[2]<1 || elemNodes[i].N[2]>numNode || 
       elemNodes[i].N[3]<1 || elemNodes[i].N[3]>numNode){
      cerr << "\tError in fem file at line: " << curLine << endl;
      return -1;
    }

    curLine++;
    if(fin.fail()){cerr <<"Error at line " << curLine <<endl;}

    if(simpleOptic){
      if(elems[i].MedIdx>numMed){
	cerr << "\tError in fem file at line: " << curLine << endl;
	cerr << "\tFEM file contains more Meds than optical parameter file: ";
	cerr << "\t" << elemNodes[i].N[0] << " " 
	     << elemNodes[i].N[1] << " " 
	     << elemNodes[i].N[2] << " " 
	     << elemNodes[i].N[3] << " " 
	     << elems[i].MedIdx << endl;
	return -1;
      }
    }else{
      elems[i].MedIdx = i;
    }
    //
    elemNodes[i].T[0] = elemNodes[i].T[1] = elemNodes[i].T[2] = elemNodes[i].T[3] = -1;
    
    // sort the Nodes by their index from small to large
    int temp;
    for (int j=0; j<=2; j++){
      for (int k=j+1; k<=3; k++){
	if(elemNodes[i].N[j]>elemNodes[i].N[k]){
	  temp = elemNodes[i].N[j];
	  elemNodes[i].N[j]=elemNodes[i].N[k];
	  elemNodes[i].N[k]=temp;
	}
      }
    }
  }    
  fin.close();
  return 0;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Update the data structures needed for MC from mesh file
//            : elems
//            : elemNodes
//            : TriNodes
//            : triangles
//            : BoundaryTrigs
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int PreProcessor(int        & numElem,
		 int        & numNode,
		 int        & numMed,
		 int        & numBoundaryTrig,
		 int        & numTrig,
		 TNode     *& nodes, 
		 TElemNode *& elemNodes, 
		 TElem     *& elems, 
		 TTriNode  *& triNodes,
		 TTriangle *& triangles,
		 int       *& boundaryTrigs){

  TSimpleTri * triangleList = new TSimpleTri [4*numElem+1]; 

  int Ar[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};

  for (int i=1; i<=numElem; i++){
    for (int j=0; j<4; j++){
      int k=(i-1)*4+j;
      triangleList[k].N[0] = elemNodes[i].N[Ar[j][0]];
      triangleList[k].N[1] = elemNodes[i].N[Ar[j][1]];
      triangleList[k].N[2] = elemNodes[i].N[Ar[j][2]];
      triangleList[k].NumConElem = 1;
      triangleList[k].ElemIdx[0] = i;
      triangleList[k].ElemIdx[1] = -1;
    }
  }
  triangleList[4*numElem].N[0] = numNode;
  triangleList[4*numElem].N[1] = numNode;
  triangleList[4*numElem].N[2] = numNode;

  sort(triangleList, triangleList+4*numElem);

  int hole   = 0;
  int curTri = 0;
  numBoundaryTrig = 0;
  do{
    if(triangleList[curTri].N[0] == triangleList[curTri+1].N[0] &&
       triangleList[curTri].N[1] == triangleList[curTri+1].N[1] &&
       triangleList[curTri].N[2] == triangleList[curTri+1].N[2]){
      triangleList[hole].N[0] = triangleList[curTri].N[0];
      triangleList[hole].N[1] = triangleList[curTri].N[1];
      triangleList[hole].N[2] = triangleList[curTri].N[2];
      triangleList[hole].NumConElem = 2;
      triangleList[hole].ElemIdx[0] = triangleList[curTri].ElemIdx[0];
      triangleList[hole].ElemIdx[1] = triangleList[curTri+1].ElemIdx[0];
      curTri+=2;
    }else{
      triangleList[hole].N[0] = triangleList[curTri].N[0];
      triangleList[hole].N[1] = triangleList[curTri].N[1];
      triangleList[hole].N[2] = triangleList[curTri].N[2];
      triangleList[hole].NumConElem = triangleList[curTri].NumConElem;
      triangleList[hole].ElemIdx[0] = triangleList[curTri].ElemIdx[0];
      triangleList[hole].ElemIdx[1] = triangleList[curTri].ElemIdx[1];
      curTri++;
      numBoundaryTrig++;
    }
    hole++;
  }while(curTri<4*numElem);

  // Now we have all the triangles. We need to update the triangle information 
  // in the elems, mainly the adject Elem information.
  //  cerr << "2: " << endl;
  int L_Idx = 1;
  numTrig  = hole;
  cerr << "\tNum_Trig: " << numTrig << endl;
  triNodes          = new TTriNode[numTrig+1];
  triangles         = new TTriangle [numTrig+1];
  boundaryTrigs     = new int [numBoundaryTrig+1];

  int tempNumBoundaryTrig = 0;

  for (int p_tri = 0; p_tri<numTrig; p_tri++){
    triNodes[L_Idx].N[0] = triangleList[p_tri].N[0];
    triNodes[L_Idx].N[1] = triangleList[p_tri].N[1];
    triNodes[L_Idx].N[2] = triangleList[p_tri].N[2];

    if(triangleList[p_tri].NumConElem==1){
      triangles[L_Idx].Num_Elem = 1;
      triangles[L_Idx].ElemIdx[0] = triangleList[p_tri].ElemIdx[0];
      triangles[L_Idx].ElemIdx[1] = -1;
      triangles[L_Idx].BoundaryIdx = (++tempNumBoundaryTrig);
      boundaryTrigs[triangles[L_Idx].BoundaryIdx] = p_tri+1;
    }else{
      triangles[L_Idx].Num_Elem = 2;
      triangles[L_Idx].ElemIdx[0] = triangleList[p_tri].ElemIdx[0];
      triangles[L_Idx].ElemIdx[1] = triangleList[p_tri].ElemIdx[1];
    }

    for (int i=0; i<triangles[L_Idx].Num_Elem; i++){
      for (int j=0; j<=3; j++){
	if(triangles[L_Idx].ElemIdx[i]<1 || triangles[L_Idx].ElemIdx[i]>numElem){
	  cerr << "\tThe finite element mesh is not correct.";
	    return -1;
	}
	if (elemNodes[triangles[L_Idx].ElemIdx[i]].T[j]==-1){
	  elemNodes[triangles[L_Idx].ElemIdx[i]].T[j] = L_Idx;
	  break;
	}
      }
    }
    L_Idx++;

  }

  //  cerr << NumBoundaryTrig << " " << T_NumBoundaryTrig << endl;

  // Update three information for each triangle:
  // 1) normal vector
  //    The nodes for the triangle is n1=<x1, y1, z1>, n2=<x2, y2, z2>, n3=<x3, y3, z3>.
  //    A = n2 - n1 = <x2-x1, y2-y1, z2-z1>, and B = n3 - n1 = = <x3-x1, y3-y1, z3-z1>
  //    The normal vector N = A X B
  //            | i  j  k|
  //    A X B = |ax ay az| = i(ay bz - az by) + j(az bx - ax bz) + k(ax by - ay bx)
  //            |bx by bz|
  //    Then the unit normal vector is N' = N/|N|
  // 2) area
  //    Since |A X B| = |A||B|sin(<AB), we can have:
  //    area = |N|/2
  // 3) plane equation a x + b y + c z + d == 0
  // for any point (x,y,z) in the plane, (x-x1) nx + (y-y1) ny + (z-z1) nz = 0
  // hence:
  // a  x + b  y + c  z + d == 0
  // nx x + ny y + nz z - (x1 nx + y1 ny + z1 nz) == 0
  // d = - (x1 nx + y1 ny + z1 nz)

  int n1, n2, n3;
  double x1, y1, z1, x2, y2, z2, x3, y3, z3;
  double ax, ay, az, bx, by, bz;
  double nx, ny, nz;
  double temp;
  for (int i=1; i<=numTrig; i++){
    n1 = triNodes[i].N[0];
    n2 = triNodes[i].N[1];
    n3 = triNodes[i].N[2];
    
    x1 = nodes[n1].X; y1 = nodes[n1].Y; z1 = nodes[n1].Z;
    x2 = nodes[n2].X; y2 = nodes[n2].Y; z2 = nodes[n2].Z;
    x3 = nodes[n3].X; y3 = nodes[n3].Y; z3 = nodes[n3].Z;

    ax = x2 - x1; ay = y2 - y1; az = z2 - z1;
    bx = x3 - x1; by = y3 - y1; bz = z3 - z1;

    nx = ay * bz - az * by;
    ny = az * bx - ax * bz;
    nz = ax * by - ay * bx;

    temp = sqrt(nx * nx + ny * ny + nz * nz);
    
    triNodes[i].Area = temp / 2.0; 

    triangles[i].X = nx/temp;
    triangles[i].Y = ny/temp;
    triangles[i].Z = nz/temp;
    triangles[i].d = -(x1*triangles[i].X + y1*triangles[i].Y + z1*triangles[i].Z);
  }

  // Update the information in elems
  // 1) The four normal vectors for four triangles (the vector needs to point in the tetrahedron
  //    The normal vector of the triangle is in the triangle structure. 
  //    But it may point to outside of the tetrahedron.
  //    For example, the three nodes of the triangle is n1, n2, n3 
  //        and the center of the tetrahedron is n4=<x4, y4, z4>.
  //    The normal vector is N=<nx, ny, nz>. 
  //    The plane function for the triangle is ax + by + cz + d = 0;
  //    The line  n4 + t N w will have a cross point on the plane:
  //    a (x4 + t nx) + b (y4 + t ny) + c (z4 + t nz) + d = 0
  //    If t < 0, then the normal vector is point inside the tetrahegron. 
  //    If t = 0. then error
  //    If t >0, then the normal vector is point outside the tetrahegron.
  // 2) Vol of the tetrahedron
  //    Based on the solution of 1), the vol of the tetrahedron can be easily computed.
  //    The area of the triangle is already computed.
  //    The height of the tetrahedron is |t|.
  //    Hence the vol of the tetrahedron is Area * |t| / 3

  // T[0] to T[3] are the triangle index list of the tetrahedron
  // T[0] should be nodes 0, 1, 2. Hence the node not in the triangle is 3.
  // T[1] should be nodes 0, 1, 3. Hence the node not in the triangle is 2.
  // T[2] should be nodes 0, 2, 3. Hence the node not in the triangle is 1.
  // T[3] should be nodes 1, 2, 3. Hence the node not in the triangle is 0.

  //  double a, b, c, d; // the plane function for the triangle
  // <nx, ny, nz> is the normal vector of the triangle
  //  double x4, y4, z4; // the center node of the tetrahedron

  //  double t;

  for (int i=1; i<=numElem; i++){
    // for each tetrahedron, there are four triangle planes
    //    elems[i].Vol = 0;

    double x4 = (nodes[elemNodes[i].N[0]].X+nodes[elemNodes[i].N[1]].X+
		 nodes[elemNodes[i].N[2]].X+nodes[elemNodes[i].N[3]].X)/4;
    double y4 = (nodes[elemNodes[i].N[0]].Y+nodes[elemNodes[i].N[1]].Y+
		 nodes[elemNodes[i].N[2]].Y+nodes[elemNodes[i].N[3]].Y)/4;
    double z4 = (nodes[elemNodes[i].N[0]].Z+nodes[elemNodes[i].N[1]].Z+
		 nodes[elemNodes[i].N[2]].Z+nodes[elemNodes[i].N[3]].Z)/4;

    // SoA TriNorm layout: [nx0,nx1,nx2,nx3, ny0,ny1,ny2,ny3, nz0,nz1,nz2,nz3, d0,d1,d2,d3]
    // Face j occupies indices [j], [j+4], [j+8], [j+12].
    // This enables sequential vld1q_f64 loads in PhotonTetrahedronIntersection.

    for (int j=0; j<=3; j++){
      if (triangles[elemNodes[i].T[j]].Num_Elem == 1){
	elems[i].AdjElemIdx[j] = -triangles[elemNodes[i].T[j]].BoundaryIdx;
      }else{
	if(triangles[elemNodes[i].T[j]].ElemIdx[0]==i){
	  elems[i].AdjElemIdx[j] = triangles[elemNodes[i].T[j]].ElemIdx[1];
	}else{
	  elems[i].AdjElemIdx[j] = triangles[elemNodes[i].T[j]].ElemIdx[0];
	}
      }

      double a = nx = triangles[elemNodes[i].T[j]].X;
      double b = ny = triangles[elemNodes[i].T[j]].Y;
      double c = nz = triangles[elemNodes[i].T[j]].Z;
      double d      = triangles[elemNodes[i].T[j]].d;

      double t = -(a*x4 + b*y4 + c*z4 + d)/(a*nx + b*ny + c*nz);

      if(t>0){
	elems[i].TriNorm[j   ] = -nx;
	elems[i].TriNorm[j+ 4] = -ny;
	elems[i].TriNorm[j+ 8] = -nz;
	elems[i].TriNorm[j+12] = -d;
	elemNodes[i].Vol += triNodes[elemNodes[i].T[j]].Area * t / 3.0;
      }else if(t<0){
	elems[i].TriNorm[j   ] = nx;
	elems[i].TriNorm[j+ 4] = ny;
	elems[i].TriNorm[j+ 8] = nz;
	elems[i].TriNorm[j+12] = d;
	elemNodes[i].Vol += -triNodes[elemNodes[i].T[j]].Area * t / 3.0;
      }
    }
  }

  // end of function
  cerr << "End Mesh prepare" << endl;

  delete [] triangleList;

  return true;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// read optical parameters from file
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int ReadOpticalParameter(char       * filename,
			 bool       & simpleOptic,
			 int        & numMed,
			 int        & uniformBoundary,
			 double     & envRefIdx,
			 TMedOptic *& medOptic){
  ifstream fin;
  fin.open(filename, ifstream::in);
  if(!fin.good()){
    cerr << "\tCould not open file: " << filename << endl;
    return -1;
  }

  int type;

  skip_comments(fin);
  fin >> type;
  if (type == 1){
    cerr << "\tOptical parameter per region " << endl;
    simpleOptic = true;
  }else if(type == 2){
    cerr << "\tOptical parameter per element" << endl;
    simpleOptic = false;
  }

  skip_comments(fin);
  fin >> numMed;

  medOptic = new TMedOptic [numMed+1];

  for (int i=1; i<=numMed; i++){
    skip_comments(fin);
    fin >> medOptic[i].mua;
    fin >> medOptic[i].mus;
    fin >> medOptic[i].g;
    fin >> medOptic[i].RefIdx;

    if(medOptic[i].mua<0 || medOptic[i].mus<0){
      cerr << "\tAbsorption and scattering coefficients should be greater than or equal to zero." << endl;
      return -1;
    }

    if(medOptic[i].g<0 ||  medOptic[i].g>1){
      cerr << "\tThe value of g should within [0, 1]." << endl;
      return -1;
    }

    if(medOptic[i].RefIdx<1){
      cerr << "\tThe value of refractive index should >= 1." << endl;
      return -1;
    }

    if(medOptic[i].mua<1e-10 && medOptic[i].mus<1e-10){
      cerr << "\tAbsorption and scattering coefficients could not equal to zero both. " << endl;
      cerr << "\tTo simulate pure glass, we can set g=1, and mua=0.1 to 0.01 and mus=1 to 10" << endl;
      cerr << "\tTIM-OS will not reduce the weight of a photon if it detect g == 1" << endl;
      return -1;
    }
    
    medOptic[i].OneMinsGG = 1 - medOptic[i].g * medOptic[i].g;  // 1-g*g
    medOptic[i].OneMinsG  = 1 - medOptic[i].g;                  // 1-g
    medOptic[i].OneAddGG  = 1 + medOptic[i].g * medOptic[i].g;  // 1+g*g
    medOptic[i].OneAddG   = 1 + medOptic[i].g;                  // 1+g
    medOptic[i].TwoG      = 2 * medOptic[i].g;                  // 2*g
    medOptic[i].pdwa      = medOptic[i].mua/(medOptic[i].mua + medOptic[i].mus); // mua/(mua+mus)
    medOptic[i].MUAMUS    = (medOptic[i].mua + medOptic[i].mus);
    medOptic[i].IMUAMUS   = 1.0/medOptic[i].MUAMUS;
  }  
  
  // read the last optical parameter: refractive index of environment.
  skip_comments(fin);
  fin >> type;
  if(type==1){
    cerr << "\tThe environment has a uniform refractive index." << endl;
    uniformBoundary = 0;
    skip_comments(fin);
    fin >> envRefIdx;
  }else if(type==2){
    cerr << "\tThe environment matches the phantom, hence, no reflection on boundary." << endl; 
    uniformBoundary = 1;
  }
  fin.close();
  return 0;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Prepare the sources
//     The sources defined in the source filename may not contain enough information.
//     This function help fill the lost information
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int Prepare_Source(int        & numNode,
		   int        & numElem,
		   int        & numTrig,
		   int        & numSource,
		   TSource   *& sources,
		   TElem     *& elems,
		   TNode     *& nodes,
		   TElemNode *& elemNodes,
		   TTriNode  *& triNodes,
		   TTriangle *& triangles){
  cerr << "Begin check source" << endl;
  // 
  int * nodeTriStartIdx = new int [numNode+1];
  for (int i=0; i<=numNode; i++){
    nodeTriStartIdx[i]=0;
  }
  for (int i=1; i<=numTrig; i++){
    if(nodeTriStartIdx[triNodes[i].N[0]]==0){
      nodeTriStartIdx[triNodes[i].N[0]]=i;
    }
  }
  nodeTriStartIdx[numNode]==numTrig;
  for (int i=numNode-1; i>=1; i--){
    if(nodeTriStartIdx[i]==0){
      nodeTriStartIdx[i] = nodeTriStartIdx[i+1];
    }
  }
  // In this version, we will do very little pre-process checking 
  // The only thing we will do is finding the element index the point source in
  double height[4];
  for (int i=0; i<numSource; i++){
    if(sources[i].SourceType==1){
      // isotropic point source
      sources[i].ElemIdx = -1;
      for (int Cur_Elem=1; Cur_Elem<=numElem; Cur_Elem++){
	double Min_H = 1000;
	// SoA: face j = [j], [j+4], [j+8], [j+12]
	height[0] = (elems[Cur_Elem].TriNorm[0]*sources[i].Position.X +
		     elems[Cur_Elem].TriNorm[4]*sources[i].Position.Y +
		     elems[Cur_Elem].TriNorm[8]*sources[i].Position.Z +
		     elems[Cur_Elem].TriNorm[12]);
	if(height[0]<0){goto L_ProcSource_Outside;}
	if(height[0]<Min_H){Min_H = height[0];}
	height[1] = (elems[Cur_Elem].TriNorm[1]*sources[i].Position.X +
		     elems[Cur_Elem].TriNorm[5]*sources[i].Position.Y +
		     elems[Cur_Elem].TriNorm[9]*sources[i].Position.Z +
		     elems[Cur_Elem].TriNorm[13]);
	if(height[1]<0){goto L_ProcSource_Outside;}
	if(height[1]<Min_H){Min_H = height[1];}
	height[2] = (elems[Cur_Elem].TriNorm[2]*sources[i].Position.X +
		     elems[Cur_Elem].TriNorm[6]*sources[i].Position.Y +
		     elems[Cur_Elem].TriNorm[10]*sources[i].Position.Z +
		     elems[Cur_Elem].TriNorm[14]);
	if(height[2]<0){goto L_ProcSource_Outside;}
	if(height[2]<Min_H){Min_H = height[2];}
	height[3] = (elems[Cur_Elem].TriNorm[3]*sources[i].Position.X +
		     elems[Cur_Elem].TriNorm[7]*sources[i].Position.Y +
		     elems[Cur_Elem].TriNorm[11]*sources[i].Position.Z +
		     elems[Cur_Elem].TriNorm[15]);
	if(height[3]<0){goto L_ProcSource_Outside;}
	if(height[3]<Min_H){Min_H = height[3];}
	sources[i].ElemIdx = Cur_Elem;
	
	if(Min_H<1e-6){
	  cerr << "\tSource: " << i+1 << " is too close to the boundary of a tetrahedra: " << Cur_Elem << endl;
	  cerr << "\tTIMOS will slightly move its location" << endl;
	  double CX = (nodes[elemNodes[Cur_Elem].N[0]].X +
		       nodes[elemNodes[Cur_Elem].N[1]].X +
		       nodes[elemNodes[Cur_Elem].N[2]].X +
		       nodes[elemNodes[Cur_Elem].N[3]].X)/4.0;
	  double CY = (nodes[elemNodes[Cur_Elem].N[0]].Y +
		       nodes[elemNodes[Cur_Elem].N[1]].Y +
		       nodes[elemNodes[Cur_Elem].N[2]].Y +
		       nodes[elemNodes[Cur_Elem].N[3]].Y)/4.0;
	  double CZ = (nodes[elemNodes[Cur_Elem].N[0]].Z +
		       nodes[elemNodes[Cur_Elem].N[1]].Z +
		       nodes[elemNodes[Cur_Elem].N[2]].Z +
		       nodes[elemNodes[Cur_Elem].N[3]].Z)/4.0;
	  double DX, DY, DZ;
	  DX = CX-sources[i].Position.X;
	  DY = CY-sources[i].Position.Y;
	  DZ = CZ-sources[i].Position.Z;
	  double t = sqrt(DX*DX+DY*DY+DZ*DZ);
	  if(t<1e-6){
	    sources[i].Position.X=CX;
	    sources[i].Position.Y=CY;
	    sources[i].Position.Z=CZ;
	  }else{
	    sources[i].Position.X+=1e-6*DX;
	    sources[i].Position.Y+=1e-6*DY;
	    sources[i].Position.Z+=1e-6*DZ;
	  }
	}
	break;
      L_ProcSource_Outside:
	;
      }
      if(sources[i].ElemIdx<=0){
	cerr << "\tSource " << i+1 << " is not within the phantom" << endl; 
	return -1;
      }
    } else if(sources[i].SourceType==11) {
      // pencil beam source
      
    } else if(sources[i].SourceType==12) {
      // surface triangle region source
      // In this type of source, we assume the light is normal to the surface triangle and point in the object
      // first, sort the nodes in case they are not in order
      for (int j=0; j<=1; j++){
	for (int k=j+1; k<=2; k++){
	  if(sources[i].SurfTriNodes[j]> sources[i].SurfTriNodes[k]){
	    int temp = sources[i].SurfTriNodes[j];
	    sources[i].SurfTriNodes[j]=sources[i].SurfTriNodes[k];
	    sources[i].SurfTriNodes[k]=temp;
	  }
	}
      }
      // Now find it in TriNodes
      int Source_TriIdx=-1;
      for (int j=nodeTriStartIdx[sources[i].SurfTriNodes[0]]; 
	   j<nodeTriStartIdx[sources[i].SurfTriNodes[0]+1]; j++){
	if(sources[i].SurfTriNodes[1]==triNodes[j].N[1] && 
	   sources[i].SurfTriNodes[2]==triNodes[j].N[2]){
	  Source_TriIdx=j;
	  break;
	}
      }
      if(Source_TriIdx==-1){
	cerr << "\tWrong source 1" << endl;
	return -1;
      }
      if(triangles[Source_TriIdx].Num_Elem==2){
	cerr << "\tWrong source 2" << endl;
	return -1;
      }
      sources[i].SurfTriIdx = Source_TriIdx;
      sources[i].ElemIdx    = triangles[Source_TriIdx].ElemIdx[0];
      for (int j=0; j<=3; j++){
	if (sources[i].SurfTriNodes[0]==elemNodes[sources[i].ElemIdx].N[0]){
	  if (sources[i].SurfTriNodes[1]==elemNodes[sources[i].ElemIdx].N[1]){
	    if (sources[i].SurfTriNodes[2]==elemNodes[sources[i].ElemIdx].N[2]){
	      // face 0: SoA indices [0],[4],[8]
	      sources[i].IncAngle.X = elems[sources[i].ElemIdx].TriNorm[0];
	      sources[i].IncAngle.Y = elems[sources[i].ElemIdx].TriNorm[4];
	      sources[i].IncAngle.Z = elems[sources[i].ElemIdx].TriNorm[8];
	    }else{
	      // face 1: SoA indices [1],[5],[9]
	      sources[i].IncAngle.X = elems[sources[i].ElemIdx].TriNorm[1];
	      sources[i].IncAngle.Y = elems[sources[i].ElemIdx].TriNorm[5];
	      sources[i].IncAngle.Z = elems[sources[i].ElemIdx].TriNorm[9];
	    }
	  }else{
	    // face 2: SoA indices [2],[6],[10]
	    sources[i].IncAngle.X = elems[sources[i].ElemIdx].TriNorm[2];
	    sources[i].IncAngle.Y = elems[sources[i].ElemIdx].TriNorm[6];
	    sources[i].IncAngle.Z = elems[sources[i].ElemIdx].TriNorm[10];
	  }
	}else{
	  // face 3: SoA indices [3],[7],[11]
	  sources[i].IncAngle.X = elems[sources[i].ElemIdx].TriNorm[3];
	  sources[i].IncAngle.Y = elems[sources[i].ElemIdx].TriNorm[7];
	  sources[i].IncAngle.Z = elems[sources[i].ElemIdx].TriNorm[11];
	}
      }
      
    }
  }
  cerr << "End check source" << endl;
  delete [] nodeTriStartIdx;
  return 0;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 
// Read source setting from file
// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
long long int ReadSource(char     * filename,
			 int      & numSource,
			 TSource *& sources){

  long long int l_totalPhoton;
  ifstream fin;
  fin.open(filename, ifstream::in);
  if(!fin.good()){
    cerr << "\tCould not open file: filename" << endl;
    return -1;
  }

  int curLine = 0;

  skip_comments(fin);
  fin >> numSource;
  curLine++;
  if(fin.fail()){cerr <<"Error in source file at line " << curLine <<endl;}

  l_totalPhoton  = 0;

  sources = new TSource [numSource];

  for (int i=0; i<numSource; i++){
    curLine++;
    skip_comments(fin);
    fin >> sources[i].SourceType;
    switch(sources[i].SourceType){
    case 1:
      // isotropic point source
      fin >> sources[i].Position.X;
      fin >> sources[i].Position.Y;
      fin >> sources[i].Position.Z;
      if(fin.fail()){cerr <<"Error in source file at line " << curLine <<endl;}
      break;
    case 2:
      // isotropic region source
      fin >> sources[i].ElemIdx;
      break;
    case 11:
      // pencil beam source
      fin >> sources[i].ElemIdx;

      fin >> sources[i].Position.X;
      fin >> sources[i].Position.Y;
      fin >> sources[i].Position.Z;
      
      fin >> sources[i].IncAngle.X;
      fin >> sources[i].IncAngle.Y;
      fin >> sources[i].IncAngle.Z;
      break;
    case 12:
      // surface triangle region source
      // In this type of source, we assume the light is normal to the surface triangle and point in the object
      fin >> sources[i].SurfTriNodes[0];
      fin >> sources[i].SurfTriNodes[1];
      fin >> sources[i].SurfTriNodes[2];
      break;
    default:
      cerr << "\tUnknow source type near line: " << i+1 << endl;
      exit;
      continue;
    }

    skip_comments(fin);
    fin >> sources[i].NumPhoton;

    l_totalPhoton += sources[i].NumPhoton;
  }
  fin.close();
  return l_totalPhoton;
}



inline int Fetch_Source(TSource & ThreadSources, int & SourceIdx, int tid){
  if(SourceIdx >= g_NumSource) return 0;

  int num_fetched_photon = 0;
  os_unfair_lock_lock(&Source_Lock);
  if(SourceIdx < g_NumSource){
    if(g_Sources[SourceIdx].NumPhoton <= 0){
      SourceIdx++;
    }
    if(SourceIdx < g_NumSource){
      ThreadSources.SourceType      = g_Sources[SourceIdx].SourceType;
      ThreadSources.ElemIdx         = g_Sources[SourceIdx].ElemIdx;
      ThreadSources.Position        = g_Sources[SourceIdx].Position;
      ThreadSources.IncAngle        = g_Sources[SourceIdx].IncAngle;
      ThreadSources.SurfTriNodes[0] = g_Sources[SourceIdx].SurfTriNodes[0];
      ThreadSources.SurfTriNodes[1] = g_Sources[SourceIdx].SurfTriNodes[1];
      ThreadSources.SurfTriNodes[2] = g_Sources[SourceIdx].SurfTriNodes[2];
      ThreadSources.SurfTriIdx      = g_Sources[SourceIdx].SurfTriIdx;
      num_fetched_photon = (g_Sources[SourceIdx].NumPhoton >= 1000)
                           ? 1000 : g_Sources[SourceIdx].NumPhoton;
      ThreadSources.NumPhoton = num_fetched_photon;
      g_Sources[SourceIdx].NumPhoton -= num_fetched_photon;
    }
  }
  os_unfair_lock_unlock(&Source_Lock);
  return num_fetched_photon;
}

int InitialPhoton_PointSource(TPhoton & Photon,
			      TSource & ThreadSources,
			      double  * randnums_xyz,
			      int     & Idx_U,
			      timos::Xoshiro256ss & rng){
  // Internal point source
  Photon.Cur_Elem = ThreadSources.ElemIdx;
  Photon.Cur_Med  = g_Elems[Photon.Cur_Elem].MedIdx;
  Photon.Weight   = 1;
  Photon.Path     = 0;
  Photon.X = ThreadSources.Position.X;
  Photon.Y = ThreadSources.Position.Y;
  Photon.Z = ThreadSources.Position.Z;

  double temp;

  do{
    Photon.UX = randnums_xyz[Idx_U++];
    Photon.UY = randnums_xyz[Idx_U++];
    Photon.UZ = randnums_xyz[Idx_U++];
    if(Idx_U>=MAXRANDNUM3){
      Idx_U=0;
      rng_fill_uniform(rng, randnums_xyz, MAXRANDNUM3+6, -1.0, 1.0);
    }
    temp = (Photon.UX * Photon.UX +
	    Photon.UY * Photon.UY +
	    Photon.UZ * Photon.UZ);
  }while(temp>1);

  temp = 1.0/sqrt(temp);
  Photon.UX = Photon.UX*temp;
  Photon.UY = Photon.UY*temp;
  Photon.UZ = Photon.UZ*temp;
  return 0;
}

inline int InitialPhoton_RegionSource(TPhoton & Photon,
				      TSource & ThreadSources,
				      double  * randnums_xyz,
				      int     & Idx_U,
				      timos::Xoshiro256ss & rng){
  // Internal region source
  Photon.Cur_Elem = ThreadSources.ElemIdx;
  Photon.Cur_Med  = g_Elems[Photon.Cur_Elem].MedIdx;
  Photon.Weight   = 1;
  Photon.Path     = 0;
  double wa, wb, wc, wd;
  do{
    wa = (randnums_xyz[Idx_U++]+1)/2.0;
    wb = (randnums_xyz[Idx_U++]+1)/2.0;
    wc = (randnums_xyz[Idx_U++]+1)/2.0;
    if(Idx_U>=MAXRANDNUM3){
      Idx_U=0;
      rng_fill_uniform(rng, randnums_xyz, MAXRANDNUM3+6, -1.0, 1.0);
    }
  }while((wa+wb+wc)>=1);
  wd = 1-wa-wb-wc;
  Photon.X = (g_Nodes[g_ElemNodes[Photon.Cur_Elem].N[0]].X*wa +
	      g_Nodes[g_ElemNodes[Photon.Cur_Elem].N[1]].X*wb +
	      g_Nodes[g_ElemNodes[Photon.Cur_Elem].N[2]].X*wc + 
	      g_Nodes[g_ElemNodes[Photon.Cur_Elem].N[3]].X*wd);
  
  Photon.Y = (g_Nodes[g_ElemNodes[Photon.Cur_Elem].N[0]].Y*wa +
	      g_Nodes[g_ElemNodes[Photon.Cur_Elem].N[1]].Y*wb +
	      g_Nodes[g_ElemNodes[Photon.Cur_Elem].N[2]].Y*wc + 
	      g_Nodes[g_ElemNodes[Photon.Cur_Elem].N[3]].Y*wd);
  
  Photon.Z = (g_Nodes[g_ElemNodes[Photon.Cur_Elem].N[0]].Z*wa +
	      g_Nodes[g_ElemNodes[Photon.Cur_Elem].N[1]].Z*wb +
	      g_Nodes[g_ElemNodes[Photon.Cur_Elem].N[2]].Z*wc +
	      g_Nodes[g_ElemNodes[Photon.Cur_Elem].N[3]].Z*wd);

  double temp;

  do{
    Photon.UX = randnums_xyz[Idx_U++];
    Photon.UY = randnums_xyz[Idx_U++];
    Photon.UZ = randnums_xyz[Idx_U++];
    if(Idx_U>=MAXRANDNUM3){
      Idx_U=0;
      rng_fill_uniform(rng, randnums_xyz, MAXRANDNUM3+6, -1.0, 1.0);
    }
    temp = (Photon.UX * Photon.UX +
	    Photon.UY * Photon.UY +
	    Photon.UZ * Photon.UZ);
  }while(temp>1);

  temp = 1.0/sqrt(temp);
  Photon.UX = Photon.UX*temp;
  Photon.UY = Photon.UY*temp;
  Photon.UZ = Photon.UZ*temp;

  return 0;
}

inline int InitialPhoton_PencilSource(TPhoton & Photon, 
				      TSource & ThreadSources){
  Photon.X = ThreadSources.Position.X;
  Photon.Y = ThreadSources.Position.Y;
  Photon.Z = ThreadSources.Position.Z;
  
  Photon.UX = ThreadSources.IncAngle.X;
  Photon.UY = ThreadSources.IncAngle.Y;
  Photon.UZ = ThreadSources.IncAngle.Z;
  
  double temp = (Photon.UX*Photon.UX + Photon.UY*Photon.UY + Photon.UZ*Photon.UZ);

  Photon.UX = Photon.UX/temp;
  Photon.UY = Photon.UY/temp;
  Photon.UZ = Photon.UZ/temp;
  
  Photon.Cur_Elem = ThreadSources.ElemIdx;
  Photon.Cur_Med  = g_Elems[Photon.Cur_Elem].MedIdx;
  if(g_UniformBoundary==0){
    temp = (g_MedOptic[Photon.Cur_Med].RefIdx-g_EnvRefIdx)/(g_MedOptic[Photon.Cur_Med].RefIdx+g_EnvRefIdx);
    Photon.Weight = 1.0 - temp*temp;
  }else if(g_UniformBoundary==1){
    Photon.Weight = 1.0;
  }
  Photon.Path = 0;
  
  return 0;
}

inline int InitialPhoton_TriangleSource(TPhoton & Photon,
					TSource & ThreadSources,
					double  * randnums_xyz,
					int     & Idx_U,
					timos::Xoshiro256ss & rng){
  double u, v, w;
  u = (randnums_xyz[Idx_U++]+1)/2.0;
  v = (randnums_xyz[Idx_U++]+1)/2.0;
  if(Idx_U>=MAXRANDNUM3){
    Idx_U=0;
    rng_fill_uniform(rng, randnums_xyz, MAXRANDNUM3+6, -1.0, 1.0);
  }
  if(u+v>1){
    u=1-u;
    v=1-v;
  }
  w = 1-u-v;


  //  cerr << "u, v, w: " << u << "\t" << v << "\t" << w << endl;

  Photon.X = (g_Nodes[ThreadSources.SurfTriNodes[0]].X * u + 
	      g_Nodes[ThreadSources.SurfTriNodes[1]].X * v + 
	      g_Nodes[ThreadSources.SurfTriNodes[2]].X * w);
  Photon.Y = (g_Nodes[ThreadSources.SurfTriNodes[0]].Y * u + 
	      g_Nodes[ThreadSources.SurfTriNodes[1]].Y * v + 
	      g_Nodes[ThreadSources.SurfTriNodes[2]].Y * w);
  Photon.Z = (g_Nodes[ThreadSources.SurfTriNodes[0]].Z * u + 
	      g_Nodes[ThreadSources.SurfTriNodes[1]].Z * v + 
	      g_Nodes[ThreadSources.SurfTriNodes[2]].Z * w);
  
  Photon.UX = ThreadSources.IncAngle.X;
  Photon.UY = ThreadSources.IncAngle.Y;
  Photon.UZ = ThreadSources.IncAngle.Z;

  double temp = (Photon.UX*Photon.UX + Photon.UY*Photon.UY + Photon.UZ*Photon.UZ);

  Photon.UX = Photon.UX/temp;
  Photon.UY = Photon.UY/temp;
  Photon.UZ = Photon.UZ/temp;
  
  Photon.Cur_Elem = ThreadSources.ElemIdx;
  Photon.Cur_Med  = g_Elems[Photon.Cur_Elem].MedIdx;
  if(g_UniformBoundary==0){
    temp = (g_MedOptic[Photon.Cur_Med].RefIdx-g_EnvRefIdx)/(g_MedOptic[Photon.Cur_Med].RefIdx+g_EnvRefIdx);
    Photon.Weight = 1.0 - temp*temp;
  }else if(g_UniformBoundary==1){
    Photon.Weight = 1.0;
  }
  Photon.Path   = 0;
  
  return 0;
}

inline void SaveLocal2Global(TPhotonInfo * PhotonInfo, int & PhotonInfo_Idx){
  os_unfair_lock_lock(&Result_Lock);
  if(g_TimeDomain){
    for(int i=1; i<=PhotonInfo_Idx; i++){
      int Type     = PhotonInfo[i].Time_Type & 1;
      int Time_Idx = PhotonInfo[i].Time_Type >> 1;
      if(Type == 0) g_TimeAbsorption[PhotonInfo[i].Idx][Time_Idx] += PhotonInfo[i].LostWeight;
      else          g_TimeSurfMeas  [PhotonInfo[i].Idx][Time_Idx] += PhotonInfo[i].LostWeight;
    }
  }else{
    for(int i=1; i<=PhotonInfo_Idx; i++){
      int Type = PhotonInfo[i].Time_Type & 1;
      if(Type == 0) g_Absorption[PhotonInfo[i].Idx] += PhotonInfo[i].LostWeight;
      else          g_SurfMeas  [PhotonInfo[i].Idx] += PhotonInfo[i].LostWeight;
    }
  }
  PhotonInfo_Idx = 0;
  os_unfair_lock_unlock(&Result_Lock);
}

inline void PhotonTetrahedronIntersection(double  & MinPos,
                                          int     & MinPos_Idx,
                                          double  & MinCos,
                                          double  * Norm,
                                          TPhoton & Photon){
  // Compute step t = -height / dot(N, U) for each of 4 triangle faces.
  // height = dot(N, P) + d   (signed distance from point to plane)
  // t is valid only when dot(N, U) < 0 (photon moving toward face).
  MinPos     = 1e10;
  MinPos_Idx = -1;

#ifdef TIMOS_NEON
  // Compute all 4 dot(N,U) and all 4 height = dot(N,P)+d in parallel.
  // SoA TriNorm: [nx0,nx1,nx2,nx3, ny0,ny1,ny2,ny3, nz0,nz1,nz2,nz3, d0,d1,d2,d3]
  // Sequential vld1q_f64 loads — fully pipelined on M1.
  float64x2_t vUX = vdupq_n_f64(Photon.UX);
  float64x2_t vUY = vdupq_n_f64(Photon.UY);
  float64x2_t vUZ = vdupq_n_f64(Photon.UZ);
  float64x2_t vX  = vdupq_n_f64(Photon.X);
  float64x2_t vY  = vdupq_n_f64(Photon.Y);
  float64x2_t vZ  = vdupq_n_f64(Photon.Z);

  float64x2_t nx01 = vld1q_f64(&Norm[0]);   float64x2_t nx23 = vld1q_f64(&Norm[2]);
  float64x2_t ny01 = vld1q_f64(&Norm[4]);   float64x2_t ny23 = vld1q_f64(&Norm[6]);
  float64x2_t nz01 = vld1q_f64(&Norm[8]);   float64x2_t nz23 = vld1q_f64(&Norm[10]);
  float64x2_t d01  = vld1q_f64(&Norm[12]);  float64x2_t d23  = vld1q_f64(&Norm[14]);

  float64x2_t t01  = vmlaq_f64(vmlaq_f64(vmulq_f64(nx01, vUX), ny01, vUY), nz01, vUZ);
  float64x2_t h01  = vmlaq_f64(vmlaq_f64(vmlaq_f64(d01,  nx01, vX),  ny01, vY),  nz01, vZ);
  float64x2_t t23  = vmlaq_f64(vmlaq_f64(vmulq_f64(nx23, vUX), ny23, vUY), nz23, vUZ);
  float64x2_t h23  = vmlaq_f64(vmlaq_f64(vmlaq_f64(d23,  nx23, vX),  ny23, vY),  nz23, vZ);

  // Extract lanes directly — avoids store-load roundtrip latency (3-5 cycles on M1)
  double t0 = vgetq_lane_f64(t01, 0), h0 = vgetq_lane_f64(h01, 0);
  if(t0 < 0){ double s = -h0/t0; if(s < MinPos){ MinCos = t0; MinPos = s; MinPos_Idx = 0; } }
  double t1 = vgetq_lane_f64(t01, 1), h1 = vgetq_lane_f64(h01, 1);
  if(t1 < 0){ double s = -h1/t1; if(s < MinPos){ MinCos = t1; MinPos = s; MinPos_Idx = 1; } }
  double t2 = vgetq_lane_f64(t23, 0), h2 = vgetq_lane_f64(h23, 0);
  if(t2 < 0){ double s = -h2/t2; if(s < MinPos){ MinCos = t2; MinPos = s; MinPos_Idx = 2; } }
  double t3 = vgetq_lane_f64(t23, 1), h3 = vgetq_lane_f64(h23, 1);
  if(t3 < 0){ double s = -h3/t3; if(s < MinPos){ MinCos = t3; MinPos = s; MinPos_Idx = 3; } }
#else
  double TTemp[4], StepTemp, height;
  TTemp[0] = Norm[0 ]*Photon.UX + Norm[1 ]*Photon.UY + Norm[2 ]*Photon.UZ;
  TTemp[1] = Norm[4 ]*Photon.UX + Norm[5 ]*Photon.UY + Norm[6 ]*Photon.UZ;
  TTemp[2] = Norm[8 ]*Photon.UX + Norm[9 ]*Photon.UY + Norm[10]*Photon.UZ;
  TTemp[3] = Norm[12]*Photon.UX + Norm[13]*Photon.UY + Norm[14]*Photon.UZ;
  for(int i = 0; i <= 3; i++){
    if(TTemp[i] < 0){
      height   = Norm[i*4]*Photon.X + Norm[i*4+1]*Photon.Y + Norm[i*4+2]*Photon.Z + Norm[i*4+3];
      StepTemp = -height / TTemp[i];
      if(StepTemp < MinPos){
        MinCos     = TTemp[i];
        MinPos     = StepTemp;
        MinPos_Idx = i;
      }
    }
  }
#endif
}

inline void ScatterPhoton(TPhoton        & Photon,
                          double         & g,
                          timos::RngPool & rng){

  double cost, sinp, cosp;
  if(g<1){
    double u = rng.get_uniform();
    if(g==0){
      cost = u*2-1;
    }else{
      double tmp = (1.0 - g*g)/(1 - g + 2*g*u);
      cost = (1.0 + g*g - tmp*tmp)/(2*g);
    }
    rng.get_sincos(sinp, cosp);

    double ux = Photon.UX, uy = Photon.UY, uz = Photon.UZ;
    double sint = sqrt(1.0 - cost*cost);

    if(fabs(uz) <= G_COS_0_D){
      double temp1 = sqrt(1.0 - uz*uz);
      double temp  = sint/temp1;
      double temp2 = uz*cosp;
      Photon.UX = (ux*temp2 - uy*sinp)*temp + ux*cost;
      Photon.UY = (uy*temp2 + ux*sinp)*temp + uy*cost;
      Photon.UZ = -sint*cosp*temp1 + uz*cost;
    }else{
      Photon.UX = sint * cosp;
      Photon.UY = sint * sinp;
      Photon.UZ = (uz>0) ? cost : -cost;
    }
  }
}


void * ThreadPhotonPropagation(void * threadid){
  unsigned int tid = (unsigned int)(uintptr_t) threadid;

  // Monte Carlo Simulation Main Function

  // Per-thread accumulators (non-time-domain): array size equals mesh element count —
  // fits entirely in L1 cache. Single lock at thread end; no per-step locking.
  double * LocalAbsorption = new double[g_NumElem + 1]();
  double * LocalSurfMeas   = new double[g_NumBoundaryTrig + 1]();

  // Time-domain still uses the batch PhotonInfo buffer (complex time-bin logic).
  const int     Thread_ArraySize = 1024*1024;
  TPhotonInfo * PhotonInfo = g_TimeDomain ? new TPhotonInfo[Thread_ArraySize+4] : nullptr;
  int           PhotonInfo_Idx = 0;
  if(g_TimeDomain) PhotonInfo[0].Idx = -1;

  // ------------------------------------------------------------
  // declare and initial variables for random number generation
  // rng      — hot-path RngPool: uses Accelerate AES-CTR + vvlog/vvsincos
  // rng_init — cold-path Xoshiro256ss: used only for photon init direction sampling
  // ------------------------------------------------------------
  double * randnums_xyz = new double [MAXRANDNUM3+6];
  int    Idx_U    = 0;
  int    Idx_uvw  = 0;

  uint64_t seed_state = ((uint64_t)(tid + g_StartRandIdx) * 0x9E3779B97F4A7C15ULL)
                        ^ (uint64_t)time(NULL);
  uint64_t seed_a = timos::splitmix64(seed_state);
  uint64_t seed_b = timos::splitmix64(seed_state);

  timos::RngPool      rng(seed_a);
  timos::Xoshiro256ss rng_init(seed_b);

  rng_fill_uniform(rng_init, randnums_xyz, MAXRANDNUM3+6, -1.0, 1.0);

  // ------------------------------------------------------------
  // declare variables for the propagation process
  // ------------------------------------------------------------
  TSource ThreadSources;
  ThreadSources.NumPhoton = 0;

  TPhoton Photon;

  double ux, uy, uz; // temp variable to store UX, UY, UZ
  double NX, NY, NZ; // variables to store the normal vector of a surface
  double Rest_Step;

  double temp, temp1;  

  long long int NumSteps = 0;
  long long int NumIntersections = 0;

  TMedOptic * ThreadMedOptic = new TMedOptic [g_NumMed+1];

  for (int i=0; i<=g_NumMed; i++){
    ThreadMedOptic[i] = g_MedOptic[i];
  }

  // -------------------------------------------------------------
  // Start the main simulation loop
  // -------------------------------------------------------------

  do{
  Label_NewPhoton:
    
    // 1> --------------------------------------------------------------
    // First check ThreadSources.NumPhoton
    // If there is no photon left in NumPhoton
    //    try to fetch some more photons in global Sources
    // EndIf
    // -----------------------------------------------------------------
    if (ThreadSources.NumPhoton <=0){
      int num_fetched_photon = Fetch_Source(ThreadSources, g_SourceIdx, tid); 
      if(num_fetched_photon==0){
	// No more photon left. Exit the main loop and finish the simulation.
	if(tid==0){ cerr << "\r100%" << endl;}
	break;
      }
      g_SimedPhoton += num_fetched_photon;
      if(tid==0){
	cerr << "\r" << (int)((g_SimedPhoton*100)/g_TotalPhoton) << "%";
      }
    }
    ThreadSources.NumPhoton--;
    // ----------------------------------------------------------------

    // ----------------------------------------------------------------
    // Initial a photon
    // ----------------------------------------------------------------
    switch (ThreadSources.SourceType) {
    case 11 : 
      // External pencil beam souce
      InitialPhoton_PencilSource(Photon, ThreadSources);
      break;
    case 12 : 
      // External surface souce
      InitialPhoton_TriangleSource(Photon, ThreadSources, randnums_xyz, Idx_U, rng_init);
      break;
    case 1 :
      // Internal point source
      InitialPhoton_PointSource(Photon, ThreadSources, randnums_xyz, Idx_U, rng_init);
      break;
    case 2 :
      // Internal region source
      InitialPhoton_RegionSource(Photon, ThreadSources, randnums_xyz, Idx_U, rng_init);
      break;
    default : 
      continue;
    }
    // --------- End launch photon -------------------------------------

    // ----------------------------------------------------------------
    // Start photon propagation
    // ----------------------------------------------------------------
    do{
      // --------------------------------------------------------------
      // Set the step size of the photon
      // --------------------------------------------------------------
      NumSteps++;

      Photon.Step = rng.get_neg_log() * ThreadMedOptic[Photon.Cur_Med].IMUAMUS;

      // check whether the photon will move out of the tetrahedron.

      double MinCos, MinPos, Vol;
      int    MinPos_Idx;

    Label_100:

      NumIntersections++;

      PhotonTetrahedronIntersection(MinPos, MinPos_Idx, MinCos, g_Elems[Photon.Cur_Elem].TriNorm, Photon);

      // The reason why we need this is that there is some possibility (~ less than 10^-20) 
      // that a photon already move out of a tetrahedron, but the system still think it is in.
      //      Vol = Vol/3.0;
      //      if(Vol > g_Elems[Photon.Cur_Elem].Vol*1.01){
      //	cerr <<"Error: " <<endl;
      //	NumError++;
      //	ThreadSources.NumPhoton++;
      //	break;
      //      }

      // If the step size of a photon will cross a tetrahedron, 

      if(Photon.Step>MinPos){

	// Precondition
	// IsBTriangle = false
	// 1> Move the photon to the hitted triangle
	//    Rest_Step is the RAW rest step
	//    T0 is the Triangle index
	// 2> If the Triangle is a boundary triangle
	//      Set the tag: IsBTriangle = true
	//      Set OtherRefIdx 
	//    Else
	//      Set OtherRefIdx
	//    End If
	// 3> If OtherRefIdx == CurRefIdx
	//      The photon will directly pass the triangle
	//      If IsBTriangle
	//         photon leave the tissue
	//         save the rest weight to the boundary triangle
	//         break;
	//      Else
	//         photon enter into a new tetrahedron
	//         set new step size
	//         set Loop_back = true
	//         goto Label_100;
	// 4> Else
	//      Photon may be reflected or bended 
	//      a) The first step is calculate the incident angle
	//      b) The second step is calculate the RefPercentage
	//      c) The third step is determine whether the photon 
	//             is reflected or transmitted
	//      d) 
	// 
	//      collect Normal vector of the triangle
	//      compute cosa
	// 
	//      If cosa == G_COS_0_D
	//         RefPercentage = Specular Reflection
	//      Else If cosa = G_COS_90_D
	//         RefPercentage = 0
	//      Else
	//         compute sina
	//         compute sinb
	//         If sinb>1
	//            Total reflection
	//            RefPercentage = 1.0
	//         Else
	//            compute cosb
	//            compute RefPercentage
	//         End
	//      End
	//      
	//      Temp = randnum
	//      If Temp<RefPercentage
	//         Reflected
	//         Compute new Step size
	//         Compute new Direction
	//         Loop_back = true
	//         Goto Label_100
	//         Why need to go back to Label_100 is because that 
	//             the rest step size may be still long enough to move to
	//             another triangle
	//      Else
	//         Transmitted
	//         If IsBTriangle
	//            photon leave the tissue
	//            save the rest weight to the boundary triangle
	//            break;
	//         Else
	//            photon enter to a new tetrahedron
	//            Find Photon.Cur_Elem and Photon.Cur_Med
	//            Compute the new step size
	//            Compute the new Direction
	//            Loop_back = true
	//            Goto Label_100
	//         End
	//      End
	//    End
	//            

	bool   IsBTriangle = false;
	int    OtherMed;
	int    OtherElem;
	double RefPercentage = 0;
	double OtherRefIdx;
	// 1> ---------------------------------------------------------------------------------------
	Photon.X += MinPos*Photon.UX;	
	Photon.Y += MinPos*Photon.UY;	
	Photon.Z += MinPos*Photon.UZ;
	Photon.Path += MinPos*ThreadMedOptic[Photon.Cur_Med].RefIdx;
	Rest_Step = (Photon.Step - MinPos)*ThreadMedOptic[Photon.Cur_Med].MUAMUS;
	// 2> ---------------------------------------------------------------------------------------
	if(g_Elems[Photon.Cur_Elem].AdjElemIdx[MinPos_Idx]>0){
	  // Not a boundary triangle
	  OtherElem = g_Elems[Photon.Cur_Elem].AdjElemIdx[MinPos_Idx];
	  OtherMed  = g_Elems[OtherElem].MedIdx;
	  OtherRefIdx = ThreadMedOptic[OtherMed].RefIdx;
	}else{
	  // Is a boundary triangle
	  IsBTriangle = true;
	  if(g_UniformBoundary==0){
	    OtherRefIdx = g_EnvRefIdx;
	    OtherMed    = 0;
	  }else if(g_UniformBoundary==1){
	    OtherMed = Photon.Cur_Med;
	    OtherRefIdx = ThreadMedOptic[Photon.Cur_Med].RefIdx;
	  }
	}
	// 3> ---------------------------------------------------------------------------------------
	if(OtherMed==Photon.Cur_Med || OtherRefIdx == ThreadMedOptic[Photon.Cur_Med].RefIdx){
	  // no reflection or refraction
	  if(IsBTriangle){
	    if(g_TimeDomain){
	      int Time_Type = (unsigned int)floor(Photon.Path*g_InvLightSpeedMutTimeStep);
	      if(Time_Type >= g_NumTimeStep){goto Label_NewPhoton;}
	      Time_Type = (Time_Type<<1)+1;
	      int Curr_Idx = -g_Elems[Photon.Cur_Elem].AdjElemIdx[MinPos_Idx];
	      if(Curr_Idx == PhotonInfo[PhotonInfo_Idx].Idx && Time_Type == PhotonInfo[PhotonInfo_Idx].Time_Type){
	        PhotonInfo[PhotonInfo_Idx].LostWeight += Photon.Weight;
	      }else{
	        PhotonInfo_Idx++;
	        PhotonInfo[PhotonInfo_Idx].Time_Type = Time_Type;
	        PhotonInfo[PhotonInfo_Idx].Idx = Curr_Idx;
	        PhotonInfo[PhotonInfo_Idx].LostWeight = Photon.Weight;
	      }
	      if(PhotonInfo_Idx >= Thread_ArraySize) SaveLocal2Global(PhotonInfo, PhotonInfo_Idx);
	    }else{
	      LocalSurfMeas[-g_Elems[Photon.Cur_Elem].AdjElemIdx[MinPos_Idx]] += Photon.Weight;
	    }
	    break;
	  }else{
	    Photon.Cur_Elem = OtherElem;
	    Photon.Cur_Med  = OtherMed;
	    // The new step size will be:
	    Photon.Step = Rest_Step * ThreadMedOptic[Photon.Cur_Med].IMUAMUS;
	    goto Label_100;
	  }
	}else{
	  // may have reflection or refraction
	  // 4> ---------------------------------------------------------------------------------------
	  double cosa, sina, cosb, sinb;

	  NX = g_Elems[Photon.Cur_Elem].TriNorm[MinPos_Idx   ];
	  NY = g_Elems[Photon.Cur_Elem].TriNorm[MinPos_Idx+ 4];
	  NZ = g_Elems[Photon.Cur_Elem].TriNorm[MinPos_Idx+ 8];
	  cosa = -MinCos;  //	  cosa = -(UX*NX + UY*NY + UZ*NZ);
	  
	  if(cosa>G_COS_0_D){
	    RefPercentage = ((ThreadMedOptic[Photon.Cur_Med].RefIdx-OtherRefIdx)
			     /(ThreadMedOptic[Photon.Cur_Med].RefIdx+OtherRefIdx));
	    RefPercentage *= RefPercentage;
	    cosb = 1;
	  }else if(cosa<G_COS_90_D){
	    RefPercentage = 1.0;
	  }else{
	    sina = sqrt(1-cosa*cosa);
	    sinb = ThreadMedOptic[Photon.Cur_Med].RefIdx * sina/OtherRefIdx;
	    if(sinb>=1){
	      RefPercentage = 1.0;
	    }else{
	      cosb = sqrt(1.0-sinb*sinb);
	      double cosacosb = cosa*cosb;
	      double sinasinb = sina*sinb;
	      double sinacosb = sina*cosb;
	      double cosasinb = cosa*sinb;
	      double cap = cosacosb - sinasinb;
	      double cam = cosacosb + sinasinb;
	      double sap = sinacosb + cosasinb;
	      double sam = sinacosb - cosasinb;
	      RefPercentage = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam); 
	    }
	  }
	  //	  cerr << "RefPercentage: " << RefPercentage << endl;
	  temp = rng.get_uniform();
	  
	  if (temp <= RefPercentage){
	    // Reflection
	    // The rest step
	    Photon.Step = Photon.Step - MinPos;
	    // new direction vector
	    // W = 2*cos(a)*N + U
	    Photon.UX = 2*cosa*NX + Photon.UX;
	    Photon.UY = 2*cosa*NY + Photon.UY;
	    Photon.UZ = 2*cosa*NZ + Photon.UZ;
	    goto Label_100;
	  }else{
	    // Transmit
	    if(IsBTriangle){
	      if(g_TimeDomain){
		int Time_Type = (unsigned int)floor(Photon.Path*g_InvLightSpeedMutTimeStep);
		if(Time_Type >= g_NumTimeStep){goto Label_NewPhoton;}
		Time_Type = (Time_Type<<1)+1;
		int Curr_Idx = -g_Elems[Photon.Cur_Elem].AdjElemIdx[MinPos_Idx];
		if(Curr_Idx == PhotonInfo[PhotonInfo_Idx].Idx && Time_Type == PhotonInfo[PhotonInfo_Idx].Time_Type){
		  PhotonInfo[PhotonInfo_Idx].LostWeight += Photon.Weight;
		}else{
		  PhotonInfo_Idx++;
		  PhotonInfo[PhotonInfo_Idx].Time_Type = Time_Type;
		  PhotonInfo[PhotonInfo_Idx].Idx = Curr_Idx;
		  PhotonInfo[PhotonInfo_Idx].LostWeight = Photon.Weight;
		}
		if(PhotonInfo_Idx >= Thread_ArraySize) SaveLocal2Global(PhotonInfo, PhotonInfo_Idx);
	      }else{
		LocalSurfMeas[-g_Elems[Photon.Cur_Elem].AdjElemIdx[MinPos_Idx]] += Photon.Weight;
	      }
	      break;
	    }else{
	      // Let assume the transmit vector is WX, WY, WZ. 
	      // We know the normal vector of the triangle is <NX, NY, NZ>
	      //    As our definition, the normal vector points in the tetrahedron
	      // We know the old photon direction is <UX, UY, UZ>
	      // 
	      // W = -cos(b)*N + (sinb/sina)*(cosa*N+U)
	      // 
	      //    
	      Photon.Step     = Rest_Step * ThreadMedOptic[OtherMed].IMUAMUS;
	      double sbdivsa = ThreadMedOptic[Photon.Cur_Med].RefIdx/OtherRefIdx;
	      // sbdivsa = sinb/sina
	      // since sina maybe a very small value, sinb/sina may introduce large error
	      // Since sinb/sina = n1/n2
	      
	      Photon.UX = -cosb*NX + sbdivsa*(cosa*NX+Photon.UX);
	      Photon.UY = -cosb*NY + sbdivsa*(cosa*NY+Photon.UY);
	      Photon.UZ = -cosb*NZ + sbdivsa*(cosa*NZ+Photon.UZ);

	      Photon.Cur_Elem = OtherElem;
	      Photon.Cur_Med  = OtherMed;

	      goto Label_100;
	    }
	  }
	}
      }     

      // Move photon

      Photon.X += Photon.UX * Photon.Step;
      Photon.Y += Photon.UY * Photon.Step;
      Photon.Z += Photon.UZ * Photon.Step;
      Photon.Path += Photon.Step*ThreadMedOptic[Photon.Cur_Med].RefIdx;

      // reduce the weight of the photon
      temp = ThreadMedOptic[Photon.Cur_Med].pdwa*Photon.Weight;

      if(g_TimeDomain){
	int Time_Type = (unsigned int)floor(Photon.Path*g_InvLightSpeedMutTimeStep);
	if(Time_Type >= g_NumTimeStep){goto Label_NewPhoton;}
	Time_Type = (Time_Type<<1);
	int Curr_Idx = Photon.Cur_Elem;
	if(Curr_Idx == PhotonInfo[PhotonInfo_Idx].Idx && Time_Type == PhotonInfo[PhotonInfo_Idx].Time_Type){
	  PhotonInfo[PhotonInfo_Idx].LostWeight += temp;
	}else{
	  PhotonInfo_Idx++;
	  PhotonInfo[PhotonInfo_Idx].Time_Type = Time_Type;
	  PhotonInfo[PhotonInfo_Idx].Idx = Curr_Idx;
	  PhotonInfo[PhotonInfo_Idx].LostWeight = temp;
	}
	if(PhotonInfo_Idx >= Thread_ArraySize) SaveLocal2Global(PhotonInfo, PhotonInfo_Idx);
      }else{
	LocalAbsorption[Photon.Cur_Elem] += temp;
      }
#ifdef NOPUREGLASS
      Photon.Weight -= temp;
#else
      if(ThreadMedOptic[Photon.Cur_Med].g<1){
        Photon.Weight -= temp;
      }
#endif

      if(Photon.Weight <= 0.00001){
        if(rng.get_uniform() < 0.1){
          Photon.Weight *= 10;
        }else{
          break;
        }
      }

      // now start scattering the photon
      // For more information on the scattering, please ref to the original MCML paper/manual.

      ScatterPhoton(Photon, ThreadMedOptic[Photon.Cur_Med].g, rng);

    }while(true);
  }while(true);

  // Flush time-domain buffer; for non-time-domain, merge local accumulators under one lock
  if(g_TimeDomain){
    if(PhotonInfo_Idx > 0) SaveLocal2Global(PhotonInfo, PhotonInfo_Idx);
  }else{
    os_unfair_lock_lock(&Result_Lock);
    for(int i = 1; i <= g_NumElem;          i++) g_Absorption[i] += LocalAbsorption[i];
    for(int i = 1; i <= g_NumBoundaryTrig;  i++) g_SurfMeas[i]   += LocalSurfMeas[i];
    os_unfair_lock_unlock(&Result_Lock);
  }
  g_NumIntersections += NumIntersections;
  g_NumSteps         += NumSteps;

  delete[] randnums_xyz;
  delete[] ThreadMedOptic;
  delete[] LocalAbsorption;
  delete[] LocalSurfMeas;
  if(PhotonInfo) delete[] PhotonInfo;

  return (void *) NULL;
}


void AbsorptionToFluence(int            numElem, 
			 int            numBoundaryTrig,
			 long long int  TotalPhoton,
			 double         sufThreshhold,
			 double         intThreshhold,
			 TTriNode     * triNodes,
			 TMedOptic    * medOptic,
			 int          * boundaryTrigs,
			 TElem        * Elems,
			 TElemNode    * ElemNodes,
			 double       * SurfMeas,
			 double       * Absorption){

  double * CAbs = new double [numElem+1];
  for (int i=0; i<=numElem; i++){
    CAbs[i] = 0;
  }
  double TotalAbs = 0;

  double maxIntFluence = 0;
  double maxSufFluence = 0;
  for (int i=1; i<=numBoundaryTrig; i++){
    SurfMeas[i] = SurfMeas[i]/triNodes[boundaryTrigs[i]].Area;
    if(maxSufFluence < SurfMeas[i]){
      maxSufFluence = SurfMeas[i];
    }
  }
  
  for (int i=1; i<=numElem; i++){
    if(medOptic[Elems[i].MedIdx].mua==0){
      Absorption[i] = -1;
    }else{
      // If g == 1, the material is pure glass/air,
      // then the absorption is for fluence quantification only.
#ifdef NOPUREGLASS
      CAbs[i] = Absorption[i];
#else
      if(medOptic[Elems[i].MedIdx].g<1){
	CAbs[i] = Absorption[i];
      }
#endif
      Absorption[i] = (Absorption[i]/ElemNodes[i].Vol)/medOptic[Elems[i].MedIdx].mua;
      if(maxIntFluence < Absorption[i]){
	maxIntFluence = Absorption[i];
      }
    }
  }
  sort(CAbs, CAbs+numElem+1);
  for (int i=0; i<=numElem; i++){
    TotalAbs += CAbs[i];
  }
  cerr <<"Absorbed Fraction: " << TotalAbs/double(TotalPhoton) << endl;

  sufThreshhold = maxSufFluence /1e20;
  intThreshhold = maxIntFluence /1e20;

  delete [] CAbs;
}

void TimeAbsorptionToFluence(int            numElem, 
			     int            numBoundaryTrig,
			     int            NumTimeStep,
			     long long int  TotalPhoton,
			     double         sufThreshhold,
			     double         intThreshhold,
			     TTriNode     * triNodes,
			     TMedOptic    * medOptic,
			     int          * boundaryTrigs,
			     TElem        * Elems,
			     TElemNode    * ElemNodes,
			     double      ** TimeSurfMeas,
			     double      ** TimeAbsorption){

  double * CAbs = new double [numElem+1];
  for (int i=0; i<=numElem; i++){
    CAbs[i] = 0;
  }
  double TotalAbs = 0;

  double maxIntFluence = 0;
  double maxSufFluence = 0;

  for (int i=1; i<=numBoundaryTrig; i++){
    for (int k=0; k<NumTimeStep; k++){
      TimeSurfMeas[i][k] = TimeSurfMeas[i][k]/triNodes[boundaryTrigs[i]].Area;
      if(maxSufFluence < TimeSurfMeas[i][k]){
	maxSufFluence = TimeSurfMeas[i][k];
      }
    }
  }
  
  for (int i=1; i<=numElem; i++){
    if(medOptic[Elems[i].MedIdx].mua==0){
      for (int k=0; k<NumTimeStep; k++){
	TimeAbsorption[i][k] = -1;
      }
    }else{
      for (int k=0; k<NumTimeStep; k++){
	if(medOptic[Elems[i].MedIdx].g<1){
	  CAbs[i]+=TimeAbsorption[i][k];
	}
	TimeAbsorption[i][k] = (TimeAbsorption[i][k]/ElemNodes[i].Vol)/medOptic[Elems[i].MedIdx].mua;
	if(maxIntFluence < TimeAbsorption[i][k]){
	  maxIntFluence = TimeAbsorption[i][k];
	}
      }
    }
  }
  sort(CAbs, CAbs+numElem+1);
  for (int i=0; i<=numElem; i++){
    TotalAbs += CAbs[i];
  }
  cerr <<"Absorbed Fraction: " << TotalAbs/double(TotalPhoton) << endl;
  sufThreshhold = maxSufFluence /1e20;
  intThreshhold = maxIntFluence /1e20;

  delete [] CAbs;
}


bool WriteResultASCII(char * optical_filename,
		      char * fem_filename,
		      char * source_filename,
		      char * output_filename,
		      long long int  TotalPhoton,
		      int            NumThread,
		      int            StartRandIdx,
		      double         sufThreshhold,
		      double         intThreshhold,
		      int            numElem,
		      int            numBoundaryTrig,
		      TTriNode     * triNodes,
		      int          * boundaryTrigs,
		      TElem        * Elems,
		      TElemNode    * ElemNodes,
		      double       * SurfMeas,
		      double       * Absorption,
		      int            output_format){
  
  ofstream fout;
  
  fout.open(output_filename, ofstream::out);
  if(!fout.good()){
    cerr << "Could not write to output_filename." << endl;
    cerr << "Please check the write permission in this folder" << endl;
    output_filename = "/tmp/timos_tmp_result.dat";
    fout.open(output_filename, ofstream::out);
    if(fout.good()){
      cerr << "Current result will be write to: /tmp/timos_tmp_result.dat" << endl;
    }else{
      return false;
    }
  }

 
  fout <<   "% Optical filename:    " << optical_filename << endl;
  fout <<   "% Fem mesh filename:   " << fem_filename << endl;
  fout <<   "% Source filename:     " << source_filename << endl;
  fout <<   "%      Num Photon:     " << TotalPhoton << endl;
  fout <<   "% Num of threads:      " << NumThread << endl;
  fout <<   "% Start random stream: " << StartRandIdx << endl;
  
  int NumTimeStep = 1;

  if(output_format==1 || output_format==3){
    // output the surface fluence
    fout << "1 " << numBoundaryTrig << " " << NumTimeStep << "\n";

    for (int i=1; i<=numBoundaryTrig; i++){
      fout << triNodes[boundaryTrigs[i]].N[0] << " \t" 
	   << triNodes[boundaryTrigs[i]].N[1] << " \t"
	   << triNodes[boundaryTrigs[i]].N[2] << " \t"
	   << triNodes[boundaryTrigs[i]].Area << " \t";

      if(SurfMeas[i]<sufThreshhold){
	fout << "0\n";
      }else{
	fout << SurfMeas[i] << endl;
      }
    }
  }
    
  if(output_format==2 || output_format==3){
    // output the volumetric fluence
    fout << "2 " << numElem << " " << NumTimeStep << "\n";
    for (int i=1; i<=numElem; i++){
      fout << ElemNodes[i].N[0] << " \t"
	   << ElemNodes[i].N[1] << " \t"
	   << ElemNodes[i].N[2] << " \t"
	   << ElemNodes[i].N[3] << " \t"
	   << ElemNodes[i].Vol  << " \t";

      if(Absorption[i]<intThreshhold){
	fout << "0\n";
      }else{
	fout << Absorption[i] << endl;
      }
    }
  }

  fout.close();
  return true;
}

bool TimeWriteResultASCII(char *         optical_filename,
			  char *         fem_filename,
			  char *         source_filename,
			  char *         output_filename,
			  int            NumTimeStep,
			  double         TimeStep,
			  long long int  TotalPhoton,
			  int            NumThread,
			  int            StartRandIdx,
			  double         sufThreshhold,
			  double         intThreshhold,
			  int            numElem,
			  int            numBoundaryTrig,
			  TTriNode     * triNodes,
			  int          * boundaryTrigs,
			  TElem        * Elems,
			  TElemNode    * ElemNodes,
			  double      ** TimeSurfMeas,
			  double      ** TimeAbsorption,
			  int            output_format){
  
  ofstream fout;
  
  fout.open(output_filename, ofstream::out);
  if(!fout.good()){
    cerr << "Could not write to output_filename." << endl;
    cerr << "Please check the write permission in this folder" << endl;
    output_filename = "/tmp/timos_tmp_result.dat";
    fout.open(output_filename, ofstream::out);
    if(fout.good()){
      cerr << "Current result will be write to: /tmp/timos_tmp_result.dat" << endl;
    }else{
      return false;
    }
  }

  fout <<   "% Optical filename:    " << optical_filename << endl;
  fout <<   "% Fem mesh filename:   " << fem_filename << endl;
  fout <<   "% Source filename:     " << source_filename << endl;
  fout <<   "%      Num Photon:     " << TotalPhoton << endl;
  fout <<   "% Time domain setting. Step size: " << TimeStep 
       <<   " ns, Num Step: " << NumTimeStep << endl;
  fout <<   "% Num of threads:      " << NumThread << endl;
  fout <<   "% Start random stream: " << StartRandIdx << endl;

  if(output_format==1 || output_format==3){
    // output the surface fluence
    fout << "1 " << numBoundaryTrig << " " << NumTimeStep << "\n";

    for (int i=1; i<=numBoundaryTrig; i++){
      fout << triNodes[boundaryTrigs[i]].N[0] << " \t" 
	   << triNodes[boundaryTrigs[i]].N[1] << " \t"
	   << triNodes[boundaryTrigs[i]].N[2] << " \t"
	   << triNodes[boundaryTrigs[i]].Area << " \t";

      for (int j=0; j<NumTimeStep; j++){
	if(TimeSurfMeas[i][j]<sufThreshhold){
	  fout << "0 \t";
	}else{
	  fout << TimeSurfMeas[i][j] << " \t";
	}
      }
      fout << endl;
    }
  }
    
  if(output_format==2 || output_format==3){
    // output the volumetric fluence
    fout << "2 " << numElem << " " << NumTimeStep << "\n";
    for (int i=1; i<=numElem; i++){
      fout << ElemNodes[i].N[0] << " \t"
	   << ElemNodes[i].N[1] << " \t"
	   << ElemNodes[i].N[2] << " \t"
	   << ElemNodes[i].N[3] << " \t"
	   << ElemNodes[i].Vol  << " \t";

      for (int j=0; j<NumTimeStep; j++){
	if(TimeAbsorption[i][j]<intThreshhold){
	  fout << "0 \t";
	}else{
	  fout << TimeAbsorption[i][j] << " \t";
	}
      }
      fout << endl;
    }
  }

  fout.close();
  return true;
}

bool WriteResultBinary(){
  return true;
}

void print_usage(void){
  cerr << "" << endl;
  cerr << "TIM-OS version 3.0" << endl;
  cerr << "Usage: " << endl;
  cerr << "timos -p optical_filename -f finite_element_filename -s source_filename -m [s|i|is|si|t] -o output_filename -T time_step_size(in ns) Num_time_step -t num_thread -r start_rand_id -k mask_filename" << endl << endl;
  cerr << "Note that: all the parameters in TIM-OS are in mm or mm^-1 for time domain simulation. " << endl;
  cerr << "The program supports several output format." << endl;
  cerr << "  -m [s|i|is|si|t]" << endl;
  cerr << "  s:     output surface fluence" << endl;
  cerr << "  i:     output internal fluence" << endl;
  cerr << "  si|is: output surface and internal fluences" << endl << endl;
  cerr << "For time domain simulation, the time_step is in nano-second " << endl;
  cerr << "The program supports multi-thread." << endl;
  cerr << "  If the host system supports multi-core, TIM-OS can launch as many threads as the number of cores. It is also possible to launch more than the number of cores. For some system, launch twice number of threads than the cores can reach maximum performance." << endl << endl;
  cerr << "It is possible to split a large problem into some small jobs. Each job can be ran on a different machine. In this case, it is better to use different random stream for each thread across all the machines. The random number generate used in this program can have at most 1024 independent streams. Since each thread has a different stream, the number of machines can be used for a problem are: 1024/(num_thread_each_machine). For example, if each machine will use 16 threads, the first machine will have a start_rand_id = 1 and the second machine will have a start_rand_id = 17." << endl;
  return;
}

int check_filename(char * filename, int MAX_Length){
  struct stat f_stat;
  int         f_stat_ret;

  if(strlen(filename)>=MAX_Length){
    cerr << "The maximum length of a filename allowed in this program is: " <<   MAX_Length << "." << endl;
    return -1;
  }
  // check whether the file is exist or not
  f_stat_ret = stat(filename, &f_stat);
  if(f_stat_ret==0){
    // check whether it is readable or is it a folder.
    if(f_stat.st_mode == S_IFDIR){
      cerr << "Not a regular file: " << filename << endl;
      return -1;
    }
  }else{
    cerr << "No such file: " << filename << endl;
    return -1;
  }
  return 0;
}

int check_argu(char * argu){
  if(argu[0]=='-'){
    switch(argu[1]){
    case 'p':
      return 1;
    case 'f':
      return 2;
    case 's':
      return 3;
    case 'm':
      return 4;
    case 'o':
      return 5;
    case 't':
      return 6;
    case 'r':
      return 7;
    case 'T':
      return 8;
    default:
      return -1;
    }
  }
  return -1;
}

bool parse_argu(int argc, char * argv[], int & MAX_Length,
		char   * optical_filename, 
		char   * fem_filename,
		char   * source_filename,
		char   * output_filename,
		int    & output_format,
		int    & NumThread,
		int    & StartRandIdx,
		bool   & TimeDomain,
		double & TimeStep,
		double & InvTimeStep,
		double & InvLightSpeedMutTimeStep,
		int    & NumTimeStep){
  int  file_status;
  int  argu_idx = 1;

  bool HasOptFile = false;
  bool HasFEMFile = false;
  bool HasSouFile = false;
  bool HasOutFile = false;
  bool HasFormat  = false;
  bool HasThread  = false;
  bool HasStartRd = false;

  do{
    int option = check_argu(argv[argu_idx]);
    switch(option){
    case 1:
      // the next argv is the optical file name;
      // -p optical_filename 
      argu_idx++;
      file_status = check_filename(argv[argu_idx], MAX_Length);;
      if(file_status != 0){ return -1; } 
      strcpy(optical_filename, argv[argu_idx]);
      argu_idx++;
      HasOptFile = true;
      break;
    case 2:
      // the next argv is the finite element file name;
      // -f finite_element_filename
      argu_idx++;
      file_status = check_filename(argv[argu_idx], MAX_Length);;
      if(file_status != 0){ return -1; } 
      strcpy(fem_filename, argv[argu_idx]);
      argu_idx++;
      HasFEMFile = true;
      break;
    case 3:
      // the next argv is the source filename
      // -s source_filename
      argu_idx++;
      file_status = check_filename(argv[argu_idx], MAX_Length);;
      if(file_status != 0){ return -1; } 
      strcpy(source_filename, argv[argu_idx]);
      argu_idx++;
      HasSouFile = true;
      break;
    case 4:
      // the next argv is the output format
      // -m [s|i|is|si|t]
      // s:     surface fluence               (internal id = 1)
      // i:     internal fluence              (internal id = 2)
      // is|si: surface and internal fluence  (internal id = 3)
      // t:     tecplot                       (internal id = 4)
      argu_idx++;
      if(strlen(argv[argu_idx])==1){
	if(argv[argu_idx][0]=='s'){
	  output_format = 1;
	}else if(argv[argu_idx][0]=='i'){
	  output_format = 2;
	}else if(argv[argu_idx][0]=='t'){
	  output_format = 4;
	}else{
	  cerr << "Unsupported output format." << endl;
	  cerr << "TIM-OS will use the default output format (surface fluence only)." << endl;
	  output_format = 1;
	}
      }else{
	if(strlen(argv[argu_idx])==2){
	  if((argv[argu_idx][0]=='s' && argv[argu_idx][1]=='i') || 
	     (argv[argu_idx][0]=='i' && argv[argu_idx][1]=='s')){
	    output_format = 3;
	  }else{
	    cerr << "Unsupported output format." << endl;
	    cerr << "TIM-OS will use the default output format (surface fluence only)." << endl;
	    output_format = 1;
	  }
	}else{
	  cerr << "Unsupported output format." << endl;
	  cerr << "TIM-OS will use the default output format (surface fluence only)." << endl;
	  output_format = 1;
	}
      }
      argu_idx++;
      HasFormat = true;
      break;
    case 5:
      // the next argv is the output filename
      // -o output_filename
      argu_idx++;
      if(strlen(argv[argu_idx])>=FILENAME_MAX){
	cerr << "The maximum length of a filename allowed in this program is: " <<   MAX_Length << "." << endl;
	cerr << "Current optical filename is longer than this." << endl;
	return -1;
      }
      strcpy(output_filename, argv[argu_idx]);
      argu_idx++;
      HasOutFile = true;
      break;
    case 6:
      // the next argv is the number of thread
      // -t num_thread 
      argu_idx++;
      NumThread = atoi(argv[argu_idx]);
      if(NumThread<=0 || NumThread>=1024){
	NumThread = 64;
      }
      argu_idx++;
      HasThread = true;
      break;
    case 7:
      argu_idx++;
      StartRandIdx = atoi(argv[argu_idx]);
      argu_idx++;
      HasStartRd = true;
      break;
    case 8:
      TimeDomain = true;
      argu_idx++;
      TimeStep = atof(argv[argu_idx]);
      argu_idx++;
      NumTimeStep = atoi(argv[argu_idx]);
      argu_idx++;
      if(NumTimeStep<=1){
	cerr << "Num_Temp_Step should > 1 for time domain simulation" << endl;
	return 0;
      }
      InvTimeStep = 1/TimeStep;
      InvLightSpeedMutTimeStep = INV_LIGHT_SPEED*InvTimeStep;
      break;
    default :
      print_usage();
      return false;
    }
  }while(argu_idx<argc);

  if(HasOptFile==false || HasFEMFile==false || HasSouFile == false || HasOutFile == false){
    print_usage();
    return false;
  }
  if(HasThread == false){
    NumThread = 1;
  }
  if(HasFormat == false){
    output_format = 1;
  }
  if(HasStartRd == false){
    StartRandIdx = 1;
  }else{
    if(StartRandIdx<=0 || StartRandIdx+NumThread>1024){
      cerr << "StartRandIdx is out of range." << endl;
      return false;
    }
  }
  return true;
}

double GetCurrentTime (void) {
  struct timeval timevalue;
  gettimeofday(&timevalue,NULL);
  return (double)timevalue.tv_sec + (double)timevalue.tv_usec/1000000;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Main Function
//----------------------------------------------------------------------------
int main(int argc, char *argv[]) 
{

  // -------------------------------------------------------------------------
  // Parse command line parameters
  // -------------------------------------------------------------------------
  cerr << endl;
  cerr << "------------------------------------------------" << endl;
  if(argc<9){ print_usage(); return -1;
  }

  int MAX_Length = 1024;
  if (FILENAME_MAX<1024){ MAX_Length = FILENAME_MAX; }

  char * optical_filename;
  char * fem_filename;
  char * source_filename;
  char * output_filename;

  optical_filename = new char [MAX_Length+1];
  fem_filename     = new char [MAX_Length+1];
  source_filename  = new char [MAX_Length+1];
  output_filename  = new char [MAX_Length+1];

  int output_format;
  // 001 surface
  // 010 internal
  // 011 surface and internal
  // 100 tecplot

  if(!parse_argu(argc, argv, MAX_Length, 
		 optical_filename, fem_filename, source_filename, output_filename,
		 output_format,
		 g_NumThread, g_StartRandIdx,
		 g_TimeDomain, g_TimeStep, g_InvTimeStep, 
		 g_InvLightSpeedMutTimeStep, g_NumTimeStep)){
    return -1;
  }
  // ---------------------------------------------------------------
  // End parse command line argument
  // ---------------------------------------------------------------

  cerr << "Optical  filename: " << optical_filename << endl;
  cerr << "Fem data filename: " << fem_filename << endl;
  cerr << "Source   filename: " << source_filename << endl;

  // Read optical parameter
  cerr << "Begin read optical parameters." << endl;
  if(ReadOpticalParameter(optical_filename, 
			  g_SimpleOptic, 
			  g_NumMed, 
			  g_UniformBoundary, 
			  g_EnvRefIdx, 
			  g_MedOptic)<0){ 
    return -1; 
  }
  cerr << "End read optical parameters." << endl;

  // Read finite element mesh
  cerr << "Begin read finite element file." << endl;
  if(fem_read(fem_filename, g_SimpleOptic, g_NumNode, g_NumElem, g_NumMed, g_Nodes, g_ElemNodes, g_Elems)<0){ 
    return -1; 
  }
  cerr << "End read finite element file." << endl;

  // Pre-processor finite element mesh
  cerr << "Begin Preprocessing the finite element mesh." << endl;
  if(PreProcessor(g_NumElem, g_NumNode, g_NumMed, g_NumBoundaryTrig, g_NumTrig,
		  g_Nodes, g_ElemNodes, g_Elems, g_TriNodes, 
		  g_Triangles, g_BoundaryTrigs)<0){ 
    return -1; 
  }
  cerr <<   "End Preprocessing the finite element mesh." << endl;

  // Read light source setting
  cerr <<   "Begin read source file." << endl;
  g_TotalPhoton = ReadSource(source_filename, g_NumSource, g_Sources);
  cerr <<   "End read source file." << endl;
  if(g_TotalPhoton<=0){
    return -1;
  }

  if(Prepare_Source(g_NumNode, g_NumElem, g_NumTrig, g_NumSource, g_Sources, g_Elems, g_Nodes, g_ElemNodes, g_TriNodes, g_Triangles)<0){
    return -1;
  }

  // Output the simulation setting on terminal.
  cerr <<   "       Num Photon: " << g_TotalPhoton << endl;
  if(g_TimeDomain){
    cerr << "Time domain setting. Step size: " << g_TimeStep << " ns, Num Step: " << g_NumTimeStep << endl;
  }
  cerr <<   "Output   filename: " << output_filename << endl;
  cerr <<   "Output   format:   ";
  if(output_format==1){
    cerr << "surface fluence" << endl;
  }else if(output_format==2){
    cerr << "internal fluence" << endl;
  }else if(output_format==3){
    cerr << "surface and internal fluence" << endl;
  }
  cerr <<   "Number of threads:   " << g_NumThread << endl;
  cerr <<   "Start random stream: " << g_StartRandIdx << endl;
  //  cerr << optical_filename << " " << fem_filename << " " << source_filename << endl;
  // -----------------------------------------------------------------  
  // Initial memory 
  // -----------------------------------------------------------------
  if(!g_TimeDomain){
    g_SurfMeas = new double [g_NumBoundaryTrig+1];
    for (int j=0; j<=g_NumBoundaryTrig; j++){ g_SurfMeas[j] = 0; }
    g_Absorption = new double [g_NumElem+1];
    for (int j=0; j<=g_NumElem; j++){ g_Absorption[j] = 0; }
  }else{
    g_TimeSurfMeas = new double * [g_NumBoundaryTrig+1];
    for (int j=0; j<=g_NumBoundaryTrig; j++){ 
      g_TimeSurfMeas[j] = new double [g_NumTimeStep];
      for (int k=0; k<g_NumTimeStep; k++){
	g_TimeSurfMeas[j][k] = 0; 
      }
    }
    g_TimeAbsorption = new double * [g_NumElem+1];
    for (int j=0; j<=g_NumElem; j++){ 
      g_TimeAbsorption[j] = new double [g_NumTimeStep];
      for (int k=0; k<g_NumTimeStep; k++){
	g_TimeAbsorption[j][k] = 0;
      }
    }
  }

  //  cerr << optical_filename << " " << fem_filename << " " << source_filename << endl;
  double StartTime = GetCurrentTime();

  g_NumIntersections = 0;
  g_NumSteps = 0;
  g_SimedPhoton = 0;

  // ----------------------------------------------------------------------
  // Start Simulation
  // ----------------------------------------------------------------------
  pthread_t thread_id[g_NumThread];
  for(int i = 0; i < g_NumThread; i++)
    pthread_create(&thread_id[i], NULL, ThreadPhotonPropagation, (void *)(uintptr_t)i);
  for(int i = 0; i < g_NumThread; i++)
    pthread_join(thread_id[i], NULL);
  

  // ----------------------------------------------------------------------
  // End simulation and collect result
  // ----------------------------------------------------------------------
  double EndTime = GetCurrentTime();

  cerr << "Total simulation time: " << EndTime-StartTime << " sec"<< endl;
  cerr << "Num of Intersection:   " << g_NumIntersections  << endl;
  cerr << "Num of Step:   "         << g_NumSteps  << endl;

  double sufThreshhold = 0;
  double intThreshhold = 0;

  cerr << "Prepare output data:\n";

  if(g_TimeDomain){
    TimeAbsorptionToFluence(g_NumElem, g_NumBoundaryTrig, g_NumTimeStep, g_TotalPhoton, 
			    sufThreshhold, intThreshhold, g_TriNodes, g_MedOptic, 
			    g_BoundaryTrigs, g_Elems, g_ElemNodes, g_TimeSurfMeas, 
			    g_TimeAbsorption); 
  }else{
    AbsorptionToFluence(g_NumElem, g_NumBoundaryTrig, g_TotalPhoton, sufThreshhold, 
			intThreshhold, g_TriNodes, g_MedOptic, g_BoundaryTrigs, 
			g_Elems, g_ElemNodes, g_SurfMeas, g_Absorption); 
  }

  cerr << "Write output file \n" ;
  // -------------------------------------------------------------
  // Output the simulation result
  // -------------------------------------------------------------

  if(g_TimeDomain){
    TimeWriteResultASCII(optical_filename, fem_filename, source_filename, 
			 output_filename, 
			 g_NumTimeStep, g_TimeStep, g_TotalPhoton, g_NumThread,
			 g_StartRandIdx, sufThreshhold, intThreshhold, g_NumElem,
			 g_NumBoundaryTrig, g_TriNodes, g_BoundaryTrigs, g_Elems, 
			 g_ElemNodes, g_TimeSurfMeas, g_TimeAbsorption,output_format);
  }else{
    WriteResultASCII(optical_filename, fem_filename, source_filename, 
		     output_filename,
		     g_TotalPhoton, g_NumThread, g_StartRandIdx, 
		     sufThreshhold, intThreshhold, g_NumElem, g_NumBoundaryTrig, 
		     g_TriNodes, g_BoundaryTrigs, g_Elems, g_ElemNodes, g_SurfMeas, 
		     g_Absorption, output_format);
  }

  cerr << "Done" << endl;
  cerr << "------------------------------------------------" << endl;

  delete [] g_Sources;
  delete [] g_MedOptic;
  delete [] g_BoundaryTrigs;
  delete [] g_ElemNodes;
  delete [] g_Elems;
  delete [] g_Nodes;
  delete [] g_TriNodes;
  delete [] g_Triangles;

  delete [] optical_filename;
  delete [] fem_filename;
  delete [] source_filename;
  delete [] output_filename;

 if(!g_TimeDomain){
   delete [] g_SurfMeas;
   delete [] g_Absorption;
  }else{
    for (int j=0; j<=g_NumBoundaryTrig; j++){ 
      delete [] g_TimeSurfMeas[j];
    }
    delete [] g_TimeSurfMeas;

    for (int j=0; j<=g_NumElem; j++){ 
      delete [] g_TimeAbsorption[j];
    }
    delete [] g_TimeAbsorption; 
  }

  cerr << endl;
  return(0);
}
