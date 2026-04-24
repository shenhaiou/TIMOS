#pragma once
// Extern declarations for all simulation global variables.
// Definitions live in globals.cpp.
// Phase 3 will replace these with a SimContext struct.

#include "timos.h"
#include "constants.h"
#include <vector>

// ------- Time-domain settings -------
extern double g_TimeStep;
extern double g_InvTimeStep;
extern double g_InvLightSpeedMutTimeStep;
extern int    g_NumTimeStep;

// ------- Diagnostic counters -------
extern long long int g_NumIntersections;
extern long long int g_NumSteps;

// ------- Optics mode -------
extern bool g_SimpleOptic;
// true  = params per region; false = params per element (not fully tested)

extern int g_UniformBoundary;
// 0 = uniform env refractive index; 1 = matched boundary (no reflection)

// ------- Simulation mode -------
extern bool g_TimeDomain;

// ------- Random stream -------
extern int g_StartRandIdx;

// ------- Mesh / result format (unused) -------
extern int g_InternalResultFormat;

// ------- Optics -------
extern double     g_EnvRefIdx;
extern int        g_NumMed;
extern TMedOptic* g_MedOptic;

// ------- Mesh counts -------
extern int g_NumNode;
extern int g_NumElem;
extern int g_NumTrig;
extern int g_NumBoundaryTrig;

// ------- Mesh arrays (1-indexed) -------
extern TNode*     g_Nodes;
extern TTriNode*  g_TriNodes;
extern TElemNode* g_ElemNodes;
extern TElem*     g_Elems;
extern int*       g_BoundaryTrigs;
extern TTriangle* g_Triangles;

// ------- CW result arrays (1-indexed) -------
extern double* g_SurfMeas;
extern double* g_Absorption;

// ------- Time-domain result arrays -------
extern std::vector<std::vector<double>> g_TimeSurfMeas;
extern std::vector<std::vector<double>> g_TimeAbsorption;

// ------- Sources -------
extern int      g_NumSource;
extern TSource* g_Sources;
extern int      g_SourceIdx;

// ------- Thread count -------
extern int g_NumThread;

// ------- Photon counters -------
extern long long int g_SimedPhoton;
extern long long int g_TotalPhoton;
