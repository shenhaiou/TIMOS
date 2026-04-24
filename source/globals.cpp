#include "globals.h"

double g_TimeStep                 = 0.1;
double g_InvTimeStep              = 1.0 / 0.1;
double g_InvLightSpeedMutTimeStep = INV_LIGHT_SPEED / 0.1;
int    g_NumTimeStep              = 50;

long long int g_NumIntersections = 0;
long long int g_NumSteps         = 0;

bool g_SimpleOptic     = true;
int  g_UniformBoundary = 0;
bool g_TimeDomain      = false;
int  g_StartRandIdx    = 1;

int g_InternalResultFormat = 0;

double     g_EnvRefIdx = 1.0;
int        g_NumMed    = 0;
TMedOptic* g_MedOptic  = nullptr;

int g_NumNode         = 0;
int g_NumElem         = 0;
int g_NumTrig         = 0;
int g_NumBoundaryTrig = 0;

TNode*     g_Nodes        = nullptr;
TTriNode*  g_TriNodes     = nullptr;
TElemNode* g_ElemNodes    = nullptr;
TElem*     g_Elems        = nullptr;
int*       g_BoundaryTrigs= nullptr;
TTriangle* g_Triangles    = nullptr;

double* g_SurfMeas   = nullptr;
double* g_Absorption = nullptr;

std::vector<std::vector<double>> g_TimeSurfMeas;
std::vector<std::vector<double>> g_TimeAbsorption;

int      g_NumSource = 0;
TSource* g_Sources   = nullptr;
int      g_SourceIdx = 0;

int g_NumThread = 1;

long long int g_SimedPhoton = 0;
long long int g_TotalPhoton = 0;
