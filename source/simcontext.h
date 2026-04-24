#pragma once
// SimContext holds all simulation state that was previously global (g_* variables).
// Pass by reference through all functions; eliminates global mutable state.
// Phase 4 will migrate raw pointer arrays to std::vector (0-indexed).

#include "timos.h"
#include "constants.h"
#include <expected>
#include <vector>
#include <string>

// Canonical result type for I/O functions: success = std::expected<void,std::string>{}
// Failure = std::unexpected("human-readable message").
using TiResult = std::expected<void, std::string>;

struct SimContext {
  // ------- Time-domain settings -------
  double timeStep                 = 0.1;    // ns
  double invTimeStep              = 10.0;
  double invLightSpeedMutTimeStep = INV_LIGHT_SPEED * 10.0;
  int    numTimeStep              = 50;

  // ------- Diagnostic counters -------
  long long int numIntersections = 0;
  long long int numSteps         = 0;

  // ------- Optics mode -------
  bool simpleOptic     = true;   // true = per-region, false = per-element
  int  uniformBoundary = 0;      // 0 = uniform env n, 1 = matched (no reflection)

  // ------- Simulation mode -------
  bool timeDomain = false;

  // ------- Random stream -------
  int startRandIdx = 1;

  // ------- Result format (future use) -------
  int internalResultFormat = 0;

  // ------- Optics -------
  double     envRefIdx = 1.0;
  int        numMed    = 0;
  TMedOptic* medOptic  = nullptr;

  // ------- Mesh counts -------
  int numNode         = 0;
  int numElem         = 0;
  int numTrig         = 0;
  int numBoundaryTrig = 0;

  // ------- Mesh arrays (1-indexed, raw pointers — Phase 4 will replace with vector) -------
  TNode*     nodes        = nullptr;
  TTriNode*  triNodes     = nullptr;
  TElemNode* elemNodes    = nullptr;
  TElem*     elems        = nullptr;
  int*       boundaryTrigs= nullptr;
  TTriangle* triangles    = nullptr;

  // ------- CW result arrays (1-indexed) -------
  double* surfMeas   = nullptr;
  double* absorption = nullptr;

  // ------- Time-domain result arrays -------
  std::vector<std::vector<double>> timeSurfMeas;
  std::vector<std::vector<double>> timeAbsorption;

  // ------- Sources -------
  int      numSource = 0;
  TSource* sources   = nullptr;
  int      sourceIdx = 0;          // protected by Source_Lock in simulation.cpp

  // ------- Thread count -------
  int numThread = 1;

  // ------- Photon counters (updated under Result_Lock) -------
  long long int simedPhoton = 0;
  long long int totalPhoton = 0;

  // ------- Destructor: free all heap arrays -------
  ~SimContext(){
    delete[] sources;   delete[] medOptic;      delete[] boundaryTrigs;
    delete[] elemNodes; delete[] elems;         delete[] nodes;
    delete[] triNodes;  delete[] triangles;
    delete[] surfMeas;  delete[] absorption;
  }

  // Non-copyable (owns raw pointers)
  SimContext(const SimContext&) = delete;
  SimContext& operator=(const SimContext&) = delete;
  SimContext() = default;
};
