#pragma once
// SimContext holds all simulation state.
// Arrays are std::vector<T> sized n+1 and accessed 1-indexed (index 0 unused).
// RAII: no manual delete[] required anywhere.

#include "timos.h"
#include "constants.h"
#include <expected>
#include <vector>
#include <string>

// Canonical result type for I/O functions.
// Success: TiResult{}   Failure: std::unexpected("message")
using TiResult = std::expected<void, std::string>;

struct SimContext {
  // ------- Time-domain settings -------
  double timeStep                 = 0.1;
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
  double                envRefIdx = 1.0;
  int                   numMed    = 0;
  std::vector<TMedOptic> medOptic;  // size numMed+1, 1-indexed

  // ------- Mesh counts -------
  int numNode         = 0;
  int numElem         = 0;
  int numTrig         = 0;
  int numBoundaryTrig = 0;

  // ------- Mesh arrays: size n+1, accessed 1-indexed (index 0 unused) -------
  std::vector<TNode>     nodes;
  std::vector<TTriNode>  triNodes;
  std::vector<TElemNode> elemNodes;
  std::vector<TElem>     elems;
  std::vector<int>       boundaryTrigs;
  std::vector<TTriangle> triangles;

  // ------- CW result arrays: size n+1, 1-indexed -------
  std::vector<double> surfMeas;
  std::vector<double> absorption;

  // ------- Time-domain result arrays -------
  std::vector<std::vector<double>> timeSurfMeas;
  std::vector<std::vector<double>> timeAbsorption;

  // ------- Sources: size numSource, 0-indexed -------
  int                 numSource = 0;
  std::vector<TSource> sources;
  int                 sourceIdx = 0;

  // ------- Thread count -------
  int numThread = 1;

  // ------- Photon counters -------
  long long int simedPhoton = 0;
  long long int totalPhoton = 0;

  SimContext() = default;
  // Vectors clean up automatically — no destructor needed.
  // Non-copyable (large data, no copy semantics defined).
  SimContext(const SimContext&) = delete;
  SimContext& operator=(const SimContext&) = delete;
};
