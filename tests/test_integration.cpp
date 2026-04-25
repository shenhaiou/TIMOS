// Integration test: run a complete simulation and check absorbed fraction.
// The one-layer slab with mua=0.05, mus=5, g=0.9 should absorb ~17.95% of photons.
#include <cassert>
#include <cmath>
#include <iostream>
#include <thread>
#include <vector>
#include "simcontext.h"
#include "mesh.h"
#include "optics.h"
#include "source_io.h"
#include "simulation.h"
#include "output.h"

int main(){
  std::cerr << "--- test_integration ---\n";

  SimContext ctx;

  // Load inputs
  assert(ReadOpticalParameter("example/onelayer/mua005_mus05.opt", ctx));
  assert(fem_read("example/onelayer/one_layer_18_18_1_2.mesh", ctx));
  assert(PreProcessor(ctx));

  auto src = ReadSource("example/onelayer/one_layer_18_18_1_2.source", ctx);
  assert(src && *src > 0);
  ctx.totalPhoton = *src;
  assert(Prepare_Source(ctx));

  // Allocate result arrays
  ctx.surfMeas  .assign(ctx.numBoundaryTrig+1, 0.0);
  ctx.absorption.assign(ctx.numElem+1,         0.0);
  ctx.numIntersections = 0; ctx.numSteps = 0; ctx.simedPhoton = 0;

  // Run with 4 threads
  ctx.numThread = 4;
  {
    std::vector<ThreadArg> args(ctx.numThread);
    std::vector<std::jthread> threads;
    threads.reserve(ctx.numThread);
    for(int i = 0; i < ctx.numThread; i++){
      args[i] = {i, &ctx};
      threads.emplace_back([&args,i]{ ThreadPhotonPropagation(&args[i]); });
    }
  }

  // Compute absorbed fraction
  double totalAbs = 0.0;
  for(int i = 1; i <= ctx.numElem; i++)
    totalAbs += ctx.absorption[i];
  double absorbedFraction = totalAbs / double(ctx.totalPhoton);

  std::cerr << "Absorbed fraction: " << absorbedFraction
            << "  (expected ~0.1795)\n";

  // Allow ±1% tolerance (Monte Carlo noise)
  assert(std::fabs(absorbedFraction - 0.1795) < 0.01
         && "absorbed fraction should be ~17.95% for mua005_mus05 one-layer");

  std::cerr << "PASS test_integration\n";
  return 0;
}
