// Unit tests for mesh I/O and preprocessing (fem_read + PreProcessor).
#include <cassert>
#include <cmath>
#include <iostream>
#include "simcontext.h"
#include "mesh.h"

static void test_fem_read_onelayer(){
  SimContext ctx;
  ctx.simpleOptic = true;
  ctx.numMed = 1; // will be set by optics, but fem_read only checks against it

  auto r = fem_read("example/onelayer/one_layer_18_18_1_2.mesh", ctx);
  assert(r && "fem_read should succeed");
  assert(ctx.numNode == 32 && "one_layer mesh has 32 nodes");
  assert(ctx.numElem == 54 && "one_layer mesh has 54 elements");

  // All nodes should have non-trivially initialised coordinates
  bool any_nonzero = false;
  for(int i = 1; i <= ctx.numNode; i++)
    if(ctx.nodes[i].X != 0.0 || ctx.nodes[i].Y != 0.0 || ctx.nodes[i].Z != 0.0)
      { any_nonzero = true; break; }
  assert(any_nonzero && "at least one node should have non-zero coordinates");

  // All element medium indices should be valid (1-indexed, within numMed)
  for(int i = 1; i <= ctx.numElem; i++)
    assert(ctx.elems[i].MedIdx == 1 && "one-layer mesh: all elements in medium 1");

  std::cerr << "PASS test_fem_read_onelayer\n";
}

static void test_preprocessor_onelayer(){
  SimContext ctx;
  ctx.simpleOptic = true;
  ctx.numMed = 1;
  auto r1 = fem_read("example/onelayer/one_layer_18_18_1_2.mesh", ctx);
  assert(r1);

  auto r2 = PreProcessor(ctx);
  assert(r2 && "PreProcessor should succeed");

  assert(ctx.numTrig > 0 && "should have triangles");
  assert(ctx.numBoundaryTrig > 0 && "should have boundary triangles");
  assert(ctx.numBoundaryTrig < ctx.numTrig && "not all triangles can be boundary");

  // All element volumes should be positive
  double totalVol = 0.0;
  for(int i = 1; i <= ctx.numElem; i++){
    assert(ctx.elemNodes[i].Vol > 0.0 && "element volume must be positive");
    totalVol += ctx.elemNodes[i].Vol;
  }
  // one_layer mesh — measured total volume ≈ 388.8 mm³ (within 1%)
  assert(std::fabs(totalVol - 388.8) / 388.8 < 0.01
         && "total mesh volume should be ~388.8 mm³");

  // SoA TriNorm: each element's inward normals should be unit vectors
  for(int i = 1; i <= ctx.numElem; i++){
    double* N = ctx.elems[i].TriNorm;
    for(int j = 0; j < 4; j++){
      double nx=N[j], ny=N[j+4], nz=N[j+8];
      double len = nx*nx+ny*ny+nz*nz;
      assert(std::fabs(len - 1.0) < 1e-10 && "face normal should be unit vector");
    }
  }

  std::cerr << "PASS test_preprocessor_onelayer  totalVol=" << totalVol << " mm³\n";
}

static void test_fem_read_bad_file(){
  SimContext ctx;
  ctx.simpleOptic = true;
  ctx.numMed = 1;
  auto r = fem_read("/nonexistent/path/file.mesh", ctx);
  assert(!r && "should fail on missing file");
  assert(!r.error().empty() && "error message should not be empty");
  std::cerr << "PASS test_fem_read_bad_file\n";
}

int main(){
  std::cerr << "--- test_mesh ---\n";
  test_fem_read_onelayer();
  test_preprocessor_onelayer();
  test_fem_read_bad_file();
  std::cerr << "All mesh tests passed.\n";
  return 0;
}
