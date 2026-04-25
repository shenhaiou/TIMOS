// Unit tests for optical parameter I/O and precomputed values.
#include <cassert>
#include <cmath>
#include <iostream>
#include "simcontext.h"
#include "optics.h"

static void check_derived(const TMedOptic& m){
  double tol = 1e-12;
  assert(std::fabs(m.MUAMUS - (m.mua + m.mus)) < tol);
  assert(std::fabs(m.IMUAMUS - 1.0/(m.mua+m.mus)) < tol);
  assert(std::fabs(m.pdwa - m.mua/(m.mua+m.mus)) < tol);
  assert(std::fabs(m.OneMinsGG - (1.0 - m.g*m.g)) < tol);
  assert(std::fabs(m.TwoG - 2.0*m.g) < tol);
}

static void test_read_onelayer_opt(){
  SimContext ctx;
  auto r = ReadOpticalParameter("example/onelayer/mua005_mus05.opt", ctx);
  assert(r && "ReadOpticalParameter should succeed");
  assert(ctx.numMed == 1 && "one-layer has 1 medium");
  assert(ctx.simpleOptic && "should be per-region");
  assert(ctx.uniformBoundary == 0 && "uniform env refractive index (boundaryFlag=1)");
  assert(std::fabs(ctx.envRefIdx - 1.0) < 1e-9 && "envRefIdx should be 1.0");

  const TMedOptic& m = ctx.medOptic[1];
  assert(std::fabs(m.mua - 0.05) < 1e-9);
  assert(std::fabs(m.mus - 5.0)  < 1e-9);
  assert(std::fabs(m.g   - 0.9)  < 1e-9);
  assert(std::fabs(m.RefIdx - 1.46) < 1e-9);
  check_derived(m);

  std::cout << "PASS test_read_onelayer_opt  pdwa=" << m.pdwa << "\n";
}

static void test_read_fourlayer_opt(){
  SimContext ctx;
  auto r = ReadOpticalParameter("example/fourlayer/FourLayer.opt", ctx);
  assert(r);
  assert(ctx.numMed == 4);
  for(int i = 1; i <= ctx.numMed; i++) check_derived(ctx.medOptic[i]);
  std::cout << "PASS test_read_fourlayer_opt\n";
}

static void test_read_mouse_opt(){
  SimContext ctx;
  auto r = ReadOpticalParameter("example/mouse/mouse.opt", ctx);
  assert(r);
  assert(ctx.numMed == 17);
  for(int i = 1; i <= ctx.numMed; i++){
    assert(ctx.medOptic[i].mua > 0.0);
    assert(ctx.medOptic[i].mus > 0.0);
    assert(ctx.medOptic[i].g >= 0.0 && ctx.medOptic[i].g <= 1.0);
    check_derived(ctx.medOptic[i]);
  }
  std::cout << "PASS test_read_mouse_opt  numMed=17\n";
}

static void test_read_bad_opt(){
  SimContext ctx;
  auto r = ReadOpticalParameter("/no/such/file.opt", ctx);
  assert(!r && "should fail on missing file");
  std::cout << "PASS test_read_bad_opt\n";
}

int main(){
  std::cout << "--- test_optics ---\n";
  test_read_onelayer_opt();
  test_read_fourlayer_opt();
  test_read_mouse_opt();
  test_read_bad_opt();
  std::cout << "All optics tests passed.\n";
  return 0;
}
