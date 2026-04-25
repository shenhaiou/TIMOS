#define _REENTRANT
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#  include <arm_neon.h>
#  define TIMOS_NEON 1
#endif

#include "simulation.h"
#include "constants.h"
#include "rng.h"

#include <iostream>
#include <cmath>
#include <ctime>
#include <os/lock.h>

using namespace std;

static os_unfair_lock Result_Lock = OS_UNFAIR_LOCK_INIT;
static os_unfair_lock Source_Lock = OS_UNFAIR_LOCK_INIT;

static inline void rng_fill_uniform(timos::Xoshiro256ss& rng, double* buf, int n,
                                     double lo, double hi){
  const double scale = (hi-lo)*0x1.0p-53;
  for(int i=0;i<n;i++) buf[i]=(rng()>>11)*scale+lo;
}

// ----------------------------------------------------------------- Source fetch
static inline int Fetch_Source(TSource& dst, SimContext& ctx, int){
  if(ctx.sourceIdx>=ctx.numSource) return 0;
  int n=0;
  os_unfair_lock_lock(&Source_Lock);
  if(ctx.sourceIdx<ctx.numSource){
    if(ctx.sources[ctx.sourceIdx].NumPhoton<=0) ctx.sourceIdx++;
    if(ctx.sourceIdx<ctx.numSource){
      dst=ctx.sources[ctx.sourceIdx];
      n=(ctx.sources[ctx.sourceIdx].NumPhoton>=SOURCE_BATCH)?SOURCE_BATCH:ctx.sources[ctx.sourceIdx].NumPhoton;
      dst.NumPhoton=n;
      ctx.sources[ctx.sourceIdx].NumPhoton-=n;
    }
  }
  os_unfair_lock_unlock(&Source_Lock);
  return n;
}

// -----------------------------------------------------------------
// Photon initializers
// -----------------------------------------------------------------

// Sample a uniformly random direction on the unit sphere using rejection
// sampling inside the unit cube:  draw (x,y,z) in [-1,1]^3, reject if
//   x^2 + y^2 + z^2 > 1,
// then normalise to obtain a point on S^2.  This avoids trigonometric
// calls and is unbiased (Marsaglia 1972).
static inline void iso_dir(TPhoton& p, double* buf, int& idx, timos::Xoshiro256ss& rng){
  double t;
  do{
    p.UX=buf[idx++]; p.UY=buf[idx++]; p.UZ=buf[idx++];
    if(idx>=RNG_BUF_XYZ3){ idx=0; rng_fill_uniform(rng,buf,RNG_BUF_XYZ3+6,-1.0,1.0); }
    t=p.UX*p.UX+p.UY*p.UY+p.UZ*p.UZ;
  }while(t>1.0);
  double inv=1.0/sqrt(t); p.UX*=inv; p.UY*=inv; p.UZ*=inv;
}

static void InitPoint(TPhoton& p, TSource& s, double* buf, int& idx,
                      timos::Xoshiro256ss& rng, SimContext& ctx){
  p.Cur_Elem=s.ElemIdx; p.Cur_Med=ctx.elems[p.Cur_Elem].MedIdx;
  p.Weight=1.0; p.Path=0.0;
  p.X=s.Position.X; p.Y=s.Position.Y; p.Z=s.Position.Z;
  iso_dir(p,buf,idx,rng);
}

static inline void InitRegion(TPhoton& p, TSource& s, double* buf, int& idx,
                               timos::Xoshiro256ss& rng, SimContext& ctx){
  p.Cur_Elem=s.ElemIdx; p.Cur_Med=ctx.elems[p.Cur_Elem].MedIdx;
  p.Weight=1.0; p.Path=0.0;
  double wa,wb,wc;
  do{
    wa=(buf[idx++]+1)/2.0; wb=(buf[idx++]+1)/2.0; wc=(buf[idx++]+1)/2.0;
    if(idx>=RNG_BUF_XYZ3){ idx=0; rng_fill_uniform(rng,buf,RNG_BUF_XYZ3+6,-1.0,1.0); }
  }while((wa+wb+wc)>=1.0);
  double wd=1.0-wa-wb-wc;
  int* N=ctx.elemNodes[p.Cur_Elem].N;
  p.X=ctx.nodes[N[0]].X*wa+ctx.nodes[N[1]].X*wb+ctx.nodes[N[2]].X*wc+ctx.nodes[N[3]].X*wd;
  p.Y=ctx.nodes[N[0]].Y*wa+ctx.nodes[N[1]].Y*wb+ctx.nodes[N[2]].Y*wc+ctx.nodes[N[3]].Y*wd;
  p.Z=ctx.nodes[N[0]].Z*wa+ctx.nodes[N[1]].Z*wb+ctx.nodes[N[2]].Z*wc+ctx.nodes[N[3]].Z*wd;
  iso_dir(p,buf,idx,rng);
}

static inline void InitPencil(TPhoton& p, TSource& s, SimContext& ctx){
  p.X=s.Position.X; p.Y=s.Position.Y; p.Z=s.Position.Z;
  p.UX=s.IncAngle.X; p.UY=s.IncAngle.Y; p.UZ=s.IncAngle.Z;
  double l2=p.UX*p.UX+p.UY*p.UY+p.UZ*p.UZ; p.UX/=l2; p.UY/=l2; p.UZ/=l2;
  p.Cur_Elem=s.ElemIdx; p.Cur_Med=ctx.elems[p.Cur_Elem].MedIdx;
  double r=(ctx.medOptic[p.Cur_Med].RefIdx-ctx.envRefIdx)/(ctx.medOptic[p.Cur_Med].RefIdx+ctx.envRefIdx);
  p.Weight=(ctx.uniformBoundary==0)?(1.0-r*r):1.0;
  p.Path=0.0;
}

static inline void InitTriangle(TPhoton& p, TSource& s, double* buf, int& idx,
                                 timos::Xoshiro256ss& rng, SimContext& ctx){
  double u=(buf[idx++]+1)/2.0, v=(buf[idx++]+1)/2.0;
  if(idx>=RNG_BUF_XYZ3){ idx=0; rng_fill_uniform(rng,buf,RNG_BUF_XYZ3+6,-1.0,1.0); }
  if(u+v>1.0){ u=1.0-u; v=1.0-v; } double w=1.0-u-v;
  p.X=ctx.nodes[s.SurfTriNodes[0]].X*u+ctx.nodes[s.SurfTriNodes[1]].X*v+ctx.nodes[s.SurfTriNodes[2]].X*w;
  p.Y=ctx.nodes[s.SurfTriNodes[0]].Y*u+ctx.nodes[s.SurfTriNodes[1]].Y*v+ctx.nodes[s.SurfTriNodes[2]].Y*w;
  p.Z=ctx.nodes[s.SurfTriNodes[0]].Z*u+ctx.nodes[s.SurfTriNodes[1]].Z*v+ctx.nodes[s.SurfTriNodes[2]].Z*w;
  p.UX=s.IncAngle.X; p.UY=s.IncAngle.Y; p.UZ=s.IncAngle.Z;
  double l2=p.UX*p.UX+p.UY*p.UY+p.UZ*p.UZ; p.UX/=l2; p.UY/=l2; p.UZ/=l2;
  p.Cur_Elem=s.ElemIdx; p.Cur_Med=ctx.elems[p.Cur_Elem].MedIdx;
  double r=(ctx.medOptic[p.Cur_Med].RefIdx-ctx.envRefIdx)/(ctx.medOptic[p.Cur_Med].RefIdx+ctx.envRefIdx);
  p.Weight=(ctx.uniformBoundary==0)?(1.0-r*r):1.0;
  p.Path=0.0;
}

// ----------------------------------------------------------------- Buffer flush
static inline void SaveLocal2Global(TPhotonInfo* buf, int& idx, SimContext& ctx){
  os_unfair_lock_lock(&Result_Lock);
  if(ctx.timeDomain){
    for(int i=1;i<=idx;i++){
      int type=buf[i].Time_Type&1, t=buf[i].Time_Type>>1;
      if(type==0) ctx.timeAbsorption[buf[i].Idx][t]+=buf[i].LostWeight;
      else        ctx.timeSurfMeas  [buf[i].Idx][t]+=buf[i].LostWeight;
    }
  }else{
    for(int i=1;i<=idx;i++){
      if((buf[i].Time_Type&1)==0) ctx.absorption[buf[i].Idx]+=buf[i].LostWeight;
      else                        ctx.surfMeas  [buf[i].Idx]+=buf[i].LostWeight;
    }
  }
  idx=0;
  os_unfair_lock_unlock(&Result_Lock);
}

// Helper: record a surface-fluence hit (time-domain or CW).
// Returns false if the photon's time index exceeds the window (abandon photon).
static inline bool record_surface(TPhoton& p, int btrig_idx,
                                  SimContext& ctx,
                                  double* LocalSurf,
                                  TPhotonInfo* PInfo, int& PIdx){
  if(ctx.timeDomain){
    int TT=(unsigned)floor(p.Path*ctx.invLightSpeedMutTimeStep);
    if(TT>=ctx.numTimeStep) return false; // abandon
    TT=(TT<<1)+1;
    int Ci=btrig_idx;
    if(Ci==PInfo[PIdx].Idx&&TT==PInfo[PIdx].Time_Type) PInfo[PIdx].LostWeight+=p.Weight;
    else{ PIdx++; PInfo[PIdx]={TT,Ci,p.Weight}; }
    if(PIdx>=TD_FLUSH_SIZE) SaveLocal2Global(PInfo,PIdx,ctx);
  }else{
    LocalSurf[btrig_idx]+=p.Weight;
  }
  return true; // photon exited domain (always abandon after recording)
}

// -----------------------------------------------------------------
// Ray–tetrahedron intersection
// -----------------------------------------------------------------
//
// For a photon at position P = (x,y,z) moving in direction U = (ux,uy,uz),
// find the smallest positive t at which the ray P + t·U crosses a face.
//
// Each face j has an inward unit normal \hat{n}_j = (nx_j, ny_j, nz_j)
// and a plane equation:
//
//   nx_j · x  +  ny_j · y  +  nz_j · z  +  d_j  =  0
//
// The signed distance ("height") from P to face j is:
//
//   h_j  =  \hat{n}_j · P  +  d_j
//
// The ray hits the plane at parameter:
//
//   t_j  =  - h_j / (\hat{n}_j · U)
//
// A crossing is valid only when  \hat{n}_j · U < 0  (photon approaching
// face from inside) and t_j > 0.  We return the face with the smallest
// valid t_j.
//
// TriNorm layout (SoA for NEON efficiency):
//   Norm[0..3]   = { nx_0, nx_1, nx_2, nx_3 }
//   Norm[4..7]   = { ny_0, ny_1, ny_2, ny_3 }
//   Norm[8..11]  = { nz_0, nz_1, nz_2, nz_3 }
//   Norm[12..15] = { d_0,  d_1,  d_2,  d_3  }
//
// On ARM (TIMOS_NEON): two float64x2_t FMA pairs compute all 4 dot
// products simultaneously; vgetq_lane_f64 avoids store-load latency.
static inline void Intersect(double& MinPos, int& MinPos_Idx, double& MinCos,
                              double* Norm, TPhoton& p){
  MinPos=NO_INTERSECTION; MinPos_Idx=-1;
#ifdef TIMOS_NEON
  float64x2_t vUX=vdupq_n_f64(p.UX),vUY=vdupq_n_f64(p.UY),vUZ=vdupq_n_f64(p.UZ);
  float64x2_t vX=vdupq_n_f64(p.X),vY=vdupq_n_f64(p.Y),vZ=vdupq_n_f64(p.Z);
  float64x2_t nx01=vld1q_f64(&Norm[0]),nx23=vld1q_f64(&Norm[2]);
  float64x2_t ny01=vld1q_f64(&Norm[4]),ny23=vld1q_f64(&Norm[6]);
  float64x2_t nz01=vld1q_f64(&Norm[8]),nz23=vld1q_f64(&Norm[10]);
  float64x2_t d01=vld1q_f64(&Norm[12]),d23=vld1q_f64(&Norm[14]);
  float64x2_t t01=vmlaq_f64(vmlaq_f64(vmulq_f64(nx01,vUX),ny01,vUY),nz01,vUZ);
  float64x2_t h01=vmlaq_f64(vmlaq_f64(vmlaq_f64(d01,nx01,vX),ny01,vY),nz01,vZ);
  float64x2_t t23=vmlaq_f64(vmlaq_f64(vmulq_f64(nx23,vUX),ny23,vUY),nz23,vUZ);
  float64x2_t h23=vmlaq_f64(vmlaq_f64(vmlaq_f64(d23,nx23,vX),ny23,vY),nz23,vZ);
  double t0=vgetq_lane_f64(t01,0),h0=vgetq_lane_f64(h01,0);if(t0<0){double s=-h0/t0;if(s<MinPos){MinCos=t0;MinPos=s;MinPos_Idx=0;}}
  double t1=vgetq_lane_f64(t01,1),h1=vgetq_lane_f64(h01,1);if(t1<0){double s=-h1/t1;if(s<MinPos){MinCos=t1;MinPos=s;MinPos_Idx=1;}}
  double t2=vgetq_lane_f64(t23,0),h2=vgetq_lane_f64(h23,0);if(t2<0){double s=-h2/t2;if(s<MinPos){MinCos=t2;MinPos=s;MinPos_Idx=2;}}
  double t3=vgetq_lane_f64(t23,1),h3=vgetq_lane_f64(h23,1);if(t3<0){double s=-h3/t3;if(s<MinPos){MinCos=t3;MinPos=s;MinPos_Idx=3;}}
#else
  double TT[4];
  TT[0]=Norm[0]*p.UX+Norm[1]*p.UY+Norm[2]*p.UZ;
  TT[1]=Norm[4]*p.UX+Norm[5]*p.UY+Norm[6]*p.UZ;
  TT[2]=Norm[8]*p.UX+Norm[9]*p.UY+Norm[10]*p.UZ;
  TT[3]=Norm[12]*p.UX+Norm[13]*p.UY+Norm[14]*p.UZ;
  for(int i=0;i<=3;i++) if(TT[i]<0){
    double h=Norm[i*4]*p.X+Norm[i*4+1]*p.Y+Norm[i*4+2]*p.Z+Norm[i*4+3];
    double s=-h/TT[i]; if(s<MinPos){MinCos=TT[i];MinPos=s;MinPos_Idx=i;}
  }
#endif
}

// -----------------------------------------------------------------
// Henyey-Greenstein scattering
// -----------------------------------------------------------------
//
// The Henyey-Greenstein phase function (HG) for anisotropy factor g in [0,1):
//
//   p(\cos\theta)  =  \frac{1 - g^2}
//                          {4\pi [1 + g^2 - 2g\cos\theta]^{3/2}}
//
// Sampling the deflection angle \theta via inversion (Wang et al. 1995):
//
//   \cos\theta  =  \frac{1}{2g}
//                  \left[ 1 + g^2 - \left(\frac{1-g^2}{1-g+2g\xi}\right)^2 \right]
//
// where \xi ~ Uniform(0,1).  For g=0 (isotropic) the formula reduces to
//   \cos\theta = 2\xi - 1.
//
// The azimuth \phi ~ Uniform(0, 2\pi) is sampled directly.
//
// Rotating the direction (ux,uy,uz) by (\theta,\phi):
//   If |uz| \approx 1 (photon nearly along z-axis) use the simple form:
//     U' = (\sin\theta\cos\phi,  \sin\theta\sin\phi,  \pm\cos\theta)
//   Otherwise apply the general rotation (MCML, Prahl et al. 1989):
//     U'_x = \frac{\sin\theta}{\sqrt{1-uz^2}}
//             (ux\cdot uz\cdot\cos\phi - uy\cdot\sin\phi) + ux\cos\theta
//     U'_y = \frac{\sin\theta}{\sqrt{1-uz^2}}
//             (uy\cdot uz\cdot\cos\phi + ux\cdot\sin\phi) + uy\cos\theta
//     U'_z = -\sin\theta\cos\phi\sqrt{1-uz^2} + uz\cos\theta
static inline void Scatter(TPhoton& p, double g, timos::RngPool& rng){
  if(g>=1.0) return;   // g == 1: pure forward scattering / glass — no deflection
  double u=rng.get_uniform(), cost;
  if(g==0.0) cost=u*2.0-1.0;   // isotropic: \cos\theta = 2\xi - 1
  else{ double tmp=(1.0-g*g)/(1.0-g+2.0*g*u); cost=(1.0+g*g-tmp*tmp)/(2.0*g); }
  double sinp,cosp; rng.get_sincos(sinp,cosp);
  double ux=p.UX,uy=p.UY,uz=p.UZ, sint=sqrt(1.0-cost*cost);
  if(fabs(uz)<=G_COS_0_D){
    // General rotation formula (photon not along z-axis)
    double t1=sqrt(1.0-uz*uz),t=sint/t1,t2=uz*cosp;
    p.UX=(ux*t2-uy*sinp)*t+ux*cost; p.UY=(uy*t2+ux*sinp)*t+uy*cost; p.UZ=-sint*cosp*t1+uz*cost;
  }else{
    // Degenerate case: uz \approx \pm 1, use simplified expression
    p.UX=sint*cosp; p.UY=sint*sinp; p.UZ=(uz>0)?cost:-cost;
  }
}

// -----------------------------------------------------------------
// Propagate one complete photon lifetime
// -----------------------------------------------------------------
//
// The simulation follows a photon through the mesh using the method of
// Wang, Jacques & Zheng (1995, Comp. Meth. Prog. Biomed. 47, 131-146).
//
// Step-size sampling (Beer-Lambert):
//   The free-path length s between scattering/absorption events obeys
//   an exponential distribution with rate \mu_t = \mu_a + \mu_s:
//
//     s = -\ln(\xi) / \mu_t     \xi ~ Uniform(0,1)
//
//   Since \xi is uniform, -\ln(\xi) ~ Exponential(1), so
//   the RNG pool provides pre-computed values of -\ln(\xi) via vvlog.
//
// Weight reduction (implicit absorption, MCML):
//   Rather than terminating at each absorption event, the photon weight W
//   is reduced by the absorbed fraction at each step:
//
//     \Delta W  =  (\mu_a / \mu_t) \cdot W
//     W  \leftarrow  W - \Delta W
//
//   This is the "delta-tracking" / continuous absorption approach.
//
// Russian roulette (Haynsworth, 1956):
//   When  W < W_threshold  the photon has negligible contribution.
//   With probability p_survive = 0.1 it survives with boosted weight
//   W \leftarrow W / p_survive; otherwise it is terminated.
//   This preserves energy expectation:
//     E[\Delta W] = p_survive \cdot (W / p_survive) = W.
//
// Returns PropResult::Abandon if the photon exceeds the time window;
//         PropResult::Done    when the photon exits or is killed.
enum class PropResult { Done, Abandon };

static PropResult propagate_photon(TPhoton& p, TMedOptic* med,
                                   SimContext& ctx,
                                   double* LocalAbs, double* LocalSurf,
                                   TPhotonInfo* PInfo, int& PIdx,
                                   timos::RngPool& rng,
                                   long long& nSteps, long long& nIsect){
  for(;;){ // outer step loop
    nSteps++;
    // s = -ln(\xi) / \mu_t  (IMUAMUS = 1/\mu_t)
    p.Step = rng.get_neg_log() * med[p.Cur_Med].IMUAMUS;

    // Inner boundary-crossing loop (replaces goto Label_100)
    for(;;){
      nIsect++;
      double MinCos, MinPos; int MinPos_Idx;
      Intersect(MinPos, MinPos_Idx, MinCos, ctx.elems[p.Cur_Elem].TriNorm, p);

      if(p.Step <= MinPos) break; // no boundary crossing this iteration

      // Move to hit point
      p.X    += MinPos * p.UX;
      p.Y    += MinPos * p.UY;
      p.Z    += MinPos * p.UZ;
      p.Path += MinPos * med[p.Cur_Med].RefIdx;
      double RestStep = (p.Step - MinPos) * med[p.Cur_Med].MUAMUS;

      // Classify adjacent element
      bool   isBT = false;
      int    otherElem = 0, otherMed = 0;
      double otherN = 0.0;

      if(ctx.elems[p.Cur_Elem].AdjElemIdx[MinPos_Idx] > 0){
        otherElem = ctx.elems[p.Cur_Elem].AdjElemIdx[MinPos_Idx];
        otherMed  = ctx.elems[otherElem].MedIdx;
        otherN    = med[otherMed].RefIdx;
      }else{
        isBT = true;
        if(ctx.uniformBoundary == 0){ otherN = ctx.envRefIdx; otherMed = 0; }
        else{ otherMed = p.Cur_Med; otherN = med[p.Cur_Med].RefIdx; }
      }

      // No refractive-index mismatch
      if(otherMed == p.Cur_Med || otherN == med[p.Cur_Med].RefIdx){
        if(isBT){
          int btIdx = -ctx.elems[p.Cur_Elem].AdjElemIdx[MinPos_Idx];
          if(!record_surface(p, btIdx, ctx, LocalSurf, PInfo, PIdx))
            return PropResult::Abandon;
          return PropResult::Done; // photon exited
        }else{
          p.Cur_Elem = otherElem; p.Cur_Med = otherMed;
          p.Step = RestStep * med[p.Cur_Med].IMUAMUS;
          continue; // re-test intersection in new element
        }
      }

      // Fresnel: compute reflection percentage
      double cosa = -MinCos, cosb = 1.0, refPct = 0.0;
      double NX = ctx.elems[p.Cur_Elem].TriNorm[MinPos_Idx  ];
      double NY = ctx.elems[p.Cur_Elem].TriNorm[MinPos_Idx+4];
      double NZ = ctx.elems[p.Cur_Elem].TriNorm[MinPos_Idx+8];

      if(cosa > G_COS_0_D){
        double r = (med[p.Cur_Med].RefIdx - otherN) / (med[p.Cur_Med].RefIdx + otherN);
        refPct = r * r;
      }else if(cosa < G_COS_90_D){
        refPct = 1.0;
      }else{
        double sina = sqrt(1.0 - cosa*cosa);
        double sinb = med[p.Cur_Med].RefIdx * sina / otherN;
        if(sinb >= 1.0){
          refPct = 1.0;
        }else{
          cosb = sqrt(1.0 - sinb*sinb);
          double cab=cosa*cosb,sab=sina*sinb,scb=sina*cosb,csb=cosa*sinb;
          double cap=cab-sab,cam=cab+sab,sap=scb+csb,sam=scb-csb;
          refPct = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam);
        }
      }

      if(rng.get_uniform() <= refPct){
        // Reflect: re-test from same element with updated direction
        p.Step = p.Step - MinPos;
        p.UX = 2*cosa*NX + p.UX;
        p.UY = 2*cosa*NY + p.UY;
        p.UZ = 2*cosa*NZ + p.UZ;
        continue; // re-test intersection (= old goto Label_100)
      }else{
        // Transmit
        if(isBT){
          int btIdx = -ctx.elems[p.Cur_Elem].AdjElemIdx[MinPos_Idx];
          if(!record_surface(p, btIdx, ctx, LocalSurf, PInfo, PIdx))
            return PropResult::Abandon;
          return PropResult::Done;
        }else{
          double sbd = med[p.Cur_Med].RefIdx / otherN;
          p.Step = RestStep * med[otherMed].IMUAMUS;
          p.UX = -cosb*NX + sbd*(cosa*NX + p.UX);
          p.UY = -cosb*NY + sbd*(cosa*NY + p.UY);
          p.UZ = -cosb*NZ + sbd*(cosa*NZ + p.UZ);
          p.Cur_Elem = otherElem; p.Cur_Med = otherMed;
          continue; // re-test in new element
        }
      }
    } // end inner boundary loop

    // Move photon (no boundary crossing)
    p.X    += p.UX * p.Step;
    p.Y    += p.UY * p.Step;
    p.Z    += p.UZ * p.Step;
    p.Path += p.Step * med[p.Cur_Med].RefIdx;

    // Absorption recording
    double temp = med[p.Cur_Med].pdwa * p.Weight;
    if(ctx.timeDomain){
      int TT = (unsigned)floor(p.Path * ctx.invLightSpeedMutTimeStep);
      if(TT >= ctx.numTimeStep) return PropResult::Abandon;
      TT = (TT<<1);
      int Ci = p.Cur_Elem;
      if(Ci==PInfo[PIdx].Idx && TT==PInfo[PIdx].Time_Type) PInfo[PIdx].LostWeight+=temp;
      else{ PIdx++; PInfo[PIdx]={TT,Ci,temp}; }
      if(PIdx >= TD_FLUSH_SIZE) SaveLocal2Global(PInfo, PIdx, ctx);
    }else{
      LocalAbs[p.Cur_Elem] += temp;
    }

#ifdef NOPUREGLASS
    p.Weight -= temp;
#else
    if(med[p.Cur_Med].g < 1.0) p.Weight -= temp;
#endif

    // Russian roulette
    if(p.Weight <= WEIGHT_THRESHOLD){
      if(rng.get_uniform() < ROULETTE_SURVIVE) p.Weight *= ROULETTE_BOOST;
      else return PropResult::Done;
    }

    Scatter(p, med[p.Cur_Med].g, rng);
  } // end outer step loop
}

// ----------------------------------------------------------------- Thread entry point
void* ThreadPhotonPropagation(void* varg){
  ThreadArg*  arg = static_cast<ThreadArg*>(varg);
  int         tid = arg->tid;
  SimContext& ctx = *arg->ctx;

  double* LocalAbs  = new double[ctx.numElem+1]();
  double* LocalSurf = new double[ctx.numBoundaryTrig+1]();

  TPhotonInfo* PInfo = ctx.timeDomain ? new TPhotonInfo[TD_FLUSH_SIZE+4] : nullptr;
  int          PIdx  = 0;
  if(ctx.timeDomain) PInfo[0].Idx = -1;

  double* xyz   = new double[RNG_BUF_XYZ3+6];
  int     Idx_U = 0;

  uint64_t ss = ((uint64_t)(tid+ctx.startRandIdx)*0x9E3779B97F4A7C15ULL) ^ (uint64_t)time(nullptr);
  uint64_t sa = timos::splitmix64(ss), sb = timos::splitmix64(ss);
  timos::RngPool      rng(sa);
  timos::Xoshiro256ss rng_init(sb);
  rng_fill_uniform(rng_init, xyz, RNG_BUF_XYZ3+6, -1.0, 1.0);

  TSource   ThreadSrc; ThreadSrc.NumPhoton = 0;
  TPhoton   p;
  long long nSteps = 0, nIsect = 0;

  TMedOptic* med = new TMedOptic[ctx.numMed+1];
  for(int i = 0; i <= ctx.numMed; i++) med[i] = ctx.medOptic[i];

  // ---- Photon loop: each iteration is one complete photon lifetime
  for(;;){
    // Fetch a batch of photons if the current batch is exhausted
    if(ThreadSrc.NumPhoton <= 0){
      int n = Fetch_Source(ThreadSrc, ctx, tid);
      if(n == 0){ if(tid==0) cout<<"\r100%\n"; break; }
      ctx.simedPhoton += n;
      if(tid==0) cout<<"\r"<<(int)(ctx.simedPhoton*100/ctx.totalPhoton)<<"%";
    }
    ThreadSrc.NumPhoton--;

    // Initialise photon according to source type
    switch(ThreadSrc.SourceType){
    case 11: InitPencil  (p, ThreadSrc, ctx);               break;
    case 12: InitTriangle(p, ThreadSrc, xyz, Idx_U, rng_init, ctx); break;
    case  1: InitPoint   (p, ThreadSrc, xyz, Idx_U, rng_init, ctx); break;
    case  2: InitRegion  (p, ThreadSrc, xyz, Idx_U, rng_init, ctx); break;
    default: continue; // unknown source type — skip
    }

    // Propagate until photon exits, is killed, or exceeds time window
    propagate_photon(p, med, ctx, LocalAbs, LocalSurf, PInfo, PIdx, rng, nSteps, nIsect);
    // PropResult::Abandon and PropResult::Done are both "start next photon"
  }

  // Merge per-thread accumulators into global arrays
  if(ctx.timeDomain){
    if(PIdx > 0) SaveLocal2Global(PInfo, PIdx, ctx);
  }else{
    os_unfair_lock_lock(&Result_Lock);
    for(int i=1; i<=ctx.numElem;         i++) ctx.absorption[i] += LocalAbs[i];
    for(int i=1; i<=ctx.numBoundaryTrig; i++) ctx.surfMeas[i]   += LocalSurf[i];
    os_unfair_lock_unlock(&Result_Lock);
  }
  ctx.numIntersections += nIsect;
  ctx.numSteps         += nSteps;

  delete[] xyz; delete[] med; delete[] LocalAbs; delete[] LocalSurf;
  if(PInfo) delete[] PInfo;
  return nullptr;
}
