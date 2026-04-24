#define _REENTRANT
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#  include <arm_neon.h>
#  define TIMOS_NEON 1
#endif

#include "simulation.h"
#include "globals.h"
#include "constants.h"
#include "rng.h"

#include <iostream>
#include <cmath>
#include <ctime>
#include <os/lock.h>

using namespace std;

// Per-file locks — only accessed from simulation functions
static os_unfair_lock Result_Lock = OS_UNFAIR_LOCK_INIT;
static os_unfair_lock Source_Lock = OS_UNFAIR_LOCK_INIT;

// Fill buf[0..n) with uniform doubles in [lo, hi)
static inline void rng_fill_uniform(timos::Xoshiro256ss& rng, double* buf, int n, double lo, double hi){
  const double scale = (hi - lo) * 0x1.0p-53;
  for(int i = 0; i < n; i++)
    buf[i] = (rng() >> 11) * scale + lo;
}

// +++ Source fetcher +++

static inline int Fetch_Source(TSource& ThreadSources, int& SourceIdx, int /*tid*/){
  if(SourceIdx >= g_NumSource) return 0;
  int n = 0;
  os_unfair_lock_lock(&Source_Lock);
  if(SourceIdx < g_NumSource){
    if(g_Sources[SourceIdx].NumPhoton <= 0) SourceIdx++;
    if(SourceIdx < g_NumSource){
      ThreadSources = g_Sources[SourceIdx];
      n = (g_Sources[SourceIdx].NumPhoton >= SOURCE_BATCH)
          ? SOURCE_BATCH : g_Sources[SourceIdx].NumPhoton;
      ThreadSources.NumPhoton = n;
      g_Sources[SourceIdx].NumPhoton -= n;
    }
  }
  os_unfair_lock_unlock(&Source_Lock);
  return n;
}

// +++ Photon initializers +++

static inline void init_isotropic_direction(TPhoton& Photon, double* buf, int& idx,
                                             timos::Xoshiro256ss& rng){
  double temp;
  do{
    Photon.UX = buf[idx++]; Photon.UY = buf[idx++]; Photon.UZ = buf[idx++];
    if(idx >= RNG_BUF_XYZ3){ idx = 0; rng_fill_uniform(rng, buf, RNG_BUF_XYZ3+6, -1.0, 1.0); }
    temp = Photon.UX*Photon.UX + Photon.UY*Photon.UY + Photon.UZ*Photon.UZ;
  }while(temp > 1.0);
  double inv = 1.0 / sqrt(temp);
  Photon.UX *= inv; Photon.UY *= inv; Photon.UZ *= inv;
}

static int InitialPhoton_PointSource(TPhoton& Photon, TSource& src,
                                     double* buf, int& idx, timos::Xoshiro256ss& rng){
  Photon.Cur_Elem = src.ElemIdx;
  Photon.Cur_Med  = g_Elems[Photon.Cur_Elem].MedIdx;
  Photon.Weight   = 1.0; Photon.Path = 0.0;
  Photon.X = src.Position.X; Photon.Y = src.Position.Y; Photon.Z = src.Position.Z;
  init_isotropic_direction(Photon, buf, idx, rng);
  return 0;
}

static inline int InitialPhoton_RegionSource(TPhoton& Photon, TSource& src,
                                             double* buf, int& idx, timos::Xoshiro256ss& rng){
  Photon.Cur_Elem = src.ElemIdx;
  Photon.Cur_Med  = g_Elems[Photon.Cur_Elem].MedIdx;
  Photon.Weight   = 1.0; Photon.Path = 0.0;
  double wa, wb, wc;
  do{
    wa = (buf[idx++]+1)/2.0; wb = (buf[idx++]+1)/2.0; wc = (buf[idx++]+1)/2.0;
    if(idx >= RNG_BUF_XYZ3){ idx = 0; rng_fill_uniform(rng, buf, RNG_BUF_XYZ3+6, -1.0, 1.0); }
  }while((wa+wb+wc) >= 1.0);
  double wd = 1.0 - wa - wb - wc;
  int* N = g_ElemNodes[Photon.Cur_Elem].N;
  Photon.X = g_Nodes[N[0]].X*wa + g_Nodes[N[1]].X*wb + g_Nodes[N[2]].X*wc + g_Nodes[N[3]].X*wd;
  Photon.Y = g_Nodes[N[0]].Y*wa + g_Nodes[N[1]].Y*wb + g_Nodes[N[2]].Y*wc + g_Nodes[N[3]].Y*wd;
  Photon.Z = g_Nodes[N[0]].Z*wa + g_Nodes[N[1]].Z*wb + g_Nodes[N[2]].Z*wc + g_Nodes[N[3]].Z*wd;
  init_isotropic_direction(Photon, buf, idx, rng);
  return 0;
}

static inline int InitialPhoton_PencilSource(TPhoton& Photon, TSource& src){
  Photon.X = src.Position.X; Photon.Y = src.Position.Y; Photon.Z = src.Position.Z;
  Photon.UX = src.IncAngle.X; Photon.UY = src.IncAngle.Y; Photon.UZ = src.IncAngle.Z;
  double len2 = Photon.UX*Photon.UX + Photon.UY*Photon.UY + Photon.UZ*Photon.UZ;
  Photon.UX /= len2; Photon.UY /= len2; Photon.UZ /= len2;
  Photon.Cur_Elem = src.ElemIdx;
  Photon.Cur_Med  = g_Elems[Photon.Cur_Elem].MedIdx;
  if(g_UniformBoundary == 0){
    double r = (g_MedOptic[Photon.Cur_Med].RefIdx - g_EnvRefIdx)
             / (g_MedOptic[Photon.Cur_Med].RefIdx + g_EnvRefIdx);
    Photon.Weight = 1.0 - r*r;
  }else{
    Photon.Weight = 1.0;
  }
  Photon.Path = 0.0;
  return 0;
}

static inline int InitialPhoton_TriangleSource(TPhoton& Photon, TSource& src,
                                               double* buf, int& idx, timos::Xoshiro256ss& rng){
  double u = (buf[idx++]+1)/2.0, v = (buf[idx++]+1)/2.0;
  if(idx >= RNG_BUF_XYZ3){ idx = 0; rng_fill_uniform(rng, buf, RNG_BUF_XYZ3+6, -1.0, 1.0); }
  if(u+v > 1.0){ u = 1.0-u; v = 1.0-v; }
  double w = 1.0 - u - v;
  Photon.X = g_Nodes[src.SurfTriNodes[0]].X*u + g_Nodes[src.SurfTriNodes[1]].X*v + g_Nodes[src.SurfTriNodes[2]].X*w;
  Photon.Y = g_Nodes[src.SurfTriNodes[0]].Y*u + g_Nodes[src.SurfTriNodes[1]].Y*v + g_Nodes[src.SurfTriNodes[2]].Y*w;
  Photon.Z = g_Nodes[src.SurfTriNodes[0]].Z*u + g_Nodes[src.SurfTriNodes[1]].Z*v + g_Nodes[src.SurfTriNodes[2]].Z*w;
  Photon.UX = src.IncAngle.X; Photon.UY = src.IncAngle.Y; Photon.UZ = src.IncAngle.Z;
  double len2 = Photon.UX*Photon.UX + Photon.UY*Photon.UY + Photon.UZ*Photon.UZ;
  Photon.UX /= len2; Photon.UY /= len2; Photon.UZ /= len2;
  Photon.Cur_Elem = src.ElemIdx;
  Photon.Cur_Med  = g_Elems[Photon.Cur_Elem].MedIdx;
  if(g_UniformBoundary == 0){
    double r = (g_MedOptic[Photon.Cur_Med].RefIdx - g_EnvRefIdx)
             / (g_MedOptic[Photon.Cur_Med].RefIdx + g_EnvRefIdx);
    Photon.Weight = 1.0 - r*r;
  }else{
    Photon.Weight = 1.0;
  }
  Photon.Path = 0.0;
  return 0;
}

// +++ Flush thread-local PhotonInfo buffer to global arrays +++

static inline void SaveLocal2Global(TPhotonInfo* PhotonInfo, int& idx){
  os_unfair_lock_lock(&Result_Lock);
  if(g_TimeDomain){
    for(int i = 1; i <= idx; i++){
      int type = PhotonInfo[i].Time_Type & 1;
      int t    = PhotonInfo[i].Time_Type >> 1;
      if(type == 0) g_TimeAbsorption[PhotonInfo[i].Idx][t] += PhotonInfo[i].LostWeight;
      else          g_TimeSurfMeas  [PhotonInfo[i].Idx][t] += PhotonInfo[i].LostWeight;
    }
  }else{
    for(int i = 1; i <= idx; i++){
      if((PhotonInfo[i].Time_Type & 1) == 0)
        g_Absorption[PhotonInfo[i].Idx] += PhotonInfo[i].LostWeight;
      else
        g_SurfMeas  [PhotonInfo[i].Idx] += PhotonInfo[i].LostWeight;
    }
  }
  idx = 0;
  os_unfair_lock_unlock(&Result_Lock);
}

// +++ Ray–tetrahedron intersection (NEON-accelerated on M1) +++

static inline void PhotonTetrahedronIntersection(double& MinPos, int& MinPos_Idx,
                                                 double& MinCos, double* Norm,
                                                 TPhoton& Photon){
  MinPos     = NO_INTERSECTION;
  MinPos_Idx = -1;

#ifdef TIMOS_NEON
  float64x2_t vUX=vdupq_n_f64(Photon.UX), vUY=vdupq_n_f64(Photon.UY), vUZ=vdupq_n_f64(Photon.UZ);
  float64x2_t vX =vdupq_n_f64(Photon.X),  vY =vdupq_n_f64(Photon.Y),  vZ =vdupq_n_f64(Photon.Z);
  float64x2_t nx01=vld1q_f64(&Norm[0]),  nx23=vld1q_f64(&Norm[2]);
  float64x2_t ny01=vld1q_f64(&Norm[4]),  ny23=vld1q_f64(&Norm[6]);
  float64x2_t nz01=vld1q_f64(&Norm[8]),  nz23=vld1q_f64(&Norm[10]);
  float64x2_t d01 =vld1q_f64(&Norm[12]), d23 =vld1q_f64(&Norm[14]);
  float64x2_t t01=vmlaq_f64(vmlaq_f64(vmulq_f64(nx01,vUX),ny01,vUY),nz01,vUZ);
  float64x2_t h01=vmlaq_f64(vmlaq_f64(vmlaq_f64(d01,nx01,vX),ny01,vY),nz01,vZ);
  float64x2_t t23=vmlaq_f64(vmlaq_f64(vmulq_f64(nx23,vUX),ny23,vUY),nz23,vUZ);
  float64x2_t h23=vmlaq_f64(vmlaq_f64(vmlaq_f64(d23,nx23,vX),ny23,vY),nz23,vZ);
  double t0=vgetq_lane_f64(t01,0),h0=vgetq_lane_f64(h01,0);
  if(t0<0){double s=-h0/t0; if(s<MinPos){MinCos=t0;MinPos=s;MinPos_Idx=0;}}
  double t1=vgetq_lane_f64(t01,1),h1=vgetq_lane_f64(h01,1);
  if(t1<0){double s=-h1/t1; if(s<MinPos){MinCos=t1;MinPos=s;MinPos_Idx=1;}}
  double t2=vgetq_lane_f64(t23,0),h2=vgetq_lane_f64(h23,0);
  if(t2<0){double s=-h2/t2; if(s<MinPos){MinCos=t2;MinPos=s;MinPos_Idx=2;}}
  double t3=vgetq_lane_f64(t23,1),h3=vgetq_lane_f64(h23,1);
  if(t3<0){double s=-h3/t3; if(s<MinPos){MinCos=t3;MinPos=s;MinPos_Idx=3;}}
#else
  double TT[4], height;
  TT[0]=Norm[0]*Photon.UX+Norm[1]*Photon.UY+Norm[2]*Photon.UZ;
  TT[1]=Norm[4]*Photon.UX+Norm[5]*Photon.UY+Norm[6]*Photon.UZ;
  TT[2]=Norm[8]*Photon.UX+Norm[9]*Photon.UY+Norm[10]*Photon.UZ;
  TT[3]=Norm[12]*Photon.UX+Norm[13]*Photon.UY+Norm[14]*Photon.UZ;
  for(int i=0;i<=3;i++) if(TT[i]<0){
    height=Norm[i*4]*Photon.X+Norm[i*4+1]*Photon.Y+Norm[i*4+2]*Photon.Z+Norm[i*4+3];
    double s=-height/TT[i];
    if(s<MinPos){MinCos=TT[i];MinPos=s;MinPos_Idx=i;}
  }
#endif
}

// +++ Henyey-Greenstein scattering +++

static inline void ScatterPhoton(TPhoton& Photon, double g, timos::RngPool& rng){
  if(g >= 1.0) return;
  double u = rng.get_uniform();
  double cost;
  if(g == 0.0){
    cost = u*2.0 - 1.0;
  }else{
    double tmp = (1.0 - g*g) / (1.0 - g + 2.0*g*u);
    cost = (1.0 + g*g - tmp*tmp) / (2.0*g);
  }
  double sinp, cosp;
  rng.get_sincos(sinp, cosp);
  double ux=Photon.UX, uy=Photon.UY, uz=Photon.UZ;
  double sint = sqrt(1.0 - cost*cost);
  if(fabs(uz) <= G_COS_0_D){
    double temp1 = sqrt(1.0 - uz*uz);
    double temp  = sint / temp1;
    double temp2 = uz * cosp;
    Photon.UX = (ux*temp2 - uy*sinp)*temp + ux*cost;
    Photon.UY = (uy*temp2 + ux*sinp)*temp + uy*cost;
    Photon.UZ = -sint*cosp*temp1 + uz*cost;
  }else{
    Photon.UX = sint * cosp;
    Photon.UY = sint * sinp;
    Photon.UZ = (uz > 0) ? cost : -cost;
  }
}

// +++ Main photon propagation thread function +++

void* ThreadPhotonPropagation(void* threadid){
  unsigned int tid = (unsigned int)(uintptr_t)threadid;

  double* LocalAbsorption = new double[g_NumElem + 1]();
  double* LocalSurfMeas   = new double[g_NumBoundaryTrig + 1]();

  TPhotonInfo* PhotonInfo   = g_TimeDomain ? new TPhotonInfo[TD_FLUSH_SIZE+4] : nullptr;
  int          PhotonInfo_Idx = 0;
  if(g_TimeDomain) PhotonInfo[0].Idx = -1;

  double*  randnums_xyz = new double[RNG_BUF_XYZ3+6];
  int      Idx_U = 0;

  uint64_t seed_state = ((uint64_t)(tid + g_StartRandIdx) * 0x9E3779B97F4A7C15ULL)
                        ^ (uint64_t)time(nullptr);
  uint64_t seed_a = timos::splitmix64(seed_state);
  uint64_t seed_b = timos::splitmix64(seed_state);
  timos::RngPool      rng(seed_a);
  timos::Xoshiro256ss rng_init(seed_b);
  rng_fill_uniform(rng_init, randnums_xyz, RNG_BUF_XYZ3+6, -1.0, 1.0);

  TSource   ThreadSources; ThreadSources.NumPhoton = 0;
  TPhoton   Photon;
  double    NX, NY, NZ, Rest_Step, temp, temp1;
  long long int NumSteps = 0, NumIntersections = 0;

  TMedOptic* ThreadMedOptic = new TMedOptic[g_NumMed+1];
  for(int i = 0; i <= g_NumMed; i++) ThreadMedOptic[i] = g_MedOptic[i];

  // ------------------------------------------------------------------ main loop
  do{
  Label_NewPhoton:
    if(ThreadSources.NumPhoton <= 0){
      int n = Fetch_Source(ThreadSources, g_SourceIdx, tid);
      if(n == 0){ if(tid==0) cerr << "\r100%\n"; break; }
      g_SimedPhoton += n;
      if(tid==0) cerr << "\r" << (int)(g_SimedPhoton*100/g_TotalPhoton) << "%";
    }
    ThreadSources.NumPhoton--;

    switch(ThreadSources.SourceType){
    case 11: InitialPhoton_PencilSource(Photon, ThreadSources); break;
    case 12: InitialPhoton_TriangleSource(Photon, ThreadSources, randnums_xyz, Idx_U, rng_init); break;
    case  1: InitialPhoton_PointSource  (Photon, ThreadSources, randnums_xyz, Idx_U, rng_init); break;
    case  2: InitialPhoton_RegionSource (Photon, ThreadSources, randnums_xyz, Idx_U, rng_init); break;
    default: continue;
    }

    // ---------------------------------------------------------- photon loop
    do{
      NumSteps++;
      Photon.Step = rng.get_neg_log() * ThreadMedOptic[Photon.Cur_Med].IMUAMUS;

      double MinCos, MinPos;
      int    MinPos_Idx;

    Label_100:
      NumIntersections++;
      PhotonTetrahedronIntersection(MinPos, MinPos_Idx, MinCos,
                                    g_Elems[Photon.Cur_Elem].TriNorm, Photon);

      if(Photon.Step > MinPos){
        bool   IsBTriangle  = false;
        int    OtherMed     = 0, OtherElem = 0;
        double OtherRefIdx  = 0.0, RefPercentage = 0.0;

        // Move to hit triangle
        Photon.X    += MinPos * Photon.UX;
        Photon.Y    += MinPos * Photon.UY;
        Photon.Z    += MinPos * Photon.UZ;
        Photon.Path += MinPos * ThreadMedOptic[Photon.Cur_Med].RefIdx;
        Rest_Step    = (Photon.Step - MinPos) * ThreadMedOptic[Photon.Cur_Med].MUAMUS;

        if(g_Elems[Photon.Cur_Elem].AdjElemIdx[MinPos_Idx] > 0){
          OtherElem  = g_Elems[Photon.Cur_Elem].AdjElemIdx[MinPos_Idx];
          OtherMed   = g_Elems[OtherElem].MedIdx;
          OtherRefIdx = ThreadMedOptic[OtherMed].RefIdx;
        }else{
          IsBTriangle = true;
          if(g_UniformBoundary == 0){ OtherRefIdx = g_EnvRefIdx; OtherMed = 0; }
          else                      { OtherMed = Photon.Cur_Med; OtherRefIdx = ThreadMedOptic[Photon.Cur_Med].RefIdx; }
        }

        if(OtherMed == Photon.Cur_Med || OtherRefIdx == ThreadMedOptic[Photon.Cur_Med].RefIdx){
          // No reflection/refraction
          if(IsBTriangle){
            if(g_TimeDomain){
              int TT = (unsigned int)floor(Photon.Path * g_InvLightSpeedMutTimeStep);
              if(TT >= g_NumTimeStep) goto Label_NewPhoton;
              TT = (TT<<1)+1;
              int Ci = -g_Elems[Photon.Cur_Elem].AdjElemIdx[MinPos_Idx];
              if(Ci == PhotonInfo[PhotonInfo_Idx].Idx && TT == PhotonInfo[PhotonInfo_Idx].Time_Type)
                PhotonInfo[PhotonInfo_Idx].LostWeight += Photon.Weight;
              else{ PhotonInfo_Idx++; PhotonInfo[PhotonInfo_Idx]={TT,Ci,Photon.Weight}; }
              if(PhotonInfo_Idx >= TD_FLUSH_SIZE) SaveLocal2Global(PhotonInfo, PhotonInfo_Idx);
            }else{
              LocalSurfMeas[-g_Elems[Photon.Cur_Elem].AdjElemIdx[MinPos_Idx]] += Photon.Weight;
            }
            break;
          }else{
            Photon.Cur_Elem = OtherElem;
            Photon.Cur_Med  = OtherMed;
            Photon.Step     = Rest_Step * ThreadMedOptic[Photon.Cur_Med].IMUAMUS;
            goto Label_100;
          }
        }else{
          // Possible reflection / refraction (Fresnel)
          double cosa, sina, cosb = 1.0, sinb;
          NX = g_Elems[Photon.Cur_Elem].TriNorm[MinPos_Idx   ];
          NY = g_Elems[Photon.Cur_Elem].TriNorm[MinPos_Idx+ 4];
          NZ = g_Elems[Photon.Cur_Elem].TriNorm[MinPos_Idx+ 8];
          cosa = -MinCos;
          if(cosa > G_COS_0_D){
            double r = (ThreadMedOptic[Photon.Cur_Med].RefIdx - OtherRefIdx)
                     / (ThreadMedOptic[Photon.Cur_Med].RefIdx + OtherRefIdx);
            RefPercentage = r*r;
          }else if(cosa < G_COS_90_D){
            RefPercentage = 1.0;
          }else{
            sina = sqrt(1.0 - cosa*cosa);
            sinb = ThreadMedOptic[Photon.Cur_Med].RefIdx * sina / OtherRefIdx;
            if(sinb >= 1.0){
              RefPercentage = 1.0;
            }else{
              cosb = sqrt(1.0 - sinb*sinb);
              double cab=cosa*cosb, sab=sina*sinb, scb=sina*cosb, csb=cosa*sinb;
              double cap=cab-sab, cam=cab+sab, sap=scb+csb, sam=scb-csb;
              RefPercentage = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam);
            }
          }

          if(rng.get_uniform() <= RefPercentage){
            // Reflect
            Photon.Step = Photon.Step - MinPos;
            Photon.UX = 2*cosa*NX + Photon.UX;
            Photon.UY = 2*cosa*NY + Photon.UY;
            Photon.UZ = 2*cosa*NZ + Photon.UZ;
            goto Label_100;
          }else{
            // Transmit
            if(IsBTriangle){
              if(g_TimeDomain){
                int TT = (unsigned int)floor(Photon.Path * g_InvLightSpeedMutTimeStep);
                if(TT >= g_NumTimeStep) goto Label_NewPhoton;
                TT = (TT<<1)+1;
                int Ci = -g_Elems[Photon.Cur_Elem].AdjElemIdx[MinPos_Idx];
                if(Ci == PhotonInfo[PhotonInfo_Idx].Idx && TT == PhotonInfo[PhotonInfo_Idx].Time_Type)
                  PhotonInfo[PhotonInfo_Idx].LostWeight += Photon.Weight;
                else{ PhotonInfo_Idx++; PhotonInfo[PhotonInfo_Idx]={TT,Ci,Photon.Weight}; }
                if(PhotonInfo_Idx >= TD_FLUSH_SIZE) SaveLocal2Global(PhotonInfo, PhotonInfo_Idx);
              }else{
                LocalSurfMeas[-g_Elems[Photon.Cur_Elem].AdjElemIdx[MinPos_Idx]] += Photon.Weight;
              }
              break;
            }else{
              double sbdivsa = ThreadMedOptic[Photon.Cur_Med].RefIdx / OtherRefIdx;
              Photon.Step = Rest_Step * ThreadMedOptic[OtherMed].IMUAMUS;
              Photon.UX = -cosb*NX + sbdivsa*(cosa*NX + Photon.UX);
              Photon.UY = -cosb*NY + sbdivsa*(cosa*NY + Photon.UY);
              Photon.UZ = -cosb*NZ + sbdivsa*(cosa*NZ + Photon.UZ);
              Photon.Cur_Elem = OtherElem;
              Photon.Cur_Med  = OtherMed;
              goto Label_100;
            }
          }
        }
      } // end if(Photon.Step > MinPos)

      // Move photon (no boundary crossing this step)
      Photon.X    += Photon.UX * Photon.Step;
      Photon.Y    += Photon.UY * Photon.Step;
      Photon.Z    += Photon.UZ * Photon.Step;
      Photon.Path += Photon.Step * ThreadMedOptic[Photon.Cur_Med].RefIdx;

      // Absorption
      temp = ThreadMedOptic[Photon.Cur_Med].pdwa * Photon.Weight;
      if(g_TimeDomain){
        int TT = (unsigned int)floor(Photon.Path * g_InvLightSpeedMutTimeStep);
        if(TT >= g_NumTimeStep) goto Label_NewPhoton;
        TT = (TT<<1);
        int Ci = Photon.Cur_Elem;
        if(Ci == PhotonInfo[PhotonInfo_Idx].Idx && TT == PhotonInfo[PhotonInfo_Idx].Time_Type)
          PhotonInfo[PhotonInfo_Idx].LostWeight += temp;
        else{ PhotonInfo_Idx++; PhotonInfo[PhotonInfo_Idx]={TT,Ci,temp}; }
        if(PhotonInfo_Idx >= TD_FLUSH_SIZE) SaveLocal2Global(PhotonInfo, PhotonInfo_Idx);
      }else{
        LocalAbsorption[Photon.Cur_Elem] += temp;
      }

#ifdef NOPUREGLASS
      Photon.Weight -= temp;
#else
      if(ThreadMedOptic[Photon.Cur_Med].g < 1.0) Photon.Weight -= temp;
#endif

      if(Photon.Weight <= WEIGHT_THRESHOLD){
        if(rng.get_uniform() < ROULETTE_SURVIVE) Photon.Weight *= ROULETTE_BOOST;
        else break;
      }

      ScatterPhoton(Photon, ThreadMedOptic[Photon.Cur_Med].g, rng);

    }while(true);
  }while(true);

  // Flush and merge
  if(g_TimeDomain){
    if(PhotonInfo_Idx > 0) SaveLocal2Global(PhotonInfo, PhotonInfo_Idx);
  }else{
    os_unfair_lock_lock(&Result_Lock);
    for(int i = 1; i <= g_NumElem;         i++) g_Absorption[i] += LocalAbsorption[i];
    for(int i = 1; i <= g_NumBoundaryTrig; i++) g_SurfMeas[i]   += LocalSurfMeas[i];
    os_unfair_lock_unlock(&Result_Lock);
  }
  g_NumIntersections += NumIntersections;
  g_NumSteps         += NumSteps;

  delete[] randnums_xyz;
  delete[] ThreadMedOptic;
  delete[] LocalAbsorption;
  delete[] LocalSurfMeas;
  if(PhotonInfo) delete[] PhotonInfo;
  return nullptr;
}
