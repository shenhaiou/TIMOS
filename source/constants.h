#pragma once
// Simulation-wide named constants. All floating-point constants use static const
// (not constexpr) because -ffast-math disables IEEE float constexpr in clang.

#include <cmath>

// Speed of light in vacuum (mm / ns)
static const double LIGHT_SPEED     = 2.99792458e2;
static const double INV_LIGHT_SPEED = 1.0 / LIGHT_SPEED;

// Henyey-Greenstein angle thresholds
static const double G_COS_0_D  = 1.0 - 1.0e-14; // treat as cos(0)
static const double G_COS_90_D = 1.0e-7;          // treat as cos(90°)

// Batch size for the init-direction RNG buffer (xyz triples)
static constexpr int RNG_BUF_XYZ  = 256;
static constexpr int RNG_BUF_XYZ3 = RNG_BUF_XYZ * 3;

// Per-thread time-domain accumulation buffer size (entries before flush)
static constexpr int TD_FLUSH_SIZE = 1024 * 1024;

// Photon batch size fetched from the source queue per lock acquisition
static constexpr int SOURCE_BATCH = 1000;

// Photon weight threshold below which Russian-roulette is applied
static const double WEIGHT_THRESHOLD = 1e-5;
// Survival probability and boost factor for Russian roulette
static const double ROULETTE_SURVIVE = 0.1;
static const double ROULETTE_BOOST   = 1.0 / ROULETTE_SURVIVE;

// Sentinel "no intersection yet" distance (larger than any mesh dimension)
static const double NO_INTERSECTION  = 1e10;
