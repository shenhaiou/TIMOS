#pragma once
#include "simcontext.h"

// Thread payload: tid is the thread index; ctx is the shared simulation state.
struct ThreadArg { int tid; SimContext* ctx; };

void* ThreadPhotonPropagation(void* arg);
