#pragma once
#include "simcontext.h"
#include <string>

long long int ReadSource(const std::string& filename, SimContext& ctx);

int Prepare_Source(SimContext& ctx);
