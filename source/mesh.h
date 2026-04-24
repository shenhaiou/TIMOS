#pragma once
#include "simcontext.h"
#include <string>

int fem_read(const std::string& filename, SimContext& ctx);

int PreProcessor(SimContext& ctx);
