#pragma once
#include "simcontext.h"
#include <string>

TiResult fem_read(const std::string& filename, SimContext& ctx);

TiResult PreProcessor(SimContext& ctx);
