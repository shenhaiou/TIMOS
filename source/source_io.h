#pragma once
#include "simcontext.h"
#include <string>

// Returns total photon count on success, or an error message.
std::expected<long long int, std::string> ReadSource(const std::string& filename, SimContext& ctx);

TiResult Prepare_Source(SimContext& ctx);
