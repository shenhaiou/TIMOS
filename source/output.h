#pragma once
#include "simcontext.h"
#include <string>

// sufThreshhold / intThreshhold are set by these functions (output params).
void AbsorptionToFluence(SimContext& ctx,
                         double& sufThreshhold, double& intThreshhold);

void TimeAbsorptionToFluence(SimContext& ctx,
                             double& sufThreshhold, double& intThreshhold);

bool WriteResultASCII(const std::string& opt_f, const std::string& fem_f,
                      const std::string& src_f, const std::string& out_f,
                      SimContext& ctx,
                      double sufThreshhold, double intThreshhold,
                      int output_format);

bool TimeWriteResultASCII(const std::string& opt_f, const std::string& fem_f,
                          const std::string& src_f, const std::string& out_f,
                          SimContext& ctx,
                          double sufThreshhold, double intThreshhold,
                          int output_format);
