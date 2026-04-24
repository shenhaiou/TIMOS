#pragma once
#include "timos.h"
#include <string>

int ReadOpticalParameter(const std::string& filename,
                         bool&       simpleOptic,
                         int&        numMed,
                         int&        uniformBoundary,
                         double&     envRefIdx,
                         TMedOptic*& medOptic);
