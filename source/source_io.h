#pragma once
#include "timos.h"
#include <string>

long long int ReadSource(const std::string& filename,
                         int&      numSource,
                         TSource*& sources);

int Prepare_Source(int&        numNode,
                   int&        numElem,
                   int&        numTrig,
                   int&        numSource,
                   TSource*&   sources,
                   TElem*&     elems,
                   TNode*&     nodes,
                   TElemNode*& elemNodes,
                   TTriNode*&  triNodes,
                   TTriangle*& triangles);
