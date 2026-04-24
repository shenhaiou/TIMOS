#pragma once
#include "timos.h"
#include <string>

int fem_read(const std::string& filename,
             bool         simpleOptic,
             int&         numNode,
             int&         numElem,
             int&         numMed,
             TNode*&      nodes,
             TElemNode*&  elemNodes,
             TElem*&      elems);

int PreProcessor(int&        numElem,
                 int&        numNode,
                 int&        numMed,
                 int&        numBoundaryTrig,
                 int&        numTrig,
                 TNode*&     nodes,
                 TElemNode*& elemNodes,
                 TElem*&     elems,
                 TTriNode*&  triNodes,
                 TTriangle*& triangles,
                 int*&       boundaryTrigs);
