#pragma once
// Shared I/O helpers used by all three file parsers.
#include <istream>
#include <string>

// Skip whitespace then any number of '#'-to-end-of-line comment lines.
static inline void skip_comments(std::istream& fin) {
  fin >> std::ws;
  while(fin.good() && fin.peek() == '#') {
    std::string line;
    std::getline(fin, line);
    fin >> std::ws;
  }
}
