/*****************************************************************************************[Encodings.cc]
Open-WBO -- Copyright (c) 2013-2015, Ruben Martins, Vasco Manquinho, Ines Lynce

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
************************************************************************************************/

#ifndef Enc_CNetworks_h
#define Enc_CNetworks_h

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "Encodings.h"
#include "core/SolverTypes.h"

namespace openwbo {

class CNetworks : public Encodings {

public:
  CNetworks() {
    current_cardinality_rhs = -1; // -1 corresponds to an unitialized value.
  }
  ~CNetworks() {}

  void encode(Solver *S, vec<Lit> &lits, int64_t rhs);
  void update(Solver *S, int64_t rhs);

  bool hasCreatedEncoding() { return hasEncoding; }

protected:
  // Auxiliary methods for the cardinality network encoding:
  //
  void CN_hmerge(Solver *S, vec<Lit> &a_s, vec<Lit> &b_s, vec<Lit> &c_s);
  void CN_hsort(Solver *S, vec<Lit> &a_s, vec<Lit> &c_s);
  void CN_smerge(Solver *S, vec<Lit> &a_s, vec<Lit> &b_s, vec<Lit> &c_s);
  void CN_encode(Solver *S, vec<Lit> &a_s, vec<Lit> &c_s, int64_t rhs);

  // Stores the current value of the rhs of the cardinality constraint.
  int64_t current_cardinality_rhs;

  // Stores the outputs of the cardinality constraint encoding
  // for incremental solving.
  vec<Lit> cardinality_outlits;
};
} // namespace openwbo

#endif