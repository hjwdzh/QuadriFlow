/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *
 * Open-WBO, Copyright (c) 2013-2017, Ruben Martins, Vasco Manquinho, Ines Lynce
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#ifndef Encoder_h
#define Encoder_h

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "MaxTypes.h"
#include "core/SolverTypes.h"

// Encodings
#include "encodings/Enc_CNetworks.h"
#include "encodings/Enc_GTE.h"
#include "encodings/Enc_Ladder.h"
#include "encodings/Enc_MTotalizer.h"
#include "encodings/Enc_SWC.h"
#include "encodings/Enc_Totalizer.h"

using NSPACE::vec;
using NSPACE::Lit;
using NSPACE::Solver;

namespace openwbo {

//=================================================================================================
class Encoder {

public:
  Encoder(int incremental = _INCREMENTAL_NONE_,
          int cardinality = _CARD_TOTALIZER_, int amo = _AMO_LADDER_,
          int pb = _PB_SWC_) {
    pb_encoding = pb;
    amo_encoding = amo;
    incremental_strategy = incremental;
    cardinality_encoding = cardinality;
    totalizer.setIncremental(incremental);
  }

  ~Encoder() {}

  // TEMP
  vec<Lit> &lits();
  vec<Lit> &outputs();

  // At-most-one encodings:
  //
  // Encode exactly-one constraint into CNF.
  void encodeAMO(Solver *S, vec<Lit> &lits);

  // Cardinality encodings:
  //
  // Encode cardinality constraint into CNF.
  void encodeCardinality(Solver *S, vec<Lit> &lits, int64_t rhs);

  // Update the rhs of an already existent cardinality constraint
  void updateCardinality(Solver *S, int64_t rhs);

  // Incremental cardinality encodings:
  //
  // Build a cardinality constraint that can count up to 'rhs'.
  // No restriction is made on the value of 'rhs'.
  // buildCardinality + updateCardinality is equivalent to encodeCardinality.
  // Useful for incremental encodings.
  void buildCardinality(Solver *S, vec<Lit> &lits, int64_t rhs);

  // Incremental update for cardinality constraints;
  void incUpdateCardinality(Solver *S, vec<Lit> &join, vec<Lit> &lits,
                            int64_t rhs, vec<Lit> &assumptions);
  void incUpdateCardinality(Solver *S, vec<Lit> &lits, int64_t rhs,
                            vec<Lit> &assumptions) {

    vec<Lit> empty;
    incUpdateCardinality(S, empty, lits, rhs, assumptions);
  }

  // Add two disjoint cardinality constraints
  void addCardinality(Solver *S, Encoder &enc, int64_t rhs);

  // PB encodings:
  //
  // Encode pseudo-Boolean constraint into CNF.
  void encodePB(Solver *S, vec<Lit> &lits, vec<uint64_t> &coeffs, uint64_t rhs);
  // Update the rhs of an already existent pseudo-Boolean constraint.
  void updatePB(Solver *S, uint64_t rhs);

  // Incremental PB encodings:
  //
  // Incremental PB encoding.
  void incEncodePB(Solver *S, vec<Lit> &lits, vec<uint64_t> &coeffs,
                   int64_t rhs, vec<Lit> &assumptions, int size);

  // Incremental update of PB encodings.
  void incUpdatePB(Solver *S, vec<Lit> &lits, vec<uint64_t> &coeffs,
                   int64_t rhs, vec<Lit> &assumptions);

  // Incremental update of assumptions.
  void incUpdatePBAssumptions(Solver *S, vec<Lit> &assumptions);

  // Incremental construction of the totalizer encoding.
  // Joins a set of new literals, x_1 + ... + x_i, to an existing encoding of
  // the type
  // y_1 + ... + y_j <= k. It also updates 'k' to 'rhs'.
  void joinEncoding(Solver *S, vec<Lit> &lits, int64_t rhs);

  // Other:
  //
  // Returns true if an encoding has been built, false otherwise.
  bool hasCardEncoding();
  bool hasPBEncoding();

  // Controls the type of encoding to be used:
  //
  void setCardEncoding(int enc) { cardinality_encoding = enc; }
  int getCardEncoding() { return cardinality_encoding; }

  void setPBEncoding(int enc) { pb_encoding = enc; }
  int getPBEncoding() { return pb_encoding; }

  void setAMOEncoding(int enc) { amo_encoding = enc; }
  int getAMOEncoding() { return amo_encoding; }

  // Controls the modulo value that is used in the modulo totalizer encoding.
  //
  void setModulo(int m) { mtotalizer.setModulo(m); }
  int getModulo() { return mtotalizer.getModulo(); }

  // Sets the incremental strategy for the totalizer encoding.
  //
  void setIncremental(int incremental) {
    incremental_strategy = incremental;
    totalizer.setIncremental(incremental);
  }

protected:
  int incremental_strategy;
  int cardinality_encoding;
  int pb_encoding;
  int amo_encoding;

  // At-most-one encodings
  Ladder ladder;

  // Cardinality encodings
  CNetworks cnetworks;
  MTotalizer mtotalizer;
  Totalizer totalizer;

  // PB encodings
  SWC swc;
  GTE gte;
};
} // namespace openwbo

#endif
