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

#include "Enc_MTotalizer.h"

using namespace openwbo;

/*_________________________________________________________________________________________________
  |
  |  encode : (S : Solver *) (lits : vec<Lit>&) (rhs : int64_t) ->  [void]
  |
  |  Description:
  |
  |     Encodes that at most 'rhs' literals can be assigned value true.
  |     Uses the Modulo Totalizer encoding for translating the cardinality
  |     constraint into CNF.
  |
  |  For further details see:
  |    *  Toru Ogawa, Yangyang Liu, Ryuzo Hasegawa, Miyuki Koshimura,
  |       Hiroshi Fujita:
  |       Modulo Based CNF Encoding of Cardinality Constraints and
  |       Its Application to MaxSAT Solvers. ICTAI 2013: 9-17
  |
  |  Pre-conditions:
  |    * Assumes that 'lits' is not empty.
  |    * Assumes that 'rhs' is larger or equal to 0.
  |
  |  Post-conditions:
  |    * 'S' is updated with the clauses that encode the cardinality constraint.
  |
  |________________________________________________________________________________________________@*/
void MTotalizer::encode(Solver *S, vec<Lit> &lits, int64_t rhs) {
  assert(lits.size() > 0);
  hasEncoding = false;

  cardinality_upoutlits.clear();
  cardinality_lwoutlits.clear();

  if (rhs == 0) {
    for (int i = 0; i < lits.size(); i++)
      addUnitClause(S, ~lits[i]);
    return;
  }

  assert(rhs >= 1 && rhs <= lits.size());

  if (rhs == lits.size()) {
    return;
  }

  hasEncoding = true;

  // Modulo that was used in the original paper.
  // int mod = ceil(sqrt(lits.size()));
  // Slightly better results than using the above modulo.
  int mod = ceil(sqrt(rhs + 1));
  if (modulo == -1)
    modulo = mod;
  else
    mod = modulo;

  for (int i = 0; i < floor(lits.size() / mod); i++) {
    Lit p = mkLit(S->nVars(), false);
    newSATVariable(S);
    cardinality_upoutlits.push(p);
  }

  for (int i = 0; i < mod - 1; i++) {
    Lit p = mkLit(S->nVars(), false);
    newSATVariable(S);
    cardinality_lwoutlits.push(p);
  }

  lits.copyTo(cardinality_inlits);
  current_cardinality_rhs = rhs + 1;

  if (cardinality_upoutlits.size() == 0)
    cardinality_upoutlits.push(h0);

  toCNF(S, mod, cardinality_upoutlits, cardinality_lwoutlits, lits.size());
  assert(cardinality_inlits.size() == 0);

  // Limit the rhs.
  update(S, rhs);
}

/*_________________________________________________________________________________________________
  |
  |  update : (S : Solver *) (rhs : int64_t) ->  [void]
  |
  |  Description:
  |
  |     Updates the 'rhs' of an already existent cardinality encoding.
  |     This method allows for all learned clauses from previous iterations to
  |     be kept in the next iteration.
  |     This procedure is denoted by incremental strengthening.
  |
  |  For further details see:
  |    *  Toru Ogawa, Yangyang Liu, Ryuzo Hasegawa, Miyuki Koshimura,
  |       Hiroshi Fujita:
  |       Modulo Based CNF Encoding of Cardinality Constraints and
  |       Its Application to MaxSAT Solvers. ICTAI 2013: 9-17
  |
  |  Pre-conditions:
  |    * Only one cardinality constraint is being used.
  |    * Assumes that 'current_cardinality_rhs' is larger than -1, i.e. assumes
  |      a cardinality constraint has already been encoded.
  |
  |  Post-conditions:
  |    * 'S' is updated with unit clauses that update the 'rhs' of the
  |      cardinality encoding.
  |    * 'current_cardinality_rhs' is updated.
  |
  |________________________________________________________________________________________________@*/
void MTotalizer::update(Solver *S, int64_t rhs) {
  assert(current_cardinality_rhs != -1);
  assert(hasEncoding);
  encode_output(S, rhs);
  current_cardinality_rhs = rhs + 1;
}

/************************************************************************************************
//
// Auxiliary methods for the cardinality encoding
//
************************************************************************************************/

void MTotalizer::encode_output(Solver *S, int64_t rhs) {

  assert(hasEncoding);
  assert(cardinality_upoutlits.size() != 0 ||
         cardinality_lwoutlits.size() != 0);

  int mod = modulo;

  int ulimit = floor((rhs + 1) / mod);
  int llimit = (rhs + 1) - ulimit * mod;

  assert(ulimit <= cardinality_upoutlits.size());
  assert(llimit <= cardinality_lwoutlits.size());

  for (int i = ulimit; i < cardinality_upoutlits.size(); i++)
    addUnitClause(S, ~cardinality_upoutlits[i]);

  if (ulimit != 0 && llimit != 0) {
    for (int i = llimit - 1; i < cardinality_lwoutlits.size(); i++)
      addBinaryClause(S, ~cardinality_upoutlits[ulimit - 1],
                      ~cardinality_lwoutlits[i]);
  } else {
    if (ulimit == 0) {
      assert(llimit != 0);
      for (int i = llimit - 1; i < cardinality_lwoutlits.size(); i++) {
        addUnitClause(S, ~cardinality_lwoutlits[i]);
      }
    } else
      addUnitClause(S, ~cardinality_upoutlits[ulimit - 1]);
  }
}

void MTotalizer::toCNF(Solver *S, int mod, vec<Lit> &ublits, vec<Lit> &lwlits,
                       int64_t rhs) {

  vec<Lit> lupper;
  vec<Lit> llower;
  vec<Lit> rupper;
  vec<Lit> rlower;

  assert(rhs > 1);
  int split = floor(rhs / 2);
  int left = 1;
  int right = 1;

  if (split == 1) {
    assert(cardinality_inlits.size() > 0);
    lupper.push(h0);
    llower.push(cardinality_inlits.last());
    cardinality_inlits.pop();
  } else {
    left = floor(split / mod);
    for (int i = 0; i < left; i++) {
      Lit p = mkLit(S->nVars(), false);
      newSATVariable(S);
      lupper.push(p);
    }
    int limit = mod - 1;
    if (left % mod == 0 && split < mod - 1) {
      limit = split;
    }

    for (int i = 0; i < limit; i++) {
      Lit p = mkLit(S->nVars(), false);
      newSATVariable(S);
      llower.push(p);
    }
  }

  if (rhs - split == 1) {
    assert(cardinality_inlits.size() > 0);
    rupper.push(h0);
    rlower.push(cardinality_inlits.last());
    cardinality_inlits.pop();
  } else {
    right = floor((rhs - split) / mod);
    for (int i = 0; i < right; i++) {
      Lit p = mkLit(S->nVars(), false);
      newSATVariable(S);
      rupper.push(p);
    }
    int limit = mod - 1;
    if (right % mod == 0 && rhs - split < mod - 1) {
      limit = rhs - split;
    }
    for (int i = 0; i < limit; i++) {
      Lit p = mkLit(S->nVars(), false);
      newSATVariable(S);
      rlower.push(p);
    }
  }

  if (lupper.size() == 0) {
    lupper.push(h0);
  }

  if (rupper.size() == 0) {
    rupper.push(h0);
  }

  adder(S, mod, ublits, lwlits, rupper, rlower, lupper, llower);
  if (left * mod + split - left * mod > 1)
    toCNF(S, mod, lupper, llower, left * mod + split - left * mod);
  if (right * mod + (rhs - split) - right * mod > 1)
    toCNF(S, mod, rupper, rlower, right * mod + (rhs - split) - right * mod);
}

void MTotalizer::adder(Solver *S, int mod, vec<Lit> &upper, vec<Lit> &lower,
                       vec<Lit> &lupper, vec<Lit> &llower, vec<Lit> &rupper,
                       vec<Lit> &rlower) {

  assert(upper.size() != 0);
  assert(lower.size() >= llower.size() && lower.size() >= rlower.size());

  Lit carry = lit_Undef;

  if (upper[0] != h0) {
    carry = mkLit(S->nVars(), false);
    newSATVariable(S);
  }

  for (int i = 0; i <= llower.size(); i++) {
    for (int j = 0; j <= rlower.size(); j++) {

      if (i + j > current_cardinality_rhs + 1 &&
          current_cardinality_rhs + 1 < modulo) {
        continue;
      }

      if (i + j < mod) {
        if (i == 0 && j != 0) {
          if (upper[0] != h0) // The sum will always be smaller than 'mod'.
            addTernaryClause(S, ~rlower[j - 1], lower[i + j - 1], carry);
          else
            addBinaryClause(S, ~rlower[j - 1], lower[i + j - 1]);
        } else if (j == 0 && i != 0) {
          if (upper[0] != h0) // The sum will always be smaller than 'mod'.
            addTernaryClause(S, ~llower[i - 1], lower[i + j - 1], carry);
          else
            addBinaryClause(S, ~llower[i - 1], lower[i + j - 1]);
        } else if (i == 0 && j == 0)
          continue;
        else {
          if (upper[0] != h0) // The sum will always be smaller than 'mod'.
            addQuaternaryClause(S, ~llower[i - 1], ~rlower[j - 1],
                                lower[i + j - 1], carry);
          else {
            assert(i + j - 1 < lower.size());
            addTernaryClause(S, ~llower[i - 1], ~rlower[j - 1],
                             lower[i + j - 1]);
          }
        }
      } else if (i + j > mod) {
        assert(i + j > 0);
        addTernaryClause(S, ~llower[i - 1], ~rlower[j - 1],
                         lower[(i + j) % mod - 1]);
      } else {
        // i + j == mod
        assert(i + j == mod);
        assert(carry != lit_Undef);
        addTernaryClause(S, ~llower[i - 1], ~rlower[j - 1], carry);
      }
    }
  }

  if (upper[0] != h0) {

    for (int i = 0; i <= lupper.size(); i++) {
      for (int j = 0; j <= rupper.size(); j++) {

        Lit a = lit_Error; // lupper
        Lit b = lit_Error; // rupper
        Lit c = lit_Error; // upper(i+j)
        Lit d = lit_Error; // upper(i+j+1)

        int close_mod = floor(current_cardinality_rhs / mod);
        if (current_cardinality_rhs % mod != 0)
          close_mod++;
        if (mod * (i + j) > close_mod * mod)
          continue;

        if (i != 0)
          a = lupper[i - 1];

        if (j != 0)
          b = rupper[j - 1];

        if (i + j != 0 && i + j - 1 < upper.size())
          c = upper[i + j - 1];

        if (i + j < upper.size()) {
          d = upper[i + j];
        }

        if (c != lit_Undef && c != lit_Error) {
          vec<Lit> clause;
          if (a != lit_Undef && a != lit_Error)
            clause.push(~a);
          if (b != lit_Undef && b != lit_Error)
            clause.push(~b);

          clause.push(c);
          if (clause.size() > 1) {
            S->addClause(clause);
          }
        }

        vec<Lit> clause;
        clause.push(~carry);
        if (a != lit_Undef && a != lit_Error)
          clause.push(~a);
        if (b != lit_Undef && b != lit_Error)
          clause.push(~b);
        if (d != lit_Error && d != lit_Undef)
          clause.push(d);

        if (clause.size() > 1) {
          S->addClause(clause);
        }
      }
    }
  }
}
