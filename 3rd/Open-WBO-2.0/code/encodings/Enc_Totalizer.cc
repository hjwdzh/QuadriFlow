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

#include "Enc_Totalizer.h"

using namespace openwbo;

void Totalizer::incremental(Solver *S, int64_t rhs) {

  for (int z = 0; z < totalizerIterative_rhs.size(); z++) {

    // We only need to count the sums up to k.
    for (int i = 0; i <= totalizerIterative_left[z].size(); i++) {
      for (int j = 0; j <= totalizerIterative_right[z].size(); j++) {

        if (i == 0 && j == 0) {
          continue;
        }

        if (i + j > rhs + 1 || i + j <= totalizerIterative_rhs[z] + 1) {
          continue;
        }

        if (i == 0) {
          addBinaryClause(S, ~(totalizerIterative_right[z])[j - 1],
                          (totalizerIterative_output[z])[j - 1], blocking);
          n_clauses++;
        } else if (j == 0) {
          addBinaryClause(S, ~(totalizerIterative_left[z])[i - 1],
                          (totalizerIterative_output[z])[i - 1], blocking);
          n_clauses++;
        } else {
          addTernaryClause(S, ~(totalizerIterative_left[z])[i - 1],
                           ~(totalizerIterative_right[z])[j - 1],
                           (totalizerIterative_output[z])[i + j - 1], blocking);
          n_clauses++;
        }
      }
    }
    totalizerIterative_rhs[z] = rhs;
  }
}

void Totalizer::join(Solver *S, vec<Lit> &lits, int64_t rhs) {

  assert(incremental_strategy == _INCREMENTAL_ITERATIVE_);

  vec<Lit> left_cardinality_outlits;
  cardinality_outlits.copyTo(left_cardinality_outlits);
  int old_cardinality = current_cardinality_rhs;

  if (lits.size() > 1) {
    build(S, lits, rhs < lits.size() ? rhs : lits.size());
  } else {
    assert(lits.size() == 1);
    cardinality_outlits.clear();
    cardinality_outlits.push(lits[0]);
  }

  vec<Lit> right_cardinality_outlits;
  cardinality_outlits.copyTo(right_cardinality_outlits);
  cardinality_outlits.clear();

  for (int i = 0;
       i < left_cardinality_outlits.size() + right_cardinality_outlits.size();
       i++) {
    Lit p = mkLit(S->nVars(), false);
    newSATVariable(S);
    n_variables++;
    cardinality_outlits.push(p);
  }

  current_cardinality_rhs = rhs;
  // TO_adder is using the 'current_cardinality_rhs' value
  adder(S, left_cardinality_outlits, right_cardinality_outlits,
        cardinality_outlits);
  current_cardinality_rhs = old_cardinality;

  for (int i = 0; i < lits.size(); i++)
    ilits.push(lits[i]);
}

void Totalizer::adder(Solver *S, vec<Lit> &left, vec<Lit> &right,
                      vec<Lit> &output) {
  assert(output.size() == left.size() + right.size());
  if (incremental_strategy == _INCREMENTAL_ITERATIVE_) {
    totalizerIterative_left.push();
    new (&totalizerIterative_left[totalizerIterative_left.size() - 1])
        vec<Lit>();
    left.copyTo(totalizerIterative_left.last());
    totalizerIterative_right.push();
    new (&totalizerIterative_right[totalizerIterative_right.size() - 1])
        vec<Lit>();
    right.copyTo(totalizerIterative_right.last());
    totalizerIterative_output.push();
    new (&totalizerIterative_output[totalizerIterative_output.size() - 1])
        vec<Lit>();
    output.copyTo(totalizerIterative_output.last());
    totalizerIterative_rhs.push(current_cardinality_rhs);
  }

  // We only need to count the sums up to k.
  for (int i = 0; i <= left.size(); i++) {
    for (int j = 0; j <= right.size(); j++) {
      if (i == 0 && j == 0)
        continue;

      if (i + j > current_cardinality_rhs + 1)
        continue;

      if (i == 0) {
        addBinaryClause(S, ~right[j - 1], output[j - 1], blocking);
        n_clauses++;
      } else if (j == 0) {
        addBinaryClause(S, ~left[i - 1], output[i - 1], blocking);
        n_clauses++;
      } else {
        addTernaryClause(S, ~left[i - 1], ~right[j - 1], output[i + j - 1],
                         blocking);
        n_clauses++;
      }
    }
  }
}

void Totalizer::toCNF(Solver *S, vec<Lit> &lits) {

  vec<Lit> left;
  vec<Lit> right;

  assert(lits.size() > 1);
  int split = floor(lits.size() / 2);

  for (int i = 0; i < lits.size(); i++) {

    if (i < split) {
      // left branch
      if (split == 1) {
        assert(cardinality_inlits.size() > 0);
        left.push(cardinality_inlits.last());
        cardinality_inlits.pop();
      } else {
        Lit p = mkLit(S->nVars(), false);
        newSATVariable(S);
        left.push(p);
      }
    } else {

      // right branch
      if (lits.size() - split == 1) {
        assert(cardinality_inlits.size() > 0);
        right.push(cardinality_inlits.last());
        cardinality_inlits.pop();
      } else {
        Lit p = mkLit(S->nVars(), false);
        newSATVariable(S);
        right.push(p);
      }
    }
  }

  if (left.size() > 1)
    toCNF(S, left);
  if (right.size() > 1)
    toCNF(S, right);
  adder(S, left, right, lits);
}

void Totalizer::update(Solver *S, int64_t rhs, vec<Lit> &lits,
                       vec<Lit> &assumptions) {

  assert(hasEncoding);

  switch (incremental_strategy) {
  case _INCREMENTAL_NONE_:
    for (int i = rhs; i < cardinality_outlits.size(); i++)
      addUnitClause(S, ~cardinality_outlits[i]);
    break;

  case _INCREMENTAL_BLOCKING_:
    assumptions.clear();

    // Disable previous iterations
    for (int i = 0; i < disable_lits.size(); i++)
      addUnitClause(S, disable_lits[i]);

    build(S, lits, rhs);
    if (blocking != lit_Undef)
      assumptions.push(~blocking);

    for (int i = rhs; i < cardinality_outlits.size(); i++)
      addUnitClause(S, ~cardinality_outlits[i]);
    break;

  case _INCREMENTAL_WEAKENING_:
    // Change the set of assumption variables
    assumptions.clear();
    for (int i = rhs; i < cardinality_outlits.size(); i++)
      assumptions.push(~cardinality_outlits[i]);
    break;

  case _INCREMENTAL_ITERATIVE_:
    incremental(S, rhs);
    assumptions.clear();
    for (int i = rhs; i < cardinality_outlits.size(); i++)
      assumptions.push(~cardinality_outlits[i]);
    break;

  default:
    printf("Error: No incremental strategy.\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }
}

void Totalizer::add(Solver *S, Totalizer &tot, int64_t rhs) {
  assert(incremental_strategy == _INCREMENTAL_ITERATIVE_ &&
         tot.incremental_strategy == _INCREMENTAL_ITERATIVE_);
  int left_idx = totalizerIterative_rhs.size() - 1;
  for (int i = 0; i < tot.totalizerIterative_rhs.size(); ++i) {
    totalizerIterative_left.push();
    new (&totalizerIterative_left[totalizerIterative_left.size() - 1])
        vec<Lit>();
    tot.totalizerIterative_left[i].copyTo(totalizerIterative_left.last());
    totalizerIterative_right.push();
    new (&totalizerIterative_right[totalizerIterative_right.size() - 1])
        vec<Lit>();
    tot.totalizerIterative_right[i].copyTo(totalizerIterative_right.last());
    totalizerIterative_output.push();
    new (&totalizerIterative_output[totalizerIterative_output.size() - 1])
        vec<Lit>();
    tot.totalizerIterative_output[i].copyTo(totalizerIterative_output.last());
    totalizerIterative_rhs.push(tot.totalizerIterative_rhs[i]);
  }
  int right_idx = totalizerIterative_rhs.size() - 1;

  vec<Lit> left, right;
  totalizerIterative_output[left_idx].copyTo(left);
  totalizerIterative_output[right_idx].copyTo(right);
  cardinality_outlits.clear();
  for (int i = 0; i < left.size() + right.size(); ++i) {
    Lit p = mkLit(S->nVars(), false);
    newSATVariable(S);
    cardinality_outlits.push(p);
  }
  current_cardinality_rhs = rhs;
  adder(S, left, right, cardinality_outlits);
}

/*_________________________________________________________________________________________________
  |
  |  build : (S : Solver *) (lits : vec<Lit>&) (int64_t : rhs)  ->  [void]
  |
  |  Description:
  |
  |    Builds a cardinality constraint of the kind x_1 + ... x_n <= k.
  |    Uses the totalizer encoding for translation the cardinality constraint
  into CNF.
  |    Does not impose any constraints on the value of 'k'.
  |    NOTE: Use method 'update' to impose a restriction on the value of 'k'.
  |
  |  For further details see:
  |    * Olivier Bailleux, Yacine Boufkhad: Efficient CNF Encoding of Boolean
  Cardinality Constraints. CP 2003: 108-122
  |
  |  Pre-conditions:
  |    * Assumes that 'lits' is not empty and 'rhs' is larger or equal to 0.
  |
  |  Post-conditions:
  |    * 'S' is updated with the clauses that encode the cardinality constraint.
  |    * hasEncoding is set to 'true'.
  |
  |________________________________________________________________________________________________@*/
void Totalizer::build(Solver *S, vec<Lit> &lits, int64_t rhs) {

  cardinality_outlits.clear();
  hasEncoding = false;

  if (rhs == 0) {
    for (int i = 0; i < lits.size(); i++)
      addUnitClause(S, ~lits[i]);
    return;
  }

  assert(rhs >= 1 && rhs <= lits.size());

  if (incremental_strategy == _INCREMENTAL_NONE_ && rhs == lits.size()) {
    return;
  }

  if (rhs == lits.size() && !joinMode)
    return;

  for (int i = 0; i < lits.size(); i++) {
    Lit p = mkLit(S->nVars(), false);
    newSATVariable(S);
    cardinality_outlits.push(p);
  }

  lits.copyTo(cardinality_inlits);
  current_cardinality_rhs = rhs;

  // If incremental blocking is enable then all clauses will contain a blocking
  // literal 'b'. Setting this literal to 'true' is the same as deletting these
  // clauses.
  if (incremental_strategy == _INCREMENTAL_BLOCKING_) {
    Lit p = mkLit(S->nVars(), false);
    newSATVariable(S);
    if (blocking != lit_Undef)
      disable_lits.push(blocking);
    blocking = p;
  }

  toCNF(S, cardinality_outlits);
  assert(cardinality_inlits.size() == 0);

  if (!joinMode)
    joinMode = true;
  hasEncoding = true;

  lits.copyTo(ilits);
}
