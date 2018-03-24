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

#include "Enc_SWC.h"

using namespace openwbo;

/*_________________________________________________________________________________________________
  |
  |  encode : (S : Solver *) (lits : vec<Lit>&) (rhs : int64_t) ->  [void]
  |
  |  Description:
  |
  |     Encodes that at most 'rhs' sum of literals (each literal has an
  |     associated weight) can be assigned value true.
  |     Uses the Sequential Weight Counter encoding for translating the
  |     pseudo-Boolean constraint into CNF.
  |
  |  For further details see:
  |    * S. Holldobler, N. Manthey, P. Steinke: A Compact Encoding of
  |      Pseudo-Boolean Constraints into SAT. KI 2012: 107-118
  |
  |  Pre-conditions:
  |    * Assumes that 'rhs' is larger or equal to 0.
  |
  |  Post-conditions:
  |    * 'S' is updated with the clauses that encode the pseudo-Boolean
  |      constraint.
  |
  |________________________________________________________________________________________________@*/
void SWC::encode(Solver *S, vec<Lit> &lits, vec<uint64_t> &coeffs,
                 uint64_t rhs) {
  // FIXME: do not change coeffs in this method. Make coeffs const.

  // If the rhs is larger than INT32_MAX is not feasible to encode this
  // pseudo-Boolean constraint to CNF.
  if (rhs >= INT32_MAX) {
    printf("c Overflow in the Encoding\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }

  hasEncoding = false;
  nb_variables = 0;
  nb_clauses = 0;

  vec<Lit> simp_lits;
  vec<uint64_t> simp_coeffs;
  lits.copyTo(simp_lits);
  coeffs.copyTo(simp_coeffs);

  lits.clear();
  coeffs.clear();

  // Fix literals that have a coeff larger than rhs.
  for (int i = 0; i < simp_lits.size(); i++) {
    if (simp_coeffs[i] == 0)
      continue;

    if (simp_coeffs[i] >= INT32_MAX) {
      printf("c Overflow in the Encoding\n");
      printf("s UNKNOWN\n");
      exit(_ERROR_);
    }

    if (simp_coeffs[i] <= (unsigned)rhs) {
      lits.push(simp_lits[i]);
      coeffs.push(simp_coeffs[i]);
    } else
      addUnitClause(S, ~simp_lits[i]);
  }

  if (lits.size() == 1) {
    // Should this be done?
    // addUnitClause(S, ~lits[0]);
    return;
  }

  if (lits.size() == 0)
    return;

  // Create auxiliary variables.
  int n = lits.size();
  vec<Lit> *seq_auxiliary = new vec<Lit>[n + 1];
  for (int i = 0; i < n + 1; i++)
    seq_auxiliary[i].growTo(rhs + 1);

  for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= (int)rhs; ++j) {
      seq_auxiliary[i][j] = mkLit(S->nVars(), false);
      newSATVariable(S);
      nb_variables++;
    }
  }

  for (int i = 1; i <= (int)rhs; ++i)
    pb_outlits.push(seq_auxiliary[n][i]);

  for (int i = 1; i <= n; i++) {
    // WARNING: wi is used as int for array indexes but as int64_t
    // for the coeffs. Dangerous if the coeffs are larger than INT32_MAX.
    // Same problem occurs with rhs.
    uint64_t wi = coeffs[i - 1];
    // assert(wi <= rhs);

    for (int j = 1; j <= (int)rhs; j++) {
      if (i >= 2 && i <= n && j <= (int)rhs) {
        addBinaryClause(S, ~seq_auxiliary[i - 1][j], seq_auxiliary[i][j]);
        nb_clauses++;
      }
      if (i <= n && j <= (int)wi) {
        addBinaryClause(S, ~lits[i - 1], seq_auxiliary[i][j]);
        nb_clauses++;
      }
      if (i >= 2 && i <= n && j <= (int)(rhs - wi)) {
        addTernaryClause(S, ~seq_auxiliary[i - 1][j], ~lits[i - 1],
                         seq_auxiliary[i][j + (int)wi]);
        nb_clauses++;
      }
    }

    // Encode rhs.
    if (i >= 2) {
      addBinaryClause(S, ~seq_auxiliary[i - 1][(int)rhs + 1 - (int)wi],
                      ~lits[i - 1]);
      nb_clauses++;
    }
  }

  current_pb_rhs = rhs;
  hasEncoding = true;
}

/*_________________________________________________________________________________________________
  |
  |  encode : (S : Solver *) (lits : vec<Lit>&) (rhs : int64_t)
  |           (assumptions: vec<Lit>&) (size: int) ->  [void]
  |
  |  Description:
  |
  |     Incremental construction of the SWC encoding.
  |
  |  For further details see:
  |    * S. Holldobler, N. Manthey, P. Steinke: A Compact Encoding of
  |      Pseudo-Boolean Constraints into SAT. KI 2012: 107-118
  |
  |  Pre-conditions:
  |    * Assumes that 'rhs' is larger or equal to 0.
  |
  |  Post-conditions:
  |    * 'S' is updated with the clauses that encode the pseudo-Boolean
  |      constraint.
  |    * 'assumptions' is updated with a new set of assumptions.
  |
  |________________________________________________________________________________________________@*/
void SWC::encode(Solver *S, vec<Lit> &lits, vec<uint64_t> &coeffs, uint64_t rhs,
                 vec<Lit> &assumptions, int size) {

  // If the rhs is larger than INT32_MAX is not feasible to encode this
  // pseudo-Boolean constraint to CNF.
  if (rhs >= INT32_MAX) {
    printf("c Overflow in the Encoding\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }
  hasEncoding = false;

  vec<Lit> simp_lits;
  vec<uint64_t> simp_coeffs;
  lits.copyTo(simp_lits);
  coeffs.copyTo(simp_coeffs);

  lits.clear();
  coeffs.clear();

  // Add literals from the fixed literals if their coeff is smaller than rhs.
  vec<Lit> simp_unit_lits;
  vec<uint64_t> simp_unit_coeffs;
  unit_lits.copyTo(simp_unit_lits);
  unit_coeffs.copyTo(simp_unit_coeffs);

  unit_lits.clear();
  unit_coeffs.clear();

  for (int i = 0; i < simp_unit_lits.size(); i++) {
    if (simp_unit_coeffs[i] >= INT32_MAX) {
      printf("c Overflow in the Encoding\n");
      printf("s UNKNOWN\n");
      exit(_ERROR_);
    }

    if (simp_unit_coeffs[i] <= (unsigned)rhs) {
      lits.push(simp_unit_lits[i]);
      coeffs.push(simp_unit_coeffs[i]);
    } else {
      unit_lits.push(simp_unit_lits[i]);
      unit_coeffs.push(simp_unit_coeffs[i]);
    }
  }

  // Fix literals that have a coeff larger than rhs.
  for (int i = 0; i < simp_lits.size(); i++) {
    if (simp_coeffs[i] <= (unsigned)rhs) {
      lits.push(simp_lits[i]);
      coeffs.push(simp_coeffs[i]);
    } else {
      unit_lits.push(simp_lits[i]);
      unit_coeffs.push(simp_coeffs[i]);
    }
  }

  if (lits.size() == 1) {
    for (int i = 0; i < unit_lits.size(); i++)
      assumptions.push(~unit_lits[i]);

    unit_lits.push(lits[0]);
    unit_coeffs.push(coeffs[0]);
    return;
  }

  if (lits.size() == 0) {
    for (int i = 0; i < unit_lits.size(); i++)
      assumptions.push(~unit_lits[i]);
    return;
  }

  // Create auxiliary variables.
  int n = lits.size();
  seq_auxiliary_inc = new vec<Lit>[size + 1];

  for (int i = 0; i <= n; i++)
    seq_auxiliary_inc[i].growTo(rhs + 1);

  for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= (int)rhs; ++j) {
      seq_auxiliary_inc[i][j] = mkLit(S->nVars(), false);
      newSATVariable(S);
      nb_variables++;
    }
  }

  Lit blocking = mkLit(S->nVars(), false);
  newSATVariable(S);
  current_lit_blocking = blocking;
  assumptions.push(~blocking);

  for (int i = 1; i <= n; i++) {
    // WARNING: wi is used as int for array indexes but as int64_t
    // for the coeffs. Dangerous if the coeffs are larger than INT32_MAX.
    // Same problem occurs with rhs.
    int64_t wi = coeffs[i - 1];
    // assert(rhs >= wi);

    for (int j = 1; j <= (int)rhs; j++) {
      if (i >= 2 && i <= n && j <= (int)rhs) {
        addBinaryClause(S, ~seq_auxiliary_inc[i - 1][j],
                        seq_auxiliary_inc[i][j]);
        nb_clauses++;
      }
      if (i <= n && j <= wi) {
        addBinaryClause(S, ~lits[i - 1], seq_auxiliary_inc[i][j]);
        nb_clauses++;
      }
      if (i >= 2 && i <= n && j <= (int)(rhs - wi)) {
        addTernaryClause(S, ~seq_auxiliary_inc[i - 1][j], ~lits[i - 1],
                         seq_auxiliary_inc[i][(int)j + (int)wi]);
        nb_clauses++;
      }
    }

    // Encode rhs.
    if (i >= 2) {
      addBinaryClause(S, ~seq_auxiliary_inc[i - 1][(int)rhs + 1 - (int)wi],
                      ~lits[i - 1], blocking);
      nb_clauses++;
    }
  }

  for (int i = 0; i < unit_lits.size(); i++)
    assumptions.push(~unit_lits[i]);

  current_pb_rhs = rhs;
  hasEncoding = true;

  lits.copyTo(lits_inc);
  coeffs.copyTo(coeffs_inc);
}

/*_________________________________________________________________________________________________
  |
  |  update : (S : Solver *) (rhs : int64_t) ->  [void]
  |
  |  Description:
  |
  |     Updates the 'rhs' of an already existent pseudo-Boolean encoding.
  |     This method allows for all learned clauses from previous iterations to
  |     be kept in the next iteration.
  |
  |  Pre-conditions:
  |    * Only one pseudo-Boolean constraint is being used.
  |    * Assumes that 'current_pb_rhs' is larger than -1, i.e. assumes a
  |      cardinality constraint has already been encoded.
  |
  |  Post-conditions:
  |    * 'S' is updated with unit clauses that update the 'rhs' of the
  |      pseudo-Boolean encoding.
  |
  |________________________________________________________________________________________________@*/
void SWC::update(Solver *S, uint64_t rhs) {
  if (rhs >= INT32_MAX) {
    printf("c Overflow in the Encoding\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }

  assert(current_pb_rhs != -1);
  for (int i = rhs; i < current_pb_rhs; i++)
    addUnitClause(S, ~pb_outlits[i]);

  current_pb_rhs = rhs;
}

/*_________________________________________________________________________________________________
  |
  |  update : (S : Solver *) (rhs : int64_t) ( ) ->  [void]
  |
  |  Description:
  |
  |     Incremental update of the SWC encoding.
  |     Encodes the necessary clauses from a_1 x_1 + ... + a_n x_n <= k to
  |     a_1 x_1 + ... + a_n x_n <= k', with k' < k.
  |
  |  Pre-conditions:
  |    * Only one pseudo-Boolean constraint is being used.
  |    * Assumes that 'current_pb_rhs' is larger than -1, i.e. assumes a
  |      cardinality constraint has already been encoded.
  |
  |  Post-conditions:
  |    * 'S' is updated with unit clauses that update the 'rhs' of the
  |      pseudo-Boolean encoding.
  |    * 'assumptions' is updated with a new set of assumptions.
  |
  |________________________________________________________________________________________________@*/
void SWC::update(Solver *S, uint64_t rhs, vec<Lit> &assumptions) {
  if (rhs >= INT32_MAX) {
    printf("c Overflow in the Encoding\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }

  // Disable previous rhs.
  if (current_lit_blocking != lit_Undef)
    addUnitClause(S, current_lit_blocking);

  // Update coeffs.
  // create auxiliary variables.
  int n = lits_inc.size();
  int offset = current_pb_rhs + 1;
  assert((unsigned)current_pb_rhs < rhs);

  for (int i = 1; i <= n; i++) {
    for (int j = offset; j <= (int)rhs; j++)
      seq_auxiliary_inc[i].push(lit_Undef);
  }

  for (int i = 1; i <= n; ++i) {
    for (int j = offset; j <= (int)rhs; ++j) {
      assert(seq_auxiliary_inc[i].size() > j);
      seq_auxiliary_inc[i][j] = mkLit(S->nVars(), false);
      newSATVariable(S);
      nb_variables++;
    }
  }

  for (int i = 1; i < lits_inc.size(); i++)
    assert(seq_auxiliary_inc[i].size() == (int)rhs + 1);

  // blocking literal
  current_lit_blocking = mkLit(S->nVars(), false);
  newSATVariable(S);

  for (int i = 1; i <= n; i++) {
    // WARNING: wi is used as int for array indexes but as int64_t
    // for the coeffs. Dangerous if the coeffs are larger than INT32_MAX.
    // Same problem occurs with rhs.
    uint64_t wi = coeffs_inc[i - 1];
    assert(wi > 0);
    assert(rhs >= wi);

    for (int j = 1; j <= (int)rhs; j++) {
      if (i >= 2 && i <= n && j <= (int)rhs && j >= offset) {
        assert(seq_auxiliary_inc[i].size() > j);
        addBinaryClause(S, ~seq_auxiliary_inc[i - 1][j],
                        seq_auxiliary_inc[i][j]);
        nb_clauses++;
      }
      if (i >= 2 && i <= n && j <= (int)(rhs - wi) && j >= offset - (int)wi) {
        addTernaryClause(S, ~seq_auxiliary_inc[i - 1][j], ~lits_inc[i - 1],
                         seq_auxiliary_inc[i][j + (int)wi]);
        nb_clauses++;
      }
    }

    // encode rhs
    if (i >= 2) {
      assert(i - 1 > 0 &&
             seq_auxiliary_inc[i - 1].size() > (int)rhs + 1 - (int)wi);
      assert(rhs + 1 - wi > 0);
      assert(i - 1 > 0 && i - 1 < lits_inc.size());
      addBinaryClause(S, ~seq_auxiliary_inc[i - 1][(int)rhs + 1 - (int)wi],
                      ~lits_inc[i - 1], current_lit_blocking);
      nb_clauses++;
    }
  }

  current_pb_rhs = rhs;
}

/*_________________________________________________________________________________________________
  |
  |  join : (S : Solver *) (lits : vec<Lit>&)
  |         (coeffs : vec<Lit>&) (assumptions: vec<Lit>&) ->  [void]
  |
  |  Description:
  |
  |     Joins two PB constraints. Given a_1 x_1 + ... + a_n x_n <= k and
  |     b_1 y_1 + ... + b_m y_m <= k, encodes the necessary clauses to
  |     express: a_1 x_1 + ... + a_n x_n + b_1 y_1 + ... + b_m y_m <= k.
  |
  |  Pre-conditions:
  |    * Only one pseudo-Boolean constraint is being used.
  |    * Assumes that 'current_pb_rhs' is larger than -1, i.e. assumes a
  |      cardinality constraint has already been encoded.
  |
  |  Post-conditions:
  |    * 'S' is updated with unit clauses that update the 'rhs' of the
  |      pseudo-Boolean encoding.
  |    * 'assumptions' is updated.
  |
  |________________________________________________________________________________________________@*/
void SWC::join(Solver *S, vec<Lit> &lits, vec<uint64_t> &coeffs,
               vec<Lit> &assumptions) {

  assert(current_lit_blocking != lit_Undef);
  int64_t rhs = current_pb_rhs;

  // If the rhs is larger than INT32_MAX is not feasible to encode this
  // pseudo-Boolean constraint to CNF.
  if (rhs >= INT32_MAX) {
    printf("c Overflow in the Encoding\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }

  // Add literals from the fixed literals if their coeff is smaller than rhs.
  vec<Lit> simp_unit_lits;
  vec<uint64_t> simp_unit_coeffs;
  unit_lits.copyTo(simp_unit_lits);
  unit_coeffs.copyTo(simp_unit_coeffs);

  unit_lits.clear();
  unit_coeffs.clear();

  int lhs_join = lits_inc.size();

  for (int i = 0; i < simp_unit_lits.size(); i++) {
    if (simp_unit_coeffs[i] <= (unsigned)rhs) {
      lits_inc.push(simp_unit_lits[i]);
      coeffs_inc.push(simp_unit_coeffs[i]);
    } else {
      unit_lits.push(simp_unit_lits[i]);
      unit_coeffs.push(simp_unit_coeffs[i]);
    }
  }

  // Fix literals that have a coeff larger than rhs.
  for (int i = 0; i < lits.size(); i++) {
    if (coeffs[i] <= (unsigned)rhs) {
      lits_inc.push(lits[i]);
      coeffs_inc.push(coeffs[i]);
    } else {
      unit_lits.push(lits[i]);
      unit_coeffs.push(coeffs[i]);
    }
  }

  // No literals have been added.
  if (lits_inc.size() == lhs_join)
    return;

  // Create auxiliary variables.
  int n = lits_inc.size();
  int offset = lhs_join;
  assert(seq_auxiliary_inc[lhs_join].size() > 0);

  for (int i = offset + 1; i <= n; i++) {
    assert(seq_auxiliary_inc[i].size() == 0);
    seq_auxiliary_inc[i].growTo(rhs + 1);
  }

  for (int i = offset + 1; i <= n; ++i) {
    for (int j = 1; j <= rhs; ++j) {
      seq_auxiliary_inc[i][j] = mkLit(S->nVars(), false);
      newSATVariable(S);
      nb_variables++;
    }
  }

  for (int i = 1; i <= n; i++)
    assert(seq_auxiliary_inc[i].size() == rhs + 1);

  for (int i = offset; i <= n; i++) {
    // WARNING: wi is used as int for array indexes but as int64_t
    // for the coeffs. Dangerous if the coeffs are larger than INT32_MAX.
    // Same problem occurs with rhs.
    int64_t wi = coeffs_inc[i - 1];
    assert(wi > 0);
    assert(wi <= rhs);

    for (int j = 1; j <= rhs; j++) {
      if (i <= n) {
        assert(seq_auxiliary_inc[i].size() > j);
        assert(seq_auxiliary_inc[i - 1].size() > j);
        addBinaryClause(S, ~seq_auxiliary_inc[i - 1][j],
                        seq_auxiliary_inc[i][j]);
        nb_clauses++;
      }

      if (i <= n && j <= wi) {
        assert(seq_auxiliary_inc[i].size() > j);
        assert(i - 1 < lits_inc.size() && i - 1 >= 0);
        addBinaryClause(S, ~lits_inc[i - 1], seq_auxiliary_inc[i][j]);
        nb_clauses++;
      }
      if (i <= n && j <= rhs - wi) {
        addTernaryClause(S, ~seq_auxiliary_inc[i - 1][j], ~lits_inc[i - 1],
                         seq_auxiliary_inc[i][j + (int)wi]);
        nb_clauses++;
      }
    }

    // encode rhs
    if (i > offset) {
      assert(rhs + 1 - wi >= 0);
      assert(seq_auxiliary_inc[i - 1].size() > rhs + 1 - wi);
      assert(i - 1 < lits_inc.size() && i - 1 >= 0);
      addBinaryClause(S, ~seq_auxiliary_inc[i - 1][(int)rhs + 1 - (int)wi],
                      ~lits_inc[i - 1], current_lit_blocking);
      nb_clauses++;
    }
  }
}
