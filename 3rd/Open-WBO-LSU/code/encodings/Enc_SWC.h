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

#ifndef Enc_SWC_h
#define Enc_SWC_h

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "Encodings.h"
#include "core/SolverTypes.h"

namespace openwbo {

class SWC : public Encodings {

public:
  SWC() {
    current_pb_rhs = -1; // -1 corresponds to an unitialized value
    previous_lit_blocking = lit_Undef;
    current_lit_blocking = lit_Undef;
    nb_clauses = 0;
    nb_variables = 0;
  }
  ~SWC() {}

  // Encode constraint.
  void encode(Solver *S, vec<Lit> &lits, vec<uint64_t> &coeffs, uint64_t rhs);
  void encode(Solver *S, vec<Lit> &lits, vec<uint64_t> &coeffs, uint64_t rhs,
              vec<Lit> &assumptions, int size);
  // Update constraint.
  void update(Solver *S, uint64_t rhs);
  void update(Solver *S, uint64_t rhs, vec<Lit> &assumptions);

  // Update assumptions.
  void updateAssumptions(Solver *S, vec<Lit> &assumptions) {
    assumptions.push(~current_lit_blocking);

    for (int i = 0; i < unit_lits.size(); i++)
      assumptions.push(~unit_lits[i]);
  }

  // Join encodings.
  void join(Solver *S, vec<Lit> &lits, vec<uint64_t> &coeffs,
            vec<Lit> &assumptions);

  // Returns true if the encoding was built, otherwise returns false;
  bool hasCreatedEncoding() { return hasEncoding; }

protected:
  vec<Lit> pb_outlits;    // Stores the outputs of the pseudo-Boolean constraint
                          // encoding for incremental solving.
  int64_t current_pb_rhs; // Stores the current value of the rhs of the
                          // pseudo-Boolean constraint.

  Lit previous_lit_blocking; // Blocking literal.
  Lit current_lit_blocking;  // Previous blocking literal.

  // Stores unit lits. Used for lits that have a coeff larger than rhs.
  vec<Lit> unit_lits;
  vec<uint64_t> unit_coeffs;

  // Stores the matrix with the auxiliary variables.
  vec<Lit> *seq_auxiliary_inc;

  // Temporary copy of lits and coeffs for incremental approach.
  vec<Lit> lits_inc;
  vec<uint64_t> coeffs_inc;

  // Number of variables and clauses for statistics.
  int nb_variables;
  int nb_clauses;
};
} // namespace openwbo

#endif
