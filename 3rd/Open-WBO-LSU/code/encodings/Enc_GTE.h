/*!
 * \author Saurabh Joshi - sbjoshi@iith.ac.in
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

#ifndef Enc_GTE_h
#define Enc_GTE_h

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "Encodings.h"
#include "core/SolverTypes.h"
#include <map>
#include <utility>
#include <vector>

namespace openwbo {
struct wlitt {
  Lit lit;
  uint64_t weight;
};
struct less_than_wlitt {
  inline bool operator()(const wlitt &wl1, const wlitt &wl2) {
    return (wl1.weight < wl2.weight);
  }
};
struct wlit_sumt {
  inline uint64_t operator()(const uint64_t &wl1, const wlitt &wl2) {
    return (wl1 + wl2.weight);
  }
};
typedef std::map<uint64_t, Lit> wlit_mapt;
typedef std::vector<wlitt> weightedlitst;
typedef std::pair<uint64_t, Lit> wlit_pairt;
class GTE : public Encodings {

public:
  GTE() {
    // current_pb_rhs = -1; // -1 corresponds to an unitialized value
    current_pb_rhs = 0;
    nb_clauses = 0;
    nb_variables = 0;
  }
  ~GTE() {}

  // Encode constraint.
  void encode(Solver *S, vec<Lit> &lits, vec<uint64_t> &coeffs, uint64_t rhs);

  // Update constraint.
  void update(Solver *S, uint64_t rhs);

  // Returns true if the encoding was built, otherwise returns false;
  bool hasCreatedEncoding() { return hasEncoding; }

protected:
  void printLit(Lit l) { printf("%s%d\n", sign(l) ? "-" : "", var(l) + 1); }

  bool encodeLeq(uint64_t k, Solver *S, const weightedlitst &iliterals,
                 wlit_mapt &oliterals);
  Lit getNewLit(Solver *S);
  Lit get_var(Solver *S, wlit_mapt &oliterals, uint64_t weight);
  vec<Lit> pb_outlits; // Stores the outputs of the pseudo-Boolean constraint
                       // encoding for incremental solving.
  uint64_t current_pb_rhs; // Stores the current value of the rhs of the
                           // pseudo-Boolean constraint.

  // Stores unit lits. Used for lits that have a coeff larger than rhs.
  wlit_mapt pb_oliterals;
  vec<Lit> unit_lits;
  vec<uint64_t> unit_coeffs;

  // Number of variables and clauses for statistics.
  int nb_variables;
  int nb_clauses;
};

} // namespace openwbo

#endif
