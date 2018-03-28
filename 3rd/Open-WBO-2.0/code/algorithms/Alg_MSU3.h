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

#ifndef Alg_MSU3_h
#define Alg_MSU3_h

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "../Encoder.h"
#include "../MaxSAT.h"
#include <algorithm>
#include <map>
#include <set>

namespace openwbo {

//=================================================================================================
class MSU3 : public MaxSAT {

public:
  MSU3(int verb = _VERBOSITY_MINIMAL_) {
    solver = NULL;
    verbosity = verb;
    incremental_strategy = _INCREMENTAL_ITERATIVE_;
    encoding = _CARD_TOTALIZER_;
    encoder.setCardEncoding(encoding);
  }
  ~MSU3() {
    if (solver != NULL)
      delete solver;
  }

  void search(); // MSU3 search.

  // Print solver configuration.
  void printConfiguration() {

    printf("c ==========================================[ Solver Settings "
           "]============================================\n");
    printf("c |                                                                "
           "                                       |\n");

    print_MSU3_configuration();
    print_Card_configuration(_CARD_TOTALIZER_);
  }

protected:
  // Print MSU3 configuration.
  void print_MSU3_configuration();

  // Rebuild MaxSAT solver
  //
  Solver *rebuildSolver(); // Rebuild MaxSAT solver.

  void MSU3_none();      // Non-incremental MSU3.
  void MSU3_blocking();  // Incremental Blocking MSU3.
  void MSU3_weakening(); // Incremental Weakening MSU3.
  void MSU3_iterative(); // Incremental Iterative Encoding MSU3.

  // Other
  void initRelaxation(); // Relaxes soft clauses.

  Solver *solver;  // SAT Solver used as a black box.
  Encoder encoder; // Interface for the encoder of constraints to CNF.

  // Controls the incremental strategy used by MSU3 algorithms.
  int incremental_strategy;
  // Controls the cardinality encoding used by MSU3 algorithms.
  int encoding;

  // Literals to be used in the constraint that excludes models.
  vec<Lit> objFunction;
  vec<int> coeffs; // Coefficients of the literals that are used in the
                   // constraint that excludes models.

  std::map<Lit, int> coreMapping; // Mapping between the assumption literal and
                                  // the respective soft clause.

  // Soft clauses that are currently in the MaxSAT formula.
  vec<bool> activeSoft;
};
} // namespace openwbo

#endif
