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

#ifndef Alg_WBO_h
#define Alg_WBO_h

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "../Encoder.h"
#include "../MaxSAT.h"
#include "../MaxTypes.h"
#include "utils/System.h"
#include <map>
#include <set>
#include <utility>

namespace openwbo {

class WBO : public MaxSAT {

public:
  // NOTE: currently the encoding is not set as an input parameter.
  WBO(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NONE_,
      bool symmetry = true, int limit = INT32_MAX) {
    solver = NULL;
    verbosity = verb;

    nbCurrentSoft = 0;
    weightStrategy = weight;

    symmetryStrategy = symmetry;
    symmetryBreakingLimit = limit;
  }

  ~WBO() {
    if (solver != NULL)
      delete solver;
  }

  void search(); // WBO search.

protected:
  // Rebuild MaxSAT solver
  //
  // Rebuild MaxSAT solver with weight-based strategy.
  Solver *rebuildWeightSolver(int strategy);
  Solver *rebuildSolver();     // Rebuild MaxSAT solver.
  Solver *rebuildHardSolver(); // Rebuild MaxSAT solver with only hard clauses.
  void updateCurrentWeight(int strategy); // Updates 'currentWeight'.
  uint64_t
  findNextWeight(uint64_t weight); // Finds the next weight for 'currentWeight'.
  uint64_t
  findNextWeightDiversity(uint64_t weight); // Finds the next weight for
                                            // 'currentWeight' using diversify
                                            // heuristic.

  // Utils for core management
  //
  void encodeEO(vec<Lit> &lits); // Encodes exactly one constraint.
  void relaxCore(vec<Lit> &conflict, uint64_t weightCore,
                 vec<Lit> &assumps);            // Relaxes a core.
  uint64_t computeCostCore(vec<Lit> &conflict); // Computes the cost of a core.

  // Symmetry breaking methods
  //
  void initSymmetry();     // Initializes symmetry breaking data structures.
  void symmetryLog(int p); // Logs symmetry breaking information.
  void symmetryBreaking(); // Adds symmetry breaking clauses.

  // WBO search
  //
  void unsatSearch();  // Search using only hard clauses.
  void weightSearch(); // Search using weight-based methods.
  void normalSearch(); // Original WBO search.

  // Other
  // Initializes assumptions and core extraction.
  void initAssumptions(vec<Lit> &assumps);

  // SAT solver
  Solver *solver;  // SAT solver used as a black box.
  Encoder encoder; // Interface for the encoder of constraints to CNF.

  // Variables used  in 'weightSearch'
  //
  int nbCurrentSoft;  // Current number of soft clauses used by the MaxSAT
                      // solver.
  int weightStrategy; // Weight strategy to be used in 'weightSearch'.

  // Core extraction
  //
  std::map<Lit, int> coreMapping; // Maps the assumption literal to the number
                                  // of the soft clause.
  vec<Lit> assumptions; // Stores the assumptions to be used in the extraction
                        // of the core.

  // Symmetry breaking
  //
  bool symmetryStrategy; // Symmetry breaking strategy.
  vec<int>
      indexSoftCore; // Indexes of soft clauses that appear in the current core.
  // Maps the soft clause with the cores where they appears.
  vec<vec<int>> softMapping;
  vec<vec<Lit>> relaxationMapping; // Maps the relaxation variables with the
                                   // soft clause where they appear.

  typedef std::pair<int, int> symmetryClause; // Stores binary symmetry clauses.
  std::set<symmetryClause> duplicatedSymmetryClauses; // Set of binary symmetry
                                                      // clauses (prevents
                                                      // duplication).
  int symmetryBreakingLimit; // Limit on the number of symmetry clauses.
};
} // namespace openwbo

#endif
