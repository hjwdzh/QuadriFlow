/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *
 * MiniSat,  Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
 *           Copyright (c) 2007-2010, Niklas Sorensson
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

#ifndef MaxSAT_h
#define MaxSAT_h

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "MaxSATFormula.h"
#include "MaxTypes.h"
#include "utils/System.h"
#include <algorithm>
#include <map>
#include <set>
#include <utility>
#include <vector>

using NSPACE::vec;
using NSPACE::Lit;
using NSPACE::lit_Undef;
using NSPACE::mkLit;
using NSPACE::lbool;
using NSPACE::Solver;
using NSPACE::cpuTime;

namespace openwbo {

class MaxSAT {

public:
  MaxSAT(MaxSATFormula *mx) {
    maxsat_formula = mx;

    // 'ubCost' will be set to the sum of the weights of soft clauses
    //  during the parsing of the MaxSAT formula.
    ubCost = 0;
    lbCost = 0;

    off_set = 0;

    // Statistics
    nbSymmetryClauses = 0;
    nbCores = 0;
    nbSatisfiable = 0;
    sumSizeCores = 0;

    print_model = false;
  }

  MaxSAT() {
    maxsat_formula = NULL;

    // 'ubCost' will be set to the sum of the weights of soft clauses
    //  during the parsing of the MaxSAT formula.
    ubCost = 0;
    lbCost = 0;

    off_set = 0;

    // Statistics
    nbSymmetryClauses = 0;
    nbCores = 0;
    nbSatisfiable = 0;
    sumSizeCores = 0;

    print_model = false;
  }

  virtual ~MaxSAT() {
    if (maxsat_formula != NULL)
      delete maxsat_formula;
  }

  void setInitialTime(double initial); // Set initial time.

  // Print configuration of the MaxSAT solver.
  // virtual void printConfiguration();
  void printConfiguration();

  // Encoding information.
  void print_AMO_configuration(int encoding);
  void print_PB_configuration(int encoding);
  void print_Card_configuration(int encoding);

  // Incremental information.
  void print_Incremental_configuration(int incremental);

  virtual void search();      // MaxSAT search.
  void printAnswer(int type); // Print the answer.

  // Tests if a MaxSAT formula has a lexicographical optimization criterion.
  bool isBMO(bool cache = true);

  void loadFormula(MaxSATFormula *maxsat) {
    maxsat_formula = maxsat;
    maxsat_formula->setInitialVars(maxsat_formula->nVars());

    if (maxsat_formula->getObjFunction() != NULL) {
      off_set = maxsat_formula->getObjFunction()->_const;
      maxsat_formula->convertPBtoMaxSAT();
    }

    ubCost = maxsat_formula->getSumWeights();
  }

  void blockModel(Solver *solver);

  // Get bounds methods
  uint64_t getUB();
  std::pair<uint64_t, int> getLB();

  Soft &getSoftClause(int i) { return maxsat_formula->getSoftClause(i); }
  Hard &getHardClause(int i) { return maxsat_formula->getHardClause(i); }
  Lit getAssumptionLit(int soft) {
    return maxsat_formula->getSoftClause(soft).assumption_var;
  }
  Lit getRelaxationLit(int soft, int i = 0) {
    return maxsat_formula->getSoftClause(soft).relaxation_vars[i];
  }

  int64_t getOffSet() { return off_set; }

  MaxSATFormula *getMaxSATFormula() { return maxsat_formula; }

  void setPrintModel(bool model) { print_model = model; }
  bool getPrintModel() { return print_model; }

protected:
  // Interface with the SAT solver
  //
  Solver *newSATSolver(); // Creates a SAT solver.
  // Solves the formula that is currently loaded in the SAT solver.
  lbool searchSATSolver(Solver *S, vec<Lit> &assumptions, bool pre = false);
  lbool searchSATSolver(Solver *S, bool pre = false);

  void newSATVariable(Solver *S); // Creates a new variable in the SAT solver.

  // Properties of the MaxSAT formula
  //
  vec<lbool> model; // Stores the best satisfying model.

  // Statistics
  //
  int nbCores;           // Number of cores.
  int nbSymmetryClauses; // Number of symmetry clauses.
  uint64_t sumSizeCores; // Sum of the sizes of cores.
  int nbSatisfiable;     // Number of satisfiable calls.

  // Bound values
  //
  uint64_t ubCost; // Upper bound value.
  uint64_t lbCost; // Lower bound value.
  int64_t off_set; // Offset of the objective function for PB solving.

  MaxSATFormula *maxsat_formula;

  // Others
  // int currentWeight;  // Initialized to the maximum weight of soft clauses.
  double initialTime; // Initial time.
  int verbosity;      // Controls the verbosity of the solver.
  bool print_model;   // Controls if the model is printed at the end.

  // Different weights that corresponds to each function in the BMO algorithm.
  std::vector<uint64_t> orderWeights;

  // Utils for model management
  //
  void saveModel(vec<lbool> &currentModel); // Saves a Model.
  // Compute the cost of a model.
  uint64_t computeCostModel(vec<lbool> &currentModel,
                            uint64_t weight = UINT64_MAX);

  // Utils for printing
  //
  void printModel(); // Print the best satisfying model.
  void printStats(); // Print search statistics.

  // Greater than comparator.
  bool static greaterThan(uint64_t i, uint64_t j) { return (i > j); }
};
} // namespace openwbo

#endif
