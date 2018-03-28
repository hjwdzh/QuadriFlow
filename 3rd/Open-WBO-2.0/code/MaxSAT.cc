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

#include "MaxSAT.h"

using namespace openwbo;

/************************************************************************************************
 //
 // Public methods
 //
 ************************************************************************************************/

void MaxSAT::search() {
  printf("Error: Invalid MaxSAT algoritm.\n");
  exit(_ERROR_);
}

void MaxSAT::setInitialTime(double initial) {
  initialTime = initial;
} // Sets the initial time.

/************************************************************************************************
 //
 // SAT solver interface
 //
 ************************************************************************************************/

// Creates an empty SAT Solver.
Solver *MaxSAT::newSATSolver() {

#ifdef SIMP
  NSPACE::SimpSolver *S = new NSPACE::SimpSolver();
#else
  Solver *S = new Solver();
#endif

  return (Solver *)S;
}

// Creates a new variable in the SAT solver.
void MaxSAT::newSATVariable(Solver *S) {

#ifdef SIMP
  ((NSPACE::SimpSolver *)S)->newVar();
#else
  S->newVar();
#endif
}

// Solve the formula that is currently loaded in the SAT solver with a set of
// assumptions and with the option to use preprocessing for 'simp'.
lbool MaxSAT::searchSATSolver(Solver *S, vec<Lit> &assumptions, bool pre) {

// Currently preprocessing is disabled by default.
// Variable elimination cannot be done on relaxation variables nor on variables
// that belong to soft clauses. To preprocessing to be used those variables
// should be frozen.

#ifdef SIMP
  lbool res = ((NSPACE::SimpSolver *)S)->solveLimited(assumptions, pre);
#else
  lbool res = S->solveLimited(assumptions);
#endif

  return res;
}

// Solve the formula without assumptions.
lbool MaxSAT::searchSATSolver(Solver *S, bool pre) {
  vec<Lit> dummy; // Empty set of assumptions.
  return searchSATSolver(S, dummy, pre);
}

/************************************************************************************************
 //
 // Utils for model management
 //
 ************************************************************************************************/

/*_________________________________________________________________________________________________
  |
  |  saveModel : (currentModel : vec<lbool>&)  ->  [void]
  |
  |  Description:
  |
  |    Saves the current model found by the SAT solver.
  |
  |  Pre-conditions:
  |    * Assumes that 'nbInitialVariables' has been initialized.
  |    * Assumes that 'currentModel' is not empty.
  |
  |  Post-conditions:
  |    * 'model' is updated to the current model.
  |
  |________________________________________________________________________________________________@*/
void MaxSAT::saveModel(vec<lbool> &currentModel) {
  assert(maxsat_formula->nInitialVars() != 0);
  assert(currentModel.size() != 0);

  model.clear();
  // Only store the value of the variables that belong to the
  // original MaxSAT formula.
  for (int i = 0; i < maxsat_formula->nInitialVars(); i++)
    model.push(currentModel[i]);
}

/*_________________________________________________________________________________________________
  |
  |  computeCostModel : (currentModel : vec<lbool>&) (weight : int) ->
  |                     [uint64_t]
  |
  |  Description:
  |
  |    Computes the cost of 'currentModel'. The cost of a model is the sum of
  |    the weights of the unsatisfied soft clauses.
  |    If a weight is specified, then it only considers the sum of the weights
  |    of the unsatisfied soft clauses with the specified weight.
  |
  |  Pre-conditions:
  |    * Assumes that 'currentModel' is not empty.
  |
  |________________________________________________________________________________________________@*/
uint64_t MaxSAT::computeCostModel(vec<lbool> &currentModel, uint64_t weight) {

  assert(currentModel.size() != 0);
  uint64_t currentCost = 0;

  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    bool unsatisfied = true;
    for (int j = 0; j < maxsat_formula->getSoftClause(i).clause.size(); j++) {

      if (weight != UINT64_MAX &&
          maxsat_formula->getSoftClause(i).weight != weight) {
        unsatisfied = false;
        continue;
      }

      assert(var(maxsat_formula->getSoftClause(i).clause[j]) <
             currentModel.size());
      if ((sign(maxsat_formula->getSoftClause(i).clause[j]) &&
           currentModel[var(maxsat_formula->getSoftClause(i).clause[j])] ==
               l_False) ||
          (!sign(maxsat_formula->getSoftClause(i).clause[j]) &&
           currentModel[var(maxsat_formula->getSoftClause(i).clause[j])] ==
               l_True)) {
        unsatisfied = false;
        break;
      }
    }

    if (unsatisfied) {
      currentCost += maxsat_formula->getSoftClause(i).weight;
    }
  }

  return currentCost;
}

/*_________________________________________________________________________________________________
  |
  |  isBMO : (cache : bool)  ->  [void]
  |
  |  Description:
  |
  |    Tests if the MaxSAT formula has lexicographical optimization criterion.
  |
  |  For further details see:
  |    * Joao Marques-Silva, Josep Argelich, Ana Graça, Inês Lynce: Boolean
  |      lexicographic optimization: algorithms & applications. Ann. Math.
  |      Artif. Intell. 62(3-4): 317-343 (2011)
  |
  |  Post-conditions:
  |    * 'orderWeights' is updated with the weights in lexicographical order if
  |      'cache' is true.
  |
  |________________________________________________________________________________________________@*/
bool MaxSAT::isBMO(bool cache) {
  assert(orderWeights.size() == 0);
  bool bmo = true;
  std::set<uint64_t> partitionWeights;
  std::map<uint64_t, uint64_t> nbPartitionWeights;

  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    partitionWeights.insert(maxsat_formula->getSoftClause(i).weight);
    nbPartitionWeights[maxsat_formula->getSoftClause(i).weight]++;
  }

  for (std::set<uint64_t>::iterator iter = partitionWeights.begin();
       iter != partitionWeights.end(); ++iter) {
    orderWeights.push_back(*iter);
  }

  std::sort(orderWeights.begin(), orderWeights.end(), greaterThan);

  uint64_t totalWeights = 0;
  for (int i = 0; i < (int)orderWeights.size(); i++)
    totalWeights += orderWeights[i] * nbPartitionWeights[orderWeights[i]];

  for (int i = 0; i < (int)orderWeights.size(); i++) {
    totalWeights -= orderWeights[i] * nbPartitionWeights[orderWeights[i]];
    if (orderWeights[i] < totalWeights) {
      bmo = false;
      break;
    }
  }

  if (!cache)
    orderWeights.clear();

  return bmo;
}

/************************************************************************************************
 //
 // Utils for printing
 //
 ************************************************************************************************/

// Prints information regarding the AMO encoding.
void MaxSAT::print_AMO_configuration(int encoding) {
  switch (encoding) {
  case _AMO_LADDER_:
    printf("c |  AMO Encoding:         %12s                      "
           "                                             |\n",
           "Ladder");
    break;

  default:
    printf("c Error: Invalid AMO encoding.\n");
    printf("s UNKNOWN\n");
    break;
  }
}

// Prints information regarding the PB encoding.
void MaxSAT::print_PB_configuration(int encoding) {
  switch (encoding) {
  case _PB_SWC_:
    printf("c |  PB Encoding:         %13s                        "
           "                                           |\n",
           "SWC");
    break;

  case _PB_GTE_:
    printf("c |  PB Encoding:         %13s                        "
           "                                           |\n",
           "GTE");
    break;

  default:
    printf("c Error: Invalid PB encoding.\n");
    printf("s UNKNOWN\n");
    break;
  }
}

// Prints information regarding the cardinality encoding.
void MaxSAT::print_Card_configuration(int encoding) {
  switch (encoding) {
  case _CARD_CNETWORKS_:
    printf("c |  Cardinality Encoding: %12s                                "
           "                                   |\n",
           "CNetworks");
    break;

  case _CARD_TOTALIZER_:
    printf("c |  Cardinality Encoding: %12s                                "
           "                                   |\n",
           "Totalizer");
    break;

  case _CARD_MTOTALIZER_:
    printf("c |  Cardinality Encoding:    %19s                             "
           "                            |\n",
           "Modulo Totalizer");
    break;

  default:
    printf("c Error: Invalid cardinality encoding.\n");
    printf("s UNKNOWN\n");
    break;
  }
}

void MaxSAT::blockModel(Solver *solver) {
  assert(model.size() != 0);

  vec<Lit> blocking;

  printf("v ");
  for (int i = 0; i < model.size(); i++) {
    indexMap::const_iterator iter = maxsat_formula->getIndexToName().find(i);
    if (iter != maxsat_formula->getIndexToName().end()) {
      if (model[i] == l_False)
        printf("-");
      printf("%s ", iter->second.c_str());
    }
  }
  printf("\n");

  for (int i = 0; i < model.size(); i++) {
    blocking.push((model[i] == l_True) ? ~mkLit(i) : mkLit(i));
  }

  solver->addClause(blocking);
}

// Prints the best satisfying model. Assumes that 'model' is not empty.
void MaxSAT::printModel() {

  assert(model.size() != 0);

  if (maxsat_formula->getFormat() == _FORMAT_PB_) {

    printf("v ");
    for (int i = 0; i < model.size(); i++) {
      indexMap::const_iterator iter = maxsat_formula->getIndexToName().find(i);
      if (iter != maxsat_formula->getIndexToName().end()) {
        if (model[i] == l_False)
          printf("-");
        printf("%s ", iter->second.c_str());
      }
    }
    printf("\n");

    // printf("v ");
    // for (int i = 0; i < model.size(); i++) {
    //   indexMap::const_iterator iter =
    //   maxsat_formula->getIndexToName().find(i); if (iter !=
    //   maxsat_formula->getIndexToName().end()) {
    //     if (model[i] == l_False) printf("+1 %s = 0
    //     ;\n",iter->second.c_str()); else printf("+1 %s = 1
    //     ;\n",iter->second.c_str());
    //   }
    // }

  } else {

    printf("v ");
    for (int i = 0; i < model.size(); i++) {
      if (model[i] == l_True)
        printf("%d ", i + 1);
      else
        printf("%d ", -(i + 1));
    }
    printf("\n");
  }
}

// Prints search statistics.
void MaxSAT::printStats() {
  double totalTime = cpuTime();
  float avgCoreSize = 0;
  if (nbCores != 0)
    avgCoreSize = (float)sumSizeCores / nbCores;

  printf("c\n");
  if (model.size() == 0)
    printf("c  Best solution:          %12s\n", "-");
  else
    printf("c  Best solution:          %12" PRIu64 "\n", ubCost);
  printf("c  Total time:             %12.2f s\n", totalTime - initialTime);
  printf("c  Nb SAT calls:           %12d\n", nbSatisfiable);
  printf("c  Nb UNSAT calls:         %12d\n", nbCores);
  printf("c  Average core size:      %12.2f\n", avgCoreSize);
  printf("c  Nb symmetry clauses:    %12d\n", nbSymmetryClauses);
  printf("c\n");
}

// Prints the corresponding answer.
void MaxSAT::printAnswer(int type) {
  if (verbosity > 0)
    printStats();

  if (type == _UNKNOWN_ && model.size() > 0)
    type = _SATISFIABLE_;

  switch (type) {
  case _SATISFIABLE_:
    printf("s SATISFIABLE\n");
    if (print_model)
      printModel();
    break;
  case _OPTIMUM_:
    printf("s OPTIMUM FOUND\n");
    if (print_model)
      printModel();
    break;
  case _UNSATISFIABLE_:
    printf("s UNSATISFIABLE\n");
    break;
  case _UNKNOWN_:
    printf("s UNKNOWN\n");
    break;
  default:
    printf("c Error: Invalid answer type.\n");
  }
}

uint64_t MaxSAT::getUB() {
  // only works for partial MaxSAT currently
  Solver *solver = newSATSolver();

  vec<Lit> relaxation_vars;
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    Lit p = mkLit(maxsat_formula->nVars() + i, false);
    relaxation_vars.push(p);
  }

  for (int i = 0; i < maxsat_formula->nVars() + maxsat_formula->nSoft(); i++)
    newSATVariable(solver);

  for (int i = 0; i < maxsat_formula->nHard(); i++)
    solver->addClause(maxsat_formula->getHardClause(i).clause);

  vec<Lit> clause;
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    clause.clear();
    maxsat_formula->getSoftClause(i).clause.copyTo(clause);

    for (int j = 0; j < maxsat_formula->getSoftClause(i).relaxation_vars.size();
         j++)
      clause.push(maxsat_formula->getSoftClause(i).relaxation_vars[j]);

    clause.push(relaxation_vars[i]);

    solver->addClause(clause);
  }

  int limit = 1000;
  solver->setConfBudget(limit);

  vec<Lit> dummy;
  lbool res = searchSATSolver(solver, dummy);
  if (res == l_True) {
    uint64_t ub = computeCostModel(solver->model);
    return ub;
  } else if (res == l_False) {
    printAnswer(_UNSATISFIABLE_);
    exit(_UNSATISFIABLE_);
  }

  return maxsat_formula->nSoft();
}

std::pair<uint64_t, int> MaxSAT::getLB() {
  // only works for partial MaxSAT currently
  Solver *solver = newSATSolver();

  vec<Lit> relaxation_vars;
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    Lit p = mkLit(maxsat_formula->nVars() + i, false);
    relaxation_vars.push(p);
  }

  for (int i = 0; i < maxsat_formula->nVars() + maxsat_formula->nSoft(); i++)
    newSATVariable(solver);

  for (int i = 0; i < maxsat_formula->nHard(); i++)
    solver->addClause(maxsat_formula->getHardClause(i).clause);

  vec<Lit> clause;
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    clause.clear();
    maxsat_formula->getSoftClause(i).clause.copyTo(clause);

    clause.push(relaxation_vars[i]);

    solver->addClause(clause);
  }

  std::map<Lit, int> core; // Mapping between the assumption literal and
                           // the respective soft clause.

  for (int i = 0; i < maxsat_formula->nSoft(); i++)
    core[relaxation_vars[i]] = i;

  int limit = 1000;
  lbool res = l_False;
  uint64_t lb = 0;

  vec<bool> active;
  active.growTo(relaxation_vars.size(), false);
  vec<Lit> assumptions;
  for (int i = 0; i < relaxation_vars.size(); i++) {
    if (!active[i]) {
      assumptions.push(~relaxation_vars[i]);
    }
  }

  while (res == l_False) {
    solver->setConfBudget(limit);
    res = searchSATSolver(solver, assumptions);
    if (res == l_False) {

      for (int i = 0; i < solver->conflict.size(); i++) {
        Lit p = solver->conflict[i];
        if (core.find(p) != core.end()) {
          assert(!active[core[p]]);
          active[core[p]] = true;
        }
      }

      assumptions.clear();
      for (int i = 0; i < relaxation_vars.size(); i++) {
        if (!active[i]) {
          assumptions.push(~relaxation_vars[i]);
        }
      }
      lb++;
    }
  }

  int nb_relaxed = 0;
  for (int i = 0; i < relaxation_vars.size(); i++) {
    if (active[i])
      nb_relaxed++;
  }

  return std::make_pair(lb, nb_relaxed);
}
