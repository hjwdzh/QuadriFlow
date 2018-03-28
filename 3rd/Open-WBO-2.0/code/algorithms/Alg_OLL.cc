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

#include "Alg_OLL.h"

using namespace openwbo;

uint64_t OLL::findNextWeight(uint64_t weight,
                             std::set<Lit> &cardinality_assumptions) {

  uint64_t nextWeight = 1;
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    if (maxsat_formula->getSoftClause(i).weight > nextWeight &&
        maxsat_formula->getSoftClause(i).weight < weight)
      nextWeight = maxsat_formula->getSoftClause(i).weight;
  }

  for (std::set<Lit>::iterator it = cardinality_assumptions.begin();
       it != cardinality_assumptions.end(); ++it) {
    assert(boundMapping.find(*it) != boundMapping.end());
    std::pair<std::pair<int, uint64_t>, uint64_t> soft_id = boundMapping[*it];
    if (soft_id.second > nextWeight && soft_id.second < weight)
      nextWeight = soft_id.second;
  }

  return nextWeight;
}

uint64_t OLL::findNextWeightDiversity(uint64_t weight,
                                      std::set<Lit> &cardinality_assumptions) {

  assert(nbSatisfiable > 0); // Assumes that unsatSearch was done before.

  uint64_t nextWeight = weight;
  int nbClauses = 0;
  std::set<uint64_t> nbWeights;
  float alpha = 1.25;

  bool findNext = false;

  for (;;) {
    if (nbSatisfiable > 1 || findNext)
      nextWeight = findNextWeight(nextWeight, cardinality_assumptions);

    nbClauses = 0;
    nbWeights.clear();
    for (int i = 0; i < maxsat_formula->nSoft(); i++) {
      if (maxsat_formula->getSoftClause(i).weight >= nextWeight) {
        nbClauses++;
        nbWeights.insert(maxsat_formula->getSoftClause(i).weight);
      }
    }

    for (std::set<Lit>::iterator it = cardinality_assumptions.begin();
         it != cardinality_assumptions.end(); ++it) {
      assert(boundMapping.find(*it) != boundMapping.end());
      std::pair<std::pair<int, uint64_t>, uint64_t> soft_id = boundMapping[*it];
      if (soft_id.second >= nextWeight) {
        nbClauses++;
        nbWeights.insert(soft_id.second);
      }
    }

    if ((float)nbClauses / nbWeights.size() > alpha ||
        (unsigned)nbClauses ==
            (unsigned)maxsat_formula->nSoft() + cardinality_assumptions.size())
      break;

    if (nbSatisfiable == 1 && !findNext)
      findNext = true;
  }

  return nextWeight;
}

void OLL::unweighted() {
  // printf("unweighted\n");

  // nbInitialVariables = nVars();
  lbool res = l_True;
  initRelaxation();
  solver = rebuildSolver();

  vec<Lit> assumptions;
  vec<Lit> joinObjFunction;
  vec<Lit> currentObjFunction;
  vec<Lit> encodingAssumptions;
  encoder.setIncremental(_INCREMENTAL_ITERATIVE_);

  activeSoft.growTo(maxsat_formula->nSoft(), false);
  for (int i = 0; i < maxsat_formula->nSoft(); i++)
    coreMapping[maxsat_formula->getSoftClause(i).assumption_var] = i;

  std::set<Lit> cardinality_assumptions;
  vec<Encoder *> soft_cardinality;

  for (;;) {

    res = searchSATSolver(solver, assumptions);
    if (res == l_True) {
      nbSatisfiable++;
      uint64_t newCost = computeCostModel(solver->model);
      saveModel(solver->model);
      if (maxsat_formula->getFormat() == _FORMAT_PB_) {
        // optimization problem
        if (maxsat_formula->getObjFunction() != NULL) {
          printf("o %" PRId64 "\n", newCost + off_set);
        }
      } else
        printf("o %" PRId64 "\n", newCost + off_set);

      ubCost = newCost;

      if (nbSatisfiable == 1) {
        if (newCost == 0) {
          if (maxsat_formula->getFormat() == _FORMAT_PB_ &&
              maxsat_formula->getObjFunction() == NULL) {
            printAnswer(_SATISFIABLE_);
            exit(_SATISFIABLE_);
          } else {
            printAnswer(_OPTIMUM_);
            exit(_OPTIMUM_);
          }
        }

        for (int i = 0; i < maxsat_formula->nSoft(); i++)
          assumptions.push(~maxsat_formula->getSoftClause(i).assumption_var);
      } else {
        assert(lbCost == newCost);
        printAnswer(_OPTIMUM_);
        exit(_OPTIMUM_);
      }
    }

    if (res == l_False) {
      lbCost++;
      nbCores++;
      if (verbosity > 0)
        printf("c LB : %-12" PRIu64 "\n", lbCost);

      if (nbSatisfiable == 0) {
        printAnswer(_UNSATISFIABLE_);
        exit(_UNSATISFIABLE_);
      }

      if (lbCost == ubCost) {
        assert(nbSatisfiable > 0);
        if (verbosity > 0)
          printf("c LB = UB\n");
        printAnswer(_OPTIMUM_);
        exit(_OPTIMUM_);
      }

      sumSizeCores += solver->conflict.size();

      vec<Lit> soft_relax;
      vec<Lit> cardinality_relax;

      for (int i = 0; i < solver->conflict.size(); i++) {
        Lit p = solver->conflict[i];
        if (coreMapping.find(p) != coreMapping.end()) {
          assert(!activeSoft[coreMapping[p]]);
          activeSoft[coreMapping[solver->conflict[i]]] = true;
          assert(p ==
                 maxsat_formula->getSoftClause(coreMapping[solver->conflict[i]])
                     .relaxation_vars[0]);
          soft_relax.push(p);
        }

        if (boundMapping.find(p) != boundMapping.end()) {
          std::set<Lit>::iterator it;
          it = cardinality_assumptions.find(p);
          assert(it != cardinality_assumptions.end());
          cardinality_assumptions.erase(it);
          cardinality_relax.push(p);

          // this is a soft cardinality -- bound must be increased
          std::pair<std::pair<int, int>, int> soft_id =
              boundMapping[solver->conflict[i]];
          // increase the bound
          assert(soft_id.first.first < soft_cardinality.size());
          assert(soft_cardinality[soft_id.first.first]->hasCardEncoding());

          joinObjFunction.clear();
          encodingAssumptions.clear();
          soft_cardinality[soft_id.first.first]->incUpdateCardinality(
              solver, joinObjFunction,
              soft_cardinality[soft_id.first.first]->lits(),
              soft_id.first.second + 1, encodingAssumptions);

          // if the bound is the same as the number of lits then no restriction
          // is applied
          if (soft_id.first.second + 1 <
              soft_cardinality[soft_id.first.first]->outputs().size()) {
            assert(soft_cardinality[soft_id.first.first]->outputs().size() >
                   soft_id.first.second + 1);
            Lit out = soft_cardinality[soft_id.first.first]
                          ->outputs()[soft_id.first.second + 1];
            boundMapping[out] = std::make_pair(
                std::make_pair(soft_id.first.first, soft_id.first.second + 1),
                1);
            cardinality_assumptions.insert(out);
          }
        }
      }

      assert(soft_relax.size() + cardinality_relax.size() > 0);

      if (soft_relax.size() == 1 && cardinality_relax.size() == 0) {
        // Unit core
        solver->addClause(soft_relax[0]);
      }

      // assert (soft_relax.size() > 0 || cardinality_relax.size() != 1);

      if (soft_relax.size() + cardinality_relax.size() > 1) {

        vec<Lit> relax_harden;
        soft_relax.copyTo(relax_harden);
        for (int i = 0; i < cardinality_relax.size(); i++)
          relax_harden.push(cardinality_relax[i]);

        /*
                for (int i = 0; i < nSoft(); i++){
                  printf("soft %d, w: %d, r:
           %d\n",i,softClauses[i].weight,var(softClauses[i].relaxationVars[0])+1);
                  for (int j = 0; j < softClauses[i].clause.size(); j++){
                    printf("%d ",var(softClauses[i].clause[j])+1);
                  }
                  printf("\n");
                }

                for (int i = 0; i < relax_harden.size(); i++)
                  printf("+1 x%d ",var(relax_harden[i])+1);
                printf(" <= 1\n");
        */

        Encoder *e = new Encoder();
        e->setIncremental(_INCREMENTAL_ITERATIVE_);
        e->buildCardinality(solver, relax_harden, 1);
        soft_cardinality.push(e);
        assert(e->outputs().size() > 1);

        Lit out = e->outputs()[1];
        boundMapping[out] =
            std::make_pair(std::make_pair(soft_cardinality.size() - 1, 1), 1);
        cardinality_assumptions.insert(out);
      }

      // reset the assumptions
      assumptions.clear();
      int active_soft = 0;
      for (int i = 0; i < maxsat_formula->nSoft(); i++) {
        if (!activeSoft[i])
          assumptions.push(~maxsat_formula->getSoftClause(i).assumption_var);
        else
          active_soft++;
      }

      for (std::set<Lit>::iterator it = cardinality_assumptions.begin();
           it != cardinality_assumptions.end(); ++it) {
        assumptions.push(~(*it));
      }

      if (verbosity > 0) {
        printf("c Relaxed soft clauses %d / %d\n", active_soft,
               maxsat_formula->nSoft());
      }
    }
  }
}

void OLL::weighted() {
  // nbInitialVariables = nVars();
  lbool res = l_True;
  initRelaxation();
  solver = rebuildSolver();

  vec<Lit> assumptions;
  vec<Lit> joinObjFunction;
  vec<Lit> currentObjFunction;
  vec<Lit> encodingAssumptions;
  encoder.setIncremental(_INCREMENTAL_ITERATIVE_);

  activeSoft.growTo(maxsat_formula->nSoft(), false);
  for (int i = 0; i < maxsat_formula->nSoft(); i++)
    coreMapping[maxsat_formula->getSoftClause(i).assumption_var] = i;

  std::set<Lit> cardinality_assumptions;
  vec<Encoder *> soft_cardinality;

  min_weight = maxsat_formula->getMaximumWeight();
  // printf("current weight %d\n",maxsat_formula->getMaximumWeight());

  for (;;) {

    res = searchSATSolver(solver, assumptions);
    if (res == l_True) {
      nbSatisfiable++;
      uint64_t newCost = computeCostModel(solver->model);
      if (newCost < ubCost || nbSatisfiable == 1) {
        saveModel(solver->model);
        if (maxsat_formula->getFormat() == _FORMAT_PB_) {
          // optimization problem
          if (maxsat_formula->getObjFunction() != NULL) {
            printf("o %" PRId64 "\n", newCost + off_set);
          }
        } else
          printf("o %" PRId64 "\n", newCost + off_set);
        ubCost = newCost;
      }

      if (nbSatisfiable == 1) {
        min_weight =
            findNextWeightDiversity(min_weight, cardinality_assumptions);
        // printf("current weight %d\n",min_weight);

        for (int i = 0; i < maxsat_formula->nSoft(); i++)
          if (maxsat_formula->getSoftClause(i).weight >= min_weight)
            assumptions.push(~maxsat_formula->getSoftClause(i).assumption_var);
      } else {
        // compute min weight in soft
        int not_considered = 0;
        for (int i = 0; i < maxsat_formula->nSoft(); i++) {
          if (maxsat_formula->getSoftClause(i).weight < min_weight)
            not_considered++;
        }

        // printf("not considered %d\n",not_considered);

        for (std::set<Lit>::iterator it = cardinality_assumptions.begin();
             it != cardinality_assumptions.end(); ++it) {
          assert(boundMapping.find(*it) != boundMapping.end());
          std::pair<std::pair<int, uint64_t>, uint64_t> soft_id =
              boundMapping[*it];
          if (soft_id.second < min_weight)
            not_considered++;
        }

        if (not_considered != 0) {
          min_weight =
              findNextWeightDiversity(min_weight, cardinality_assumptions);

          // printf("currentWeight %d\n",currentWeight);

          // reset the assumptions
          assumptions.clear();
          int active_soft = 0;
          for (int i = 0; i < maxsat_formula->nSoft(); i++) {
            if (!activeSoft[i] &&
                maxsat_formula->getSoftClause(i).weight >= min_weight) {
              assumptions.push(
                  ~maxsat_formula->getSoftClause(i).assumption_var);
              // printf("s assumption
              // %d\n",var(softClauses[i].assumptionVar)+1);
            } else
              active_soft++;
          }

          // printf("assumptions %d\n",assumptions.size());
          for (std::set<Lit>::iterator it = cardinality_assumptions.begin();
               it != cardinality_assumptions.end(); ++it) {
            assert(boundMapping.find(*it) != boundMapping.end());
            std::pair<std::pair<int, uint64_t>, uint64_t> soft_id =
                boundMapping[*it];
            if (soft_id.second >= min_weight)
              assumptions.push(~(*it));
            // printf("c assumption %d\n",var(*it)+1);
          }

        } else {
          assert(lbCost == newCost);
          printAnswer(_OPTIMUM_);
          exit(_OPTIMUM_);
        }
      }
    }

    if (res == l_False) {

      // reduce the weighted to the unweighted case
      uint64_t min_core = UINT64_MAX;
      for (int i = 0; i < solver->conflict.size(); i++) {
        Lit p = solver->conflict[i];
        if (coreMapping.find(p) != coreMapping.end()) {
          assert(!activeSoft[coreMapping[p]]);
          if (maxsat_formula->getSoftClause(coreMapping[solver->conflict[i]])
                  .weight < min_core)
            min_core =
                maxsat_formula->getSoftClause(coreMapping[solver->conflict[i]])
                    .weight;
        }

        if (boundMapping.find(p) != boundMapping.end()) {
          std::pair<std::pair<int, uint64_t>, uint64_t> soft_id =
              boundMapping[solver->conflict[i]];
          if (soft_id.second < min_core)
            min_core = soft_id.second;
        }
      }

      // printf("MIN CORE %d\n",min_core);

      lbCost += min_core;
      nbCores++;
      if (verbosity > 0)
        printf("c LB : %-12" PRIu64 "\n", lbCost);

      if (nbSatisfiable == 0) {
        printAnswer(_UNSATISFIABLE_);
        exit(_UNSATISFIABLE_);
      }

      if (lbCost == ubCost) {
        assert(nbSatisfiable > 0);
        if (verbosity > 0)
          printf("c LB = UB\n");
        printAnswer(_OPTIMUM_);
        exit(_OPTIMUM_);
      }

      sumSizeCores += solver->conflict.size();

      vec<Lit> soft_relax;
      vec<Lit> cardinality_relax;

      for (int i = 0; i < solver->conflict.size(); i++) {
        Lit p = solver->conflict[i];
        if (coreMapping.find(p) != coreMapping.end()) {
          if (maxsat_formula->getSoftClause(coreMapping[p]).weight > min_core) {
            // printf("SPLIT THE CLAUSE\n");
            assert(!activeSoft[coreMapping[p]]);
            // SPLIT THE CLAUSE
            int indexSoft = coreMapping[p];
            assert(maxsat_formula->getSoftClause(indexSoft).weight - min_core >
                   0);

            // Update the weight of the soft clause.
            maxsat_formula->getSoftClause(indexSoft).weight -= min_core;

            vec<Lit> clause;
            vec<Lit> vars;

            maxsat_formula->getSoftClause(indexSoft).clause.copyTo(clause);
            // Since cardinality constraints are added the variables are not
            // in sync...
            while (maxsat_formula->nVars() < solver->nVars())
              maxsat_formula->newLiteral();

            Lit l = maxsat_formula->newLiteral();
            vars.push(l);

            // Add a new soft clause with the weight of the core.
            // addSoftClause(softClauses[indexSoft].weight-min_core, clause,
            // vars);
            maxsat_formula->addSoftClause(min_core, clause, vars);
            activeSoft.push(true);

            // Add information to the SAT solver
            newSATVariable(solver);
            clause.push(l);
            solver->addClause(clause);

            assert(clause.size() - 1 ==
                   maxsat_formula->getSoftClause(indexSoft).clause.size());
            assert(maxsat_formula->getSoftClause(maxsat_formula->nSoft() - 1)
                       .relaxation_vars.size() == 1);

            // Create a new assumption literal.

            maxsat_formula->getSoftClause(maxsat_formula->nSoft() - 1)
                .assumption_var = l;
            assert(maxsat_formula->getSoftClause(maxsat_formula->nSoft() - 1)
                       .assumption_var ==
                   maxsat_formula->getSoftClause(maxsat_formula->nSoft() - 1)
                       .relaxation_vars[0]);
            coreMapping[l] =
                maxsat_formula->nSoft() -
                1; // Map the new soft clause to its assumption literal.

            soft_relax.push(l);
            assert(maxsat_formula->getSoftClause(coreMapping[l]).weight ==
                   min_core);
            assert(activeSoft.size() == maxsat_formula->nSoft());

          } else {
            // printf("NOT SPLITTING\n");
            assert(
                maxsat_formula->getSoftClause(coreMapping[solver->conflict[i]])
                    .weight == min_core);
            soft_relax.push(p);
            // printf("ASSERT %d\n",var(p)+1);
            assert(!activeSoft[coreMapping[p]]);
            activeSoft[coreMapping[p]] = true;
          }
        }

        if (boundMapping.find(p) != boundMapping.end()) {
          // printf("CARD IN CORE\n");

          std::set<Lit>::iterator it;
          it = cardinality_assumptions.find(p);
          assert(it != cardinality_assumptions.end());

          // this is a soft cardinality -- bound must be increased
          std::pair<std::pair<int, uint64_t>, uint64_t> soft_id =
              boundMapping[solver->conflict[i]];
          // increase the bound
          assert(soft_id.first.first < soft_cardinality.size());
          assert(soft_cardinality[soft_id.first.first]->hasCardEncoding());

          if (soft_id.second == min_core) {

            cardinality_assumptions.erase(it);
            cardinality_relax.push(p);

            joinObjFunction.clear();
            encodingAssumptions.clear();
            soft_cardinality[soft_id.first.first]->incUpdateCardinality(
                solver, joinObjFunction,
                soft_cardinality[soft_id.first.first]->lits(),
                soft_id.first.second + 1, encodingAssumptions);

            // if the bound is the same as the number of lits then no
            // restriction is applied
            if (soft_id.first.second + 1 <
                (unsigned)soft_cardinality[soft_id.first.first]
                    ->outputs()
                    .size()) {
              assert((unsigned)soft_cardinality[soft_id.first.first]
                         ->outputs()
                         .size() > soft_id.first.second + 1);
              Lit out = soft_cardinality[soft_id.first.first]
                            ->outputs()[soft_id.first.second + 1];
              boundMapping[out] = std::make_pair(
                  std::make_pair(soft_id.first.first, soft_id.first.second + 1),
                  min_core);
              cardinality_assumptions.insert(out);
            }

          } else {

#if 1
            // duplicate cardinality constraint???
            Encoder *e = new Encoder();
            e->setIncremental(_INCREMENTAL_ITERATIVE_);
            e->buildCardinality(solver,
                                soft_cardinality[soft_id.first.first]->lits(),
                                soft_id.first.second);

            assert((unsigned)e->outputs().size() > soft_id.first.second);
            Lit out = e->outputs()[soft_id.first.second];
            soft_cardinality.push(e);

            boundMapping[out] =
                std::make_pair(std::make_pair(soft_cardinality.size() - 1,
                                              soft_id.first.second),
                               min_core);
            cardinality_relax.push(out);

            // Update value of the previous cardinality constraint
            assert(soft_id.second - min_core > 0);
            boundMapping[p] = std::make_pair(
                std::make_pair(soft_id.first.first, soft_id.first.second),
                soft_id.second - min_core);

            // Update bound as usual...

            std::pair<std::pair<int, int>, int> soft_core_id =
                boundMapping[out];

            joinObjFunction.clear();
            encodingAssumptions.clear();
            soft_cardinality[soft_core_id.first.first]->incUpdateCardinality(
                solver, joinObjFunction,
                soft_cardinality[soft_core_id.first.first]->lits(),
                soft_core_id.first.second + 1, encodingAssumptions);

            // if the bound is the same as the number of lits then no
            // restriction is applied
            if (soft_core_id.first.second + 1 <
                soft_cardinality[soft_core_id.first.first]->outputs().size()) {
              assert(
                  soft_cardinality[soft_core_id.first.first]->outputs().size() >
                  soft_core_id.first.second + 1);
              Lit out = soft_cardinality[soft_core_id.first.first]
                            ->outputs()[soft_core_id.first.second + 1];
              boundMapping[out] =
                  std::make_pair(std::make_pair(soft_core_id.first.first,
                                                soft_core_id.first.second + 1),
                                 min_core);
              cardinality_assumptions.insert(out);
            }
#else

            // FIXME: improve the OLL algorithm by reusing the cardinality
            // constraint This current alternative does not work!
            while (nVars() < solver->nVars())
              newLiteral();

            Lit l = newLiteral();

            // Add information to the SAT solver
            newSATVariable(solver);
            vec<Lit> clause;
            clause.push(l);
            clause.push(~p);
            solver->addClause(clause);

            clause.clear();
            clause.push(~l);
            clause.push(p);
            solver->addClause(clause);

            boundMapping[l] = std::make_pair(
                std::make_pair(soft_id.first.first, soft_id.first.second),
                min_core);
            cardinality_relax.push(l);

            // Update bound as usual...

            std::pair<std::pair<int, int>, int> soft_core_id = boundMapping[p];

            joinObjFunction.clear();
            encodingAssumptions.clear();
            soft_cardinality[soft_core_id.first.first]->incUpdateCardinality(
                solver, joinObjFunction,
                soft_cardinality[soft_core_id.first.first]->lits(),
                soft_core_id.first.second + 1, encodingAssumptions);

            // if the bound is the same as the number of lits then no
            // restriction is applied
            if (soft_core_id.first.second + 1 <
                soft_cardinality[soft_core_id.first.first]->outputs().size()) {
              assert(
                  soft_cardinality[soft_core_id.first.first]->outputs().size() >
                  soft_core_id.first.second + 1);
              Lit out = soft_cardinality[soft_core_id.first.first]
                            ->outputs()[soft_core_id.first.second + 1];
              boundMapping[out] =
                  std::make_pair(std::make_pair(soft_core_id.first.first,
                                                soft_core_id.first.second + 1),
                                 min_core);
              cardinality_assumptions.insert(out);
            }

            // Update value of the previous cardinality constraint
            assert(soft_id.second - min_core > 0);
            boundMapping[p] = std::make_pair(
                std::make_pair(soft_id.first.first, soft_id.first.second),
                soft_id.second - min_core);

#endif
          }
        }
      }

      assert(soft_relax.size() + cardinality_relax.size() > 0);

      if (soft_relax.size() == 1 && cardinality_relax.size() == 0) {
        // Unit core
        // printf("UNIT CORE\n");
        solver->addClause(soft_relax[0]);
      }

      // assert (soft_relax.size() > 0 || cardinality_relax.size() != 1);

      if (soft_relax.size() + cardinality_relax.size() > 1) {

        vec<Lit> relax_harden;
        soft_relax.copyTo(relax_harden);
        for (int i = 0; i < cardinality_relax.size(); i++)
          relax_harden.push(cardinality_relax[i]);

        /*
        for (int i = 0; i < nSoft(); i++){
          printf("soft %d, w: %d, r:
        %d\n",i,softClauses[i].weight,var(softClauses[i].relaxationVars[0])+1);
          for (int j = 0; j < softClauses[i].clause.size(); j++){
            printf("%d ",var(softClauses[i].clause[j])+1);
          }
          printf("\n");
        }

        for (int i = 0; i < relax_harden.size(); i++)
          printf("+1 x%d ",var(relax_harden[i])+1);
        printf(" <= 1\n");
        */
        Encoder *e = new Encoder();
        e->setIncremental(_INCREMENTAL_ITERATIVE_);
        e->buildCardinality(solver, relax_harden, 1);
        soft_cardinality.push(e);
        assert(e->outputs().size() > 1);

        // printf("outputs %d\n",e->outputs().size());
        Lit out = e->outputs()[1];
        boundMapping[out] = std::make_pair(
            std::make_pair(soft_cardinality.size() - 1, 1), min_core);
        cardinality_assumptions.insert(out);
      }

      // reset the assumptions
      assumptions.clear();
      int active_soft = 0;
      for (int i = 0; i < maxsat_formula->nSoft(); i++) {
        if (!activeSoft[i] &&
            maxsat_formula->getSoftClause(i).weight >= min_weight) {
          assumptions.push(~maxsat_formula->getSoftClause(i).assumption_var);
          // printf("s assumption %d\n",var(softClauses[i].assumptionVar)+1);
        } else
          active_soft++;
      }

      // printf("assumptions %d\n",assumptions.size());
      for (std::set<Lit>::iterator it = cardinality_assumptions.begin();
           it != cardinality_assumptions.end(); ++it) {
        assert(boundMapping.find(*it) != boundMapping.end());
        std::pair<std::pair<int, uint64_t>, uint64_t> soft_id =
            boundMapping[*it];
        if (soft_id.second >= min_weight)
          assumptions.push(~(*it));
        // printf("c assumption %d\n",var(*it)+1);
      }

      // printf("card assumptions %d\n",assumptions.size());

      if (verbosity > 0) {
        printf("c Relaxed soft clauses %d / %d\n", active_soft,
               maxsat_formula->nSoft());
      }
    }
  }
}

void OLL::search() {

  if (encoding != _CARD_TOTALIZER_) {
    printf("Error: Currently algorithm MSU3 with iterative encoding only "
           "supports the totalizer encoding.\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }

  printConfiguration();

  if (maxsat_formula->getProblemType() == _WEIGHTED_) {
    // FIXME: consider lexicographical optimization for weighted problems
    weighted();
  } else
    unweighted();
}

/************************************************************************************************
 //
 // Rebuild MaxSAT solver
 //
 ************************************************************************************************/

/*_________________________________________________________________________________________________
  |
  |  rebuildSolver : [void]  ->  [Solver *]
  |
  |  Description:
  |
  |    Rebuilds a SAT solver with the current MaxSAT formula.
  |
  |________________________________________________________________________________________________@*/
Solver *OLL::rebuildSolver() {

  Solver *S = newSATSolver();

  for (int i = 0; i < maxsat_formula->nVars(); i++)
    newSATVariable(S);

  for (int i = 0; i < maxsat_formula->nHard(); i++)
    S->addClause(maxsat_formula->getHardClause(i).clause);

  // printf("c #PB: %d\n", maxsat_formula->nPB());
  for (int i = 0; i < maxsat_formula->nPB(); i++) {
    Encoder *enc = new Encoder(_INCREMENTAL_NONE_, _CARD_MTOTALIZER_,
                               _AMO_LADDER_, _PB_GTE_);

    // Make sure the PB is on the form <=
    if (!maxsat_formula->getPBConstraint(i)->_sign)
      maxsat_formula->getPBConstraint(i)->changeSign();

    enc->encodePB(S, maxsat_formula->getPBConstraint(i)->_lits,
                  maxsat_formula->getPBConstraint(i)->_coeffs,
                  maxsat_formula->getPBConstraint(i)->_rhs);

    // maxsat_formula->getPBConstraint(i)->print();

    delete enc;
  }

  // printf("c #Card: %d\n", maxsat_formula->nCard());
  for (int i = 0; i < maxsat_formula->nCard(); i++) {
    Encoder *enc = new Encoder(_INCREMENTAL_NONE_, _CARD_MTOTALIZER_,
                               _AMO_LADDER_, _PB_GTE_);

    if (maxsat_formula->getCardinalityConstraint(i)->_rhs == 1) {
      enc->encodeAMO(S, maxsat_formula->getCardinalityConstraint(i)->_lits);
    } else {
      enc->encodeCardinality(S,
                             maxsat_formula->getCardinalityConstraint(i)->_lits,
                             maxsat_formula->getCardinalityConstraint(i)->_rhs);
    }

    delete enc;
  }

  vec<Lit> clause;
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    clause.clear();
    maxsat_formula->getSoftClause(i).clause.copyTo(clause);
    for (int j = 0; j < maxsat_formula->getSoftClause(i).relaxation_vars.size();
         j++)
      clause.push(maxsat_formula->getSoftClause(i).relaxation_vars[j]);

    S->addClause(clause);
  }

  return S;
}

/************************************************************************************************
 //
 // Other protected methods
 //
 ************************************************************************************************/

/*_________________________________________________________________________________________________
  |
  |  initRelaxation : [void] ->  [void]
  |
  |  Description:
  |
  |    Initializes the relaxation variables by adding a fresh variable to the
  |    'relaxationVars' of each soft clause.
  |
  |  Post-conditions:
  |    * 'objFunction' contains all relaxation variables that were added to soft
  |      clauses.
  |
  |________________________________________________________________________________________________@*/
void OLL::initRelaxation() {
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    Lit l = maxsat_formula->newLiteral();
    maxsat_formula->getSoftClause(i).relaxation_vars.push(l);
    maxsat_formula->getSoftClause(i).assumption_var = l;
  }
}
