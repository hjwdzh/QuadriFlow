/*!
 * \author Vasco Manquinho - vmm@sat.inesc-id.pt
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

#include "core/SolverTypes.h"

#include "MaxSAT_Partition.h"
#include "graph/Graph.h"
#include "graph/Graph_Communities.h"

using namespace openwbo;

#define _EDGE_LIMIT_ 50000000

MaxSAT_Partition::MaxSAT_Partition() {
  _solver = NULL;

  _nRandomPartitions = 16;
  _nPartitions = 0;
  _randomSeed = 0;

  _graph = NULL;
}

MaxSAT_Partition::~MaxSAT_Partition() {
  if (_graph != NULL)
    delete _graph;
}

void MaxSAT_Partition::init() {
  static bool executed = false;

  if (!executed) {
    srand(_randomSeed);
    executed = true;
  }

  if (_graph != NULL)
    delete _graph;
  if (_solver != NULL)
    delete _solver;
  _solver = newSATSolver();

  for (int i = 0; i < maxsat_formula->nVars(); i++)
    newSATVariable(_solver);

  for (int i = 0; i < maxsat_formula->nHard(); i++)
    _solver->addClause(maxsat_formula->getHardClause(i).clause);

  _graphMappingVar.clear();
  _graphMappingHard.clear();
  _graphMappingSoft.clear();

  _graphMappingVar.growTo(maxsat_formula->nVars());
  _graphMappingHard.growTo(maxsat_formula->nHard());
  _graphMappingSoft.growTo(maxsat_formula->nSoft());

  _partitions.clear(true);

  _nPartitions = 0;
}

void MaxSAT_Partition::split(int mode, int graphType) {
  init();

  if (!_solver->okay()) {
    delete _solver;
    return;
  }

  if (mode == RAND_MODE)
    splitRandom();
  else {
    _graph = buildGraph(true, graphType);

    if (_graph == NULL) {
      // Graph was not built because of the edge limit...
      buildSinglePartition();
    } else {
      // printf("c Graph: #V: %d\t#E: %d\n", _graph->nVertexes(),
      // _graph->nEdges());

      _gc.findCommunities(mode, _graph);
      // printf("c %d Communities found\n", _gc.nCommunities());

      buildPartitions(graphType);
    }
  }

  delete _solver;
}

void MaxSAT_Partition::splitRandom() {
  _nPartitions = _nRandomPartitions;
  _partitions.growTo(_nPartitions);

  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    // Satisfied or unsatisfied clauses are not assigned a community
    if (!unassignedLiterals(maxsat_formula->getSoftClause(i).clause))
      _graphMappingSoft[i] = -1;
    else {
      int c = rand() % _nPartitions;
      _partitions[c].sclauses.push(i);
      _graphMappingSoft[i] = c;
    }
  }
}

int MaxSAT_Partition::unassignedLiterals(vec<Lit> &sc) {
  int u = 0;
  for (int i = 0; i < sc.size(); i++)
    if (_solver->value(sc[i]) == l_True)
      return 0;
    else if (_solver->value(sc[i]) == l_Undef)
      u++;
  return u;
}

bool MaxSAT_Partition::isUnsatisfied(vec<Lit> &sc) {
  for (int i = 0; i < sc.size(); i++)
    if (_solver->value(sc[i]) != l_False)
      return false;
  return true;
}

void MaxSAT_Partition::printClause(vec<Lit> &sc) {
  for (int i = 0; i < sc.size(); i++)
    printf("%d ", (sign(sc[i]) ? -(var(sc[i]) + 1) : (var(sc[i]) + 1)));
}

void MaxSAT_Partition::buildPartitions(int graphType) {
  _nPartitions = _gc.nCommunities();
  _partitions.growTo(_nPartitions);

  if (graphType == VIG_GRAPH)
    buildVIGPartitions();
  else if (graphType == CVIG_GRAPH)
    buildCVIGPartitions();
  else if (graphType == RES_GRAPH)
    buildRESPartitions();
}

void MaxSAT_Partition::buildSinglePartition() {
  _nPartitions = 1;
  _partitions.growTo(_nPartitions);

  for (int i = 0; i < maxsat_formula->nVars(); i++) {
    // Put all unassigned variables in single partition
    if (_solver->value(i) != l_Undef)
      _graphMappingVar[i] = -1;
    else {
      _graphMappingVar[i] = 0;
      _partitions[0].vars.push(i);
    }
  }

  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    // Put all unresolved clauses in single partition
    if (unassignedLiterals(maxsat_formula->getSoftClause(i).clause)) {
      _graphMappingSoft[i] = 0;
      _partitions[0].sclauses.push(i);
    } else
      _graphMappingSoft[i] = -1;
  }

  for (int ci = 0; ci < maxsat_formula->nHard(); ci++) {
    // Put all unresolved clauses in single partition
    if (unassignedLiterals(maxsat_formula->getHardClause(ci).clause)) {
      _graphMappingHard[ci] = 0;
      _partitions[0].hclauses.push(ci);
    } else
      _graphMappingHard[ci] = -1;
  }
}

void MaxSAT_Partition::buildVIGPartitions() {
  vec<int> w;
  w.growTo(_nPartitions);
  for (int i = 0; i < _nPartitions; i++)
    w[i] = 0;

  for (int i = 0; i < maxsat_formula->nVars(); i++) {
    if (_graphMappingVar[i] != -1) {
      _graphMappingVar[i] = _gc.vertexCommunity(_graphMappingVar[i]);
      _partitions[_graphMappingVar[i]].vars.push(i);
    }
  }

  // Use soft clauses to define incidence function of variables
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    // Compute which partition the soft clause belongs to...

    if (unassignedLiterals(maxsat_formula->getSoftClause(i).clause)) {
      int d = 0;
      for (int j = 0; j < maxsat_formula->getSoftClause(i).clause.size(); j++) {
        if (_solver->value(maxsat_formula->getSoftClause(i).clause[j]) ==
            l_Undef) {
          int v = var(maxsat_formula->getSoftClause(i).clause[j]);
          int p = _graphMappingVar[v];
          w[p]++;
          if (w[p] > d) {
            d = w[p];
            _graphMappingSoft[i] = p;
          }
        }
      }
      _partitions[_graphMappingSoft[i]].sclauses.push(i);
      // vmm - Reseting the counters could be done more effectively.
      for (int j = 0; j < _nPartitions; j++)
        w[j] = 0;
    }
  }

  for (int ci = 0; ci < maxsat_formula->nHard(); ci++) {
    // Compute which partition the hard clause belongs to...
    vec<Lit> &c = maxsat_formula->getHardClause(ci).clause;
    if (unassignedLiterals(c) == 0)
      continue;

    int d = 0;
    for (int i = 0; i < c.size(); i++) {
      if (_solver->value(c[i]) != l_Undef)
        continue;
      int v = var(c[i]);
      int p = _graphMappingVar[v];
      w[p]++;
      if (w[p] > d) {
        d = w[p];
        _graphMappingHard[ci] = p;
      }
    }
    _partitions[_graphMappingHard[ci]].hclauses.push(ci);
    // vmm - Reseting the counters could be done more effectively.
    for (int i = 0; i < _nPartitions; i++)
      w[i] = 0;
  }
}

void MaxSAT_Partition::buildCVIGPartitions() {
  for (int i = 0; i < maxsat_formula->nVars(); i++) {
    if (_graphMappingVar[i] != -1) {
      _graphMappingVar[i] = _gc.vertexCommunity(_graphMappingVar[i]);
      // if (_graphMappingVar[i] >= _nPartitions)
      // 	printf("c Invalid partition %d (%d) for Var %d\n",
      // _graphMappingVar[i], _nPartitions, i);
      _partitions[_graphMappingVar[i]].vars.push(i);
    }
  }

  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    if (_graphMappingSoft[i] != -1) {
      _graphMappingSoft[i] = _gc.vertexCommunity(_graphMappingSoft[i]);
      // if (_graphMappingSoft[i] >= _nPartitions)
      // 	printf("c Invalid partition %d (%d) for Soft clause %d\n",
      // _graphMappingSoft[i], _nPartitions, i);
      _partitions[_graphMappingSoft[i]].sclauses.push(i);
    }
  }

  for (int i = 0; i < maxsat_formula->nHard(); i++) {
    if (_graphMappingHard[i] != -1) {
      _graphMappingHard[i] = _gc.vertexCommunity(_graphMappingHard[i]);
      // if (_graphMappingHard[i] >= _nPartitions)
      // 	printf("c Invalid partition %d (%d) for Hard clause %d\n",
      // _graphMappingHard[i], _nPartitions, i);
      _partitions[_graphMappingHard[i]].hclauses.push(i);
    }
  }
}

void MaxSAT_Partition::buildRESPartitions() {
  for (int i = 0; i < maxsat_formula->nVars(); i++) {
    // Compute which partition the variable belongs to...
  }

  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    if (_graphMappingSoft[i] != -1) {
      _graphMappingSoft[i] = _gc.vertexCommunity(_graphMappingSoft[i]);
      _partitions[_graphMappingSoft[i]].sclauses.push(i);
    }
  }

  for (int i = 0; i < maxsat_formula->nHard(); i++) {
    if (_graphMappingHard[i] != -1) {
      _graphMappingHard[i] = _gc.vertexCommunity(_graphMappingHard[i]);
      _partitions[_graphMappingHard[i]].hclauses.push(i);
    }
  }
}

Graph *MaxSAT_Partition::buildGraph(bool weighted, int graphType) {
  if (graphType == VIG_GRAPH)
    return buildVIGGraph(weighted);
  else if (graphType == CVIG_GRAPH)
    return buildCVIGGraph(weighted);
  else if (graphType == RES_GRAPH)
    return buildRESGraph(weighted);
  else
    return NULL;
}

Graph *MaxSAT_Partition::buildVIGGraph(bool weighted) {
  int gVars = 0;
  double *graphWeight = new double[maxsat_formula->nVars()];

  for (int i = 0; i < maxsat_formula->nVars(); i++) {
    if (_solver->value(i) != l_Undef)
      _graphMappingVar[i] = -1;
    else {
      _graphMappingVar[i] = gVars++;
      graphWeight[i] = 1;
    }
  }

  // Use soft clauses to define incidence function of variables
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    // Only considers unresolved soft clauses
    int ul;
    if (ul = unassignedLiterals(maxsat_formula->getSoftClause(i).clause)) {
      for (int j = 0; j < maxsat_formula->getSoftClause(i).clause.size(); j++) {
        if (_solver->value(maxsat_formula->getSoftClause(i).clause[j]) ==
            l_Undef)
          graphWeight[var(maxsat_formula->getSoftClause(i).clause[j])] +=
              ((double)maxsat_formula->getSoftClause(i).weight) / ul;
      }
    }
    _graphMappingSoft[i] = -1;
  }

  // Initialize hard clauses...
  for (int i = 0; i < maxsat_formula->nHard(); i++)
    _graphMappingHard[i] = -1;

  Graph *g = new Graph(gVars);

  int nEdges = 0;
  for (int ci = 0; ci < maxsat_formula->nHard(); ci++) {
    vec<Lit> &c = maxsat_formula->getHardClause(ci).clause;
    int ul = unassignedLiterals(c); // returns 0 if c is satisfied
    if (ul == 0)
      continue;

    double w = (weighted ? (2.0 / (ul * (ul - 1))) : 1.0);
    for (int i = 0; i < c.size(); i++) {
      if (_solver->value(c[i]) != l_Undef)
        continue;

      for (int j = i + 1; j < c.size(); j++) {
        if (_solver->value(c[j]) != l_Undef)
          continue;

        int u = var(c[i]), v = var(c[j]);
        g->addEdge(_graphMappingVar[u], _graphMappingVar[v],
                   graphWeight[u] * graphWeight[v] * w);
        g->addEdge(_graphMappingVar[v], _graphMappingVar[u],
                   graphWeight[u] * graphWeight[v] * w);
        nEdges++;
      }
    }

    if (nEdges >= _EDGE_LIMIT_) {
      // cout << "c Graph is too large." << endl;
      delete[] graphWeight;
      delete g;
      return NULL;
    }
  }

  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    // Only adds soft clauses that are being considered in the working formula
    int ul = unassignedLiterals(maxsat_formula->getSoftClause(i).clause);
    if (ul == 0)
      continue;

    double w = (weighted ? (2.0 / (ul * (ul - 1))) : 1.0);
    for (int j = 0; j < maxsat_formula->getSoftClause(i).clause.size(); j++) {
      if (_solver->value(maxsat_formula->getSoftClause(i).clause[j]) != l_Undef)
        continue;

      for (int k = j + 1; k < maxsat_formula->getSoftClause(i).clause.size();
           k++) {
        if (_solver->value(maxsat_formula->getSoftClause(i).clause[k]) !=
            l_Undef)
          continue;

        int u = var(maxsat_formula->getSoftClause(i).clause[j]),
            v = var(maxsat_formula->getSoftClause(i).clause[k]);
        g->addEdge(_graphMappingVar[u], _graphMappingVar[v],
                   graphWeight[u] * graphWeight[v] * w);
        g->addEdge(_graphMappingVar[v], _graphMappingVar[u],
                   graphWeight[u] * graphWeight[v] * w);
        nEdges++;
      }
    }

    if (nEdges >= _EDGE_LIMIT_) {
      // cout << "c Graph is too large." << endl;
      delete[] graphWeight;
      delete g;
      return NULL;
    }
  }

  g->mergeDuplicatedEdges();
  delete[] graphWeight;
  return g;
}

Graph *MaxSAT_Partition::buildCVIGGraph(bool weighted) {
  int gVars = 0, sVars = 0, hVars = 0;
  double *graphWeight = new double[maxsat_formula->nVars()];

  for (int i = 0; i < maxsat_formula->nVars(); i++) {
    if (_solver->value(i) != l_Undef)
      _graphMappingVar[i] = -1;
    else {
      _graphMappingVar[i] = gVars++;
      graphWeight[i] = 1;
    }
  }

  // Use soft clauses to define incidence function of variables
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    // Only considers unresolved soft clauses
    int ul;
    if (ul = unassignedLiterals(maxsat_formula->getSoftClause(i).clause)) {
      _graphMappingSoft[i] = gVars + sVars++;
      for (int j = 0; j < maxsat_formula->getSoftClause(i).clause.size(); j++) {
        if (_solver->value(maxsat_formula->getSoftClause(i).clause[j]) ==
            l_Undef)
          graphWeight[var(maxsat_formula->getSoftClause(i).clause[j])] +=
              ((double)maxsat_formula->getSoftClause(i).weight) / ul;
      }
    } else {
      // if (isUnsatisfied(softClauses[i].clause))
      // 	cout << "c soft clause is unsatisfied..." << endl;
      _graphMappingSoft[i] = -1;
    }
  }

  // Initialize graphMappingHard
  for (int i = 0; i < maxsat_formula->nHard(); i++) {
    int ul = unassignedLiterals(maxsat_formula->getHardClause(i).clause);
    if (ul == 0)
      _graphMappingHard[i] = -1;
    else
      _graphMappingHard[i] = gVars + sVars + hVars++;
  }

  Graph *g = new Graph(gVars + sVars + hVars);
  int nEdges = 0;

  for (int ci = 0; ci < maxsat_formula->nHard(); ci++) {
    if (_graphMappingHard[ci] != -1) { // -1 if it is not unresolved
      vec<Lit> &c = maxsat_formula->getHardClause(ci).clause;
      int ul = unassignedLiterals(c);

      // printf("c Clause %d is unresolved\n", ci);

      // double w = (weighted ? (2.0 / (ul * (ul-1))) : 1.0);
      for (int i = 0; i < c.size(); i++) {
        if (_solver->value(c[i]) != l_Undef)
          continue;

        int u = var(c[i]);
        g->addEdge(_graphMappingVar[u], _graphMappingHard[ci],
                   ((double)graphWeight[u]) / ul);
        g->addEdge(_graphMappingHard[ci], _graphMappingVar[u],
                   ((double)graphWeight[u]) / ul);
        nEdges++;

        // printf("c Adding edge! #E: %d\n", g->nEdges());
      }

      if (nEdges >= _EDGE_LIMIT_) {
        printf("c Graph is too large.\n");
        delete[] graphWeight;
        delete g;
        return NULL;
      }
    }
  }

  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    // Only adds unresolved soft clauses
    if (_graphMappingSoft[i] != -1) { // -1 if it is not unresolved
      int ul = unassignedLiterals(maxsat_formula->getSoftClause(i).clause);
      // double w = (weighted ? (2.0 / (ul * (ul-1))) : 1.0);

      for (int j = 0; j < maxsat_formula->getSoftClause(i).clause.size(); j++) {
        if (_solver->value(maxsat_formula->getSoftClause(i).clause[j]) !=
            l_Undef)
          continue;

        int u = var(maxsat_formula->getSoftClause(i).clause[j]);
        g->addEdge(_graphMappingVar[u], _graphMappingSoft[i],
                   ((double)graphWeight[u]) / ul);
        g->addEdge(_graphMappingSoft[i], _graphMappingVar[u],
                   ((double)graphWeight[u]) / ul);
        nEdges++;
      }

      if (nEdges >= _EDGE_LIMIT_) {
        printf("c Graph is too large.\n");
        delete[] graphWeight;
        delete g;
        return NULL;
      }
    }
  }

  g->mergeDuplicatedEdges();
  delete[] graphWeight;
  return g;
}

int MaxSAT_Partition::markUnassignedLiterals(vec<Lit> &c, int *markedLits,
                                             bool v) {
  int u = 0;
  for (int i = 0; i < c.size(); i++) {
    if (_solver->value(c[i]) != l_Undef)
      continue;
    markedLits[toInt(c[i])] = v;
    u++;
  }
  return u;
}

Graph *MaxSAT_Partition::buildRESGraph(bool weighted) {
  int sVars = 0, hVars = 0;
  int nLits = maxsat_formula->nVars() * 2;
  double *graphWeight = new double[maxsat_formula->nVars()];
  vec<int> *litClauses = new vec<int>[nLits];
  int *markedLits = new int[nLits];

  for (int i = 0; i < nLits; ++i) {
    markedLits[i] = false;
  }
  for (int i = 0; i < maxsat_formula->nVars(); i++) {
    _graphMappingVar[i] = -1;
    graphWeight[i] = 1;
  }

  // Use soft clauses to define incidence function of variables
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    // Only considers unresolved soft clauses
    int ul;
    if (ul = unassignedLiterals(maxsat_formula->getSoftClause(i).clause)) {
      _graphMappingSoft[i] = sVars++;
      for (int j = 0; j < maxsat_formula->getSoftClause(i).clause.size(); j++) {
        if (_solver->value(maxsat_formula->getSoftClause(i).clause[j]) ==
            l_Undef)
          graphWeight[var(maxsat_formula->getSoftClause(i).clause[j])] +=
              ((double)maxsat_formula->getSoftClause(i).weight) / ul;
      }
    } else
      _graphMappingSoft[i] = -1;
  }

  // Initialize graphMappingHard
  for (int i = 0; i < maxsat_formula->nHard(); i++) {
    int ul = unassignedLiterals(maxsat_formula->getHardClause(i).clause);
    if (ul == 0)
      _graphMappingHard[i] = -1;
    else
      _graphMappingHard[i] = sVars + hVars++;
  }

  Graph *g = new Graph(sVars + hVars);
  int nEdges = 0;

  for (int ci = 0; ci < maxsat_formula->nHard(); ci++) {
    if (_graphMappingHard[ci] != -1) { // -1 if it is not unresolved
      vec<Lit> &c = maxsat_formula->getHardClause(ci).clause;

      for (int i = 0; i < c.size(); i++) {
        if (_solver->value(c[i]) != l_Undef)
          continue;
        litClauses[toInt(c[i])].push(ci);
      }
    }
  }

  for (int ci = 0; ci < maxsat_formula->nHard(); ci++) {
    if (_graphMappingHard[ci] != -1) { // -1 if it is not unresolved
      vec<Lit> &c = maxsat_formula->getHardClause(ci).clause;

      // Mark clause literals - returns number of unassigned literals
      int mrk = markUnassignedLiterals(c, markedLits, true);

      for (int i = 0; i < c.size(); i++) {
        if (_solver->value(c[i]) != l_Undef)
          continue;

        int li = toInt(~c[i]);
        for (int iter = 0; iter < litClauses[li].size(); iter++) {
          int ri = litClauses[li][iter];
          if (ri <= ci)
            continue; // avoid duplication checks

          vec<Lit> &rc = maxsat_formula->getHardClause(ri).clause;
          int rl = 0, ul = mrk - 1;

          for (int j = 0; j < rc.size(); j++) {
            if (_solver->value(rc[j]) != l_Undef)
              continue;

            // Counts number of resolution literals l and ~l
            if (markedLits[toInt(~rc[j])] == true)
              rl++;
            // Counts number of different literals in resulting resolution
            // clause
            else if (markedLits[toInt(rc[j])] != true)
              ul++;
          }

          if (rl == 0)
            printf("No way!! There must be at least one!!\n");
          if (rl == 1) {
            if (!weighted)
              ul = 1;
            g->addEdge(_graphMappingHard[ci], _graphMappingHard[ri], 1.0 / ul);
            g->addEdge(_graphMappingHard[ri], _graphMappingHard[ci], 1.0 / ul);
            nEdges++;
          }
        }

        // printf("%d Edges\n", nEdges);
        if (nEdges >= _EDGE_LIMIT_) {
          printf("c Graph is too large.\n");
          for (int i = 0; i < nLits; i++)
            litClauses[i].clear();
          delete[] litClauses;
          delete[] markedLits;
          delete[] graphWeight;
          delete g;
          return NULL;
        }
      }

      // Clear marked literals
      markUnassignedLiterals(c, markedLits, false);
    }
  }

  // Connect soft clauses with hard clauses!!!
  for (int ci = 0; ci < maxsat_formula->nSoft(); ci++) {
    if (_graphMappingSoft[ci] != -1) { // -1 if it is not unresolved
      vec<Lit> &c = maxsat_formula->getSoftClause(ci).clause;

      // Mark clause literals
      int mrk = markUnassignedLiterals(c, markedLits, true);

      for (int i = 0; i < c.size(); i++) {
        if (_solver->value(c[i]) != l_Undef)
          continue;

        int li = toInt(~c[i]);
        for (int iter = 0; iter < litClauses[li].size(); iter++) {
          int ri = litClauses[li][iter];
          // if (ri <= ci) continue; //avoid duplication checks

          vec<Lit> &rc = maxsat_formula->getHardClause(ri).clause;
          int rl = 0, ul = mrk - 1;

          for (int j = 0; j < rc.size(); j++) {
            if (_solver->value(rc[j]) != l_Undef)
              continue;

            // Counts number of resolution literals l and ~l
            if (markedLits[toInt(~rc[j])] == true)
              rl++;
            // Counts number of different literals in resulting resolution
            // clause
            else if (markedLits[toInt(rc[j])] != true)
              ul++;
          }

          if (rl == 0)
            printf("No way!! There must be at least one!!\n");
          if (rl == 1) {
            if (!weighted)
              ul = 1;
            g->addEdge(_graphMappingSoft[ci], _graphMappingHard[ri], 1.0 / ul);
            g->addEdge(_graphMappingHard[ri], _graphMappingSoft[ci], 1.0 / ul);
            nEdges++;
          }
        }
        if (nEdges >= _EDGE_LIMIT_) {
          printf("c Graph is too large.\n");
          for (int i = 0; i < nLits; i++)
            litClauses[i].clear();
          delete[] litClauses;
          delete[] markedLits;
          delete[] graphWeight;
          delete g;
          return NULL;
        }
      }

      // Clear marked literals
      markUnassignedLiterals(c, markedLits, false);
    }
  }

  // litClauses cleaning
  for (int i = 0; i < nLits; i++)
    litClauses[i].clear();
  delete[] litClauses;
  delete[] markedLits;

  g->mergeDuplicatedEdges();
  delete[] graphWeight;
  return g;
}
