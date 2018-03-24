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

#ifndef MAXSAT_PARTITION_H
#define MAXSAT_PARTITION_H

#include "MaxSAT.h"

#include "graph/Graph.h"
#include "graph/Graph_Communities.h"

#include <gmpxx.h>

using NSPACE::Var;

namespace openwbo {

enum graphType_ { VIG_GRAPH = 0, CVIG_GRAPH = 1, RES_GRAPH = 2 };

typedef struct {
  vec<int> vars;
  vec<int> sclauses;
  vec<int> hclauses;
} Partition;

class MaxSAT_Partition : public MaxSAT {

public:
  MaxSAT_Partition();
  ~MaxSAT_Partition();

  void split(int mode, int graphType = RES_GRAPH); // Default Value

  // Set number of Random Partitions
  void setRandomPartitions(int n) { _nRandomPartitions = n; }
  int getRandomPartitions() { return _nRandomPartitions; }

  // Set random seed
  void setRandomSeed(int n) { _randomSeed = n; }
  int getRandomSeed() { return _randomSeed; }

  double getModularity() { return _gc.getModularity(); }
  int nPartitions() { return _nPartitions; }
  int varPartition(Var v) { return _graphMappingVar[v]; }
  int hardClausePartition(int index) { return _graphMappingHard[index]; }
  int softClausePartition(int index) {
    if (index >= maxsat_formula->nSoft())
      return 0;
    else
      return _graphMappingSoft[index];
  }

  int nPartitionVars(int index) { return _partitions[index].vars.size(); }
  int nPartitionSoft(int index) { return _partitions[index].sclauses.size(); }
  int nPartitionHard(int index) { return _partitions[index].hclauses.size(); }

  const vec<int> &communityVars(int index) { return _partitions[index].vars; }
  const vec<int> &communitySoft(int index) {
    return _partitions[index].sclauses;
  }
  const vec<int> &communityHard(int index) {
    return _partitions[index].hclauses;
  }

  const vec<int> &adjacentPartitions(int index) {
    return _gc.adjCommunities(index);
  }
  const vec<double> &adjacentPartitionWeights(int index) {
    return _gc.adjCommunityWeights(index);
  }

  mpq_class *computeSparsity() {
    mpq_class *h_val_pointer = new mpq_class("0", 10);

    for (int i = 0; i < nPartitions(); ++i) {
      *h_val_pointer += adjacentPartitions(i).size();
    }
    *h_val_pointer /= nPartitions() * nPartitions();

    return h_val_pointer;
  }

  int nVertexes() { return _graph->nVertexes(); }
  int nEdges() { return _graph->nEdges(); }

protected:
  void init();

  void splitRandom();

  void buildPartitions(int graphType);
  void buildSinglePartition();
  void buildVIGPartitions();
  void buildCVIGPartitions();
  void buildRESPartitions();

  Graph *buildGraph(bool weighted, int graphType);
  Graph *buildVIGGraph(bool weighted);
  Graph *buildCVIGGraph(bool weighted);
  Graph *buildRESGraph(bool weighted);

  int unassignedLiterals(vec<Lit> &sc);
  bool isUnsatisfied(vec<Lit> &sc);

  int markUnassignedLiterals(vec<Lit> &c, int *markedLits, bool v);

  void printClause(vec<Lit> &sc);

protected:
  Solver *_solver;

  vec<int> _graphMappingVar;
  vec<int> _graphMappingHard;
  vec<int> _graphMappingSoft;

  int _randomSeed;
  int _nRandomPartitions;
  int _nPartitions;
  vec<Partition> _partitions;

  Graph *_graph;
  Graph_Communities _gc;
};

} // namespace openwbo

#endif // MAXSAT_PARTITION_H
