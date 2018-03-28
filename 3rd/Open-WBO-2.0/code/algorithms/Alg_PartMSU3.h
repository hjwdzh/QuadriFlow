/*!
 * \author Miguel Neves - miguel.neves@ist.utl.pt
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

#ifndef PARTMSU3_H_
#define PARTMSU3_H_

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "../Encoder.h"
#include "../MaxSAT_Partition.h"
#include <algorithm>
#include <deque>
#include <map>
#include <set>

namespace openwbo {

class PartMSU3 : public MaxSAT_Partition {

protected:
  struct TreeNode {
  private:
    vec<int> parts;
    TreeNode *parent;
    int64_t lb;
    Encoder *encoder;
    vec<Lit> *encoding_assumptions;

  public:
    inline TreeNode(TreeNode *parent = NULL)
        : parent(parent), lb(0), encoder(NULL), encoding_assumptions(NULL) {}
    inline TreeNode(vec<int> &parts, TreeNode *parent = NULL)
        : parent(parent), lb(0), encoder(NULL), encoding_assumptions(NULL) {
      parts.copyTo(this->parts);
    }

    inline void addPartition(int part) { parts.push(part); }
    inline void addPartitions(vec<int> &new_parts) {
      for (int i = 0; i < new_parts.size(); ++i) {
        parts.push(new_parts[i]);
      }
    }
    inline void setParent(TreeNode *node) { parent = node; }
    inline void setEncoder(Encoder *enc) { encoder = enc; }
    inline void setEncodingAssumptions(vec<Lit> *assumpts) {
      encoding_assumptions = assumpts;
    }

    inline vec<int> &getPartitions() { return parts; }
    inline TreeNode *getParent() { return parent; }
    inline int64_t getLowerBound() { return lb; }
    inline Encoder *getEncoder() { return encoder; }
    inline vec<Lit> *getEncodingAssumptions() { return encoding_assumptions; }

    inline bool hasParent() { return parent != NULL; }
    inline bool hasEncoder() { return encoder != NULL; }
    inline bool hasEncodingAssumptions() {
      return encoding_assumptions != NULL;
    }

    inline void incrementLowerBound(int64_t inc = 1) { lb += inc; }
  };

public:
  PartMSU3(int verb = _VERBOSITY_MINIMAL_, int merge = _PART_BINARY_,
           int graph = RES_GRAPH, int enc = _CARD_TOTALIZER_) {
    solver = NULL;
    verbosity = verb;
    merge_strategy = merge;
    graph_type = graph;
    incremental_strategy = _INCREMENTAL_ITERATIVE_;
    encoding = enc;
  }
  virtual ~PartMSU3() {
    if (this->solver != NULL) {
      delete this->solver;
    }
  }

  void search();

  // Print solver configuration.
  void printConfiguration() {

    printf("c ==========================================[ Solver Settings "
           "]============================================\n");
    printf("c |                                                                "
           "                                       |\n");

    print_PartMSU3_configuration();
    print_Card_configuration(encoding);
  }

  void createGraph() {
    if (nPartitions() == 0) {
      split(UNFOLDING_MODE, graph_type);
    }
  }

  int chooseAlgorithm();

protected:
  // Print PartMSU3 configuration.
  void print_PartMSU3_configuration();

  // Rebuild MaxSAT solver
  //
  Solver *rebuildSolver(); // Rebuild MaxSAT solver.

  void PartMSU3_sequential(); // MSU3 that merges partitions sequentially into a
                              // single partition
  void PartMSU3_binary(); // MSU3 that uses a binary tree to guide the partition
                          // merging process

  // Heuristics
  mpq_class *computeSparsity();
  int sparsityHeuristic();

  // Other
  void initRelaxation(); // Relaxes soft clauses.
  void computeGuideTree(deque<TreeNode *> &out_tree);
  void dumpGuideTree(vec<TreeNode *> &tree);
  void sortPartitions(vec<int> &out_parts);

  Solver *solver; // SAT Solver used as a black box.

  // Controls the type of graph that will be used in the partitioning algorithm
  int graph_type;
  // Controls the partition merging strategy used by the algorithm
  int merge_strategy;
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

#endif /* PARTMSU3_H_ */
