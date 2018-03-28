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

#include "Alg_PartMSU3.h"

#include <gmpxx.h>
#include <iostream>

#include <algorithm>
#include <list>
#include <unordered_map>

#define NO_PAIR -1
#define ERASE -2

#define SPARSITY_HEURISTIC 0.04
#define CLAUSE_LIMIT 1000000
#define PARTITION_RATIO 0.8

using namespace openwbo;

template <typename T>
static inline void addVector(vec<T> &dest, vec<T> &source) {
  for (int i = 0; i < source.size(); ++i) {
    dest.push(source[i]);
  }
}

// Heuristic selection used in MaxSAT Evaluation 2015
mpq_class *PartMSU3::computeSparsity() {
  mpq_class *h_val_pointer = new mpq_class("0", 10);

  for (int i = 0; i < nPartitions(); ++i) {
    *h_val_pointer += adjacentPartitions(i).size();
  }
  *h_val_pointer /= nPartitions() * nPartitions();

  return h_val_pointer;
}

int PartMSU3::sparsityHeuristic() {
  int algorithm = _ALGORITHM_PART_MSU3_;
  mpq_class *h_val_pointer = computeSparsity();

  if (*h_val_pointer < SPARSITY_HEURISTIC) {
    algorithm = _ALGORITHM_MSU3_;
  }
  delete h_val_pointer;

  return algorithm;
}

// Decides if it uses the partition-based algorithm
int PartMSU3::chooseAlgorithm() {

  int algorithm = _ALGORITHM_PART_MSU3_;

  if (maxsat_formula->nSoft() + maxsat_formula->nHard() >= CLAUSE_LIMIT)
    algorithm = _ALGORITHM_MSU3_;
  else {
    createGraph();
    if (nPartitions() <= 1) {
      algorithm = _ALGORITHM_MSU3_;
    } else {
      if ((float)nPartitions() / maxsat_formula->nSoft() > PARTITION_RATIO)
        algorithm = _ALGORITHM_MSU3_;
      else
        algorithm = sparsityHeuristic();
    }
  }

  return algorithm;
}

/*
 * Computes the tree used to guided the partition merging process in the binary
 * partition algorithm. 'out_tree' is an output parameter containing the leaves
 * of the tree. Unit partitions are all merged in a single leaf. Partitions
 * with no soft clauses are excluded from the tree.
 * WARNING: assumes that the 'split' method has been invoked already.
 * WARNING: the output tree is not ordered
 */
void PartMSU3::computeGuideTree(std::deque<TreeNode *> &out_tree) {
  assert(nPartitions() > 0);
  assert(out_tree.empty());
  out_tree.clear();

  // If there is a single partition, return a single root node
  if (nPartitions() == 1) // FIXME: not doing this will cause a seg fault in the
                          // MaxSATPartition code
  {
    TreeNode *node = new TreeNode();
    node->addPartition(0);
    out_tree.push_back(node);
    return;
  }

  // Build initial adjacency matrix
  vec<std::unordered_map<int, double>> adj_matrix(nPartitions());
  std::list<int> active_parts;
  for (int i = 0; i < nPartitions(); ++i) {
    for (int j = 0; j < adjacentPartitions(i).size(); ++j) {
      if (i != adjacentPartitions(i)[j]) {
        adj_matrix[i][adjacentPartitions(i)[j]] =
            adjacentPartitionWeights(i)[j];
      }
    }
    active_parts.push_back(i);
  }

  // Discard partitions with 0 soft clauses and merge unit partitions
  std::list<int>::iterator active_it = active_parts.begin();
  vec<TreeNode *> tree(nPartitions(), NULL);
  vec<int> no_adj_parts; // unit partitions with no adjacencies
  while (active_it != active_parts.end()) {
    int part = *active_it;
    if (nPartitionSoft(part) <= 1) {
      // Select best active partition to merge with current partition
      int best = -1;
      double best_max = 0;
      std::unordered_map<int, double>::iterator it;
      active_it = active_parts.erase(active_it);
      for (it = adj_matrix[part].begin(); it != adj_matrix[part].end(); ++it) {
        if (best == -1 || best_max < adj_matrix[part][it->first]) {
          best = it->first;
          best_max = it->second;
        }
      }
      if (best == -1) // partition has no adjacencies
      {
        if (nPartitionSoft(part) == 1) {
          no_adj_parts.push(part);
        }
        if (tree[part] != NULL) {
          for (int i = 0; i < tree[part]->getPartitions().size(); ++i) {
            no_adj_parts.push(tree[part]->getPartitions()[i]);
          }
          delete tree[part];
          tree[part] = NULL;
        }
        continue;
      }

      // Update weights
      for (it = adj_matrix[part].begin(); it != adj_matrix[part].end(); ++it) {
        if (it->first != best) {
          std::unordered_map<int, double>::iterator best_entry =
              adj_matrix[best].find(it->first);
          if (best_entry == adj_matrix[best].end()) {
            adj_matrix[best][it->first] = it->second;
            adj_matrix[it->first][best] = it->second;
          } else {
            best_entry->second += it->second;
            adj_matrix[it->first][best] = best_entry->second;
          }
          adj_matrix[it->first].erase(part);
        }
      }
      adj_matrix[part].clear();
      adj_matrix[best].erase(part);

      // Update tree nodes
      if (tree[part] != NULL) {
        if (tree[best] == NULL) {
          tree[best] = tree[part];
        } else {
          tree[best]->addPartitions(tree[part]->getPartitions());
          delete tree[part];
        }
        tree[part] = NULL;
      }
      if (nPartitionSoft(part) == 1) {
        if (tree[best] == NULL) {
          TreeNode *new_leaf = new TreeNode();
          new_leaf->addPartition(part);
          tree[best] = new_leaf;
          // out_tree.push_back(new_leaf);
        } else {
          tree[best]->addPartition(part);
        }
      }
    } else {
      if (tree[part] == NULL) {
        TreeNode *new_leaf = new TreeNode();
        new_leaf->addPartition(part);
        tree[part] = new_leaf;
      } else {
        tree[part]->addPartition(part);
      }
      ++active_it;
    }
  }

  // If there were unit partitions with no adjacencies, merge them into a leaf
  if (no_adj_parts.size() > 0) {
    if (no_adj_parts.size() == 1 && active_parts.size() > 0) {
      tree[active_parts.front()]->addPartition(no_adj_parts[0]);
    } else {
      tree[no_adj_parts[0]] = new TreeNode(no_adj_parts);
      active_parts.push_back(no_adj_parts[0]);
    }
  }

  // Output the leaves
  for (int i = 0; i < tree.size(); ++i) {
    if (tree[i] != NULL) {
      out_tree.push_back(tree[i]);
    }
  }

  // Build guide tree by merging partitions
  vec<int> best_pairs(nPartitions());
  while (active_parts.size() > 1) {
    vec<int> delayed_parts_stack;
    for (std::list<int>::reverse_iterator it = active_parts.rbegin();
         it != active_parts.rend(); ++it) {
      best_pairs[*it] = NO_PAIR;
      delayed_parts_stack.push(*it);
    }
    while (delayed_parts_stack.size() > 0) {
      int part = delayed_parts_stack.last();
      delayed_parts_stack.pop();

      // Select best pair for partition 'part'
      double best_max = -1;
      int best_pair = NO_PAIR;
      for (std::unordered_map<int, double>::iterator it =
               adj_matrix[part].begin();
           it != adj_matrix[part].end(); ++it) {
        int it_best = best_pairs[it->first];
        if (it->second > best_max &&
            (it_best == NO_PAIR || it_best == part ||
             adj_matrix[it_best].find(it->first) == adj_matrix[it_best].end() ||
             adj_matrix[it_best][it->first] < it->second)) {
          best_max = it->second;
          best_pair = it->first;
        }
      }

      if (best_pair == NO_PAIR) // possible if number of partitions is odd or
                                // partition has no adjacencies
      {
        for (std::list<int>::iterator it = active_parts.begin();
             it != active_parts.end() && best_pairs[part] == NO_PAIR; ++it) {
          if (best_pairs[*it] == NO_PAIR && *it != part) {
            best_pairs[part] = *it;
            best_pairs[*it] = part;
          }
        }
      } else if (best_pairs[best_pair] != part) {
        // Check if the chosen 'best_pair' is conflicting with another
        // partition's choice
        if (best_pairs[best_pair] >= 0) {
          best_pairs[best_pairs[best_pair]] = NO_PAIR;
          delayed_parts_stack.push(best_pairs[best_pair]);
        }
        if (best_pairs[part] >= 0) {
          best_pairs[best_pairs[part]] = NO_PAIR;
          delayed_parts_stack.push(best_pairs[part]);
        }

        best_pairs[part] = best_pair;
        best_pairs[best_pair] = part;
      }
    }

    // Merge paired partitions
    for (std::list<int>::iterator active_it = active_parts.begin();
         active_it != active_parts.end();) {
      int part = *active_it;
      if (best_pairs[part] >= 0) {
        assert(best_pairs[best_pairs[part]] == part);
        int pair_part = best_pairs[part];
        best_pairs[pair_part] = ERASE;

        // Update partition weights
        for (std::unordered_map<int, double>::iterator it =
                 adj_matrix[pair_part].begin();
             it != adj_matrix[pair_part].end(); ++it) {
          if (it->first != part) {
            std::unordered_map<int, double>::iterator part_entry =
                adj_matrix[part].find(it->first);
            if (part_entry == adj_matrix[part].end()) {
              adj_matrix[part][it->first] = it->second;
              adj_matrix[it->first][part] = it->second;
            } else {
              part_entry->second += it->second;
              adj_matrix[it->first][part] = part_entry->second;
            }
            adj_matrix[it->first].erase(pair_part);
          }
        }
        adj_matrix[pair_part].clear();
        adj_matrix[part].erase(pair_part);

        // Create new parent node
        TreeNode *parent_node = new TreeNode();
        tree[part]->setParent(parent_node);
        tree[pair_part]->setParent(parent_node);
        parent_node->addPartitions(tree[part]->getPartitions());
        parent_node->addPartitions(tree[pair_part]->getPartitions());
        tree[part] = parent_node;

        ++active_it;
      } else if (best_pairs[part] == ERASE) {
        active_it = active_parts.erase(active_it);
      } else {
        ++active_it;
      }
    }
  }
}

void PartMSU3::dumpGuideTree(vec<TreeNode *> &tree) {
  vec<TreeNode *> tree_level;
  tree.copyTo(tree_level);

  printf("c Dumping guide tree:\n");
  int level = 0;
  while (tree_level.size() > 1) {
    printf("c\tLevel %d\n", level);
    vec<TreeNode *> tmp_tree;
    for (int i = 0; i < tree_level.size(); ++i) {
      for (int j = 0; j < i; ++j) {
        if (tree_level[i]->getParent() == tree_level[j]->getParent()) {
          printf("c\t\t");
          for (int k = 0; k < tree_level[j]->getPartitions().size(); ++k) {
            printf(" %d", tree_level[j]->getPartitions()[k]);
          }
          printf("\n");
          printf("c\t\t");
          for (int k = 0; k < tree_level[i]->getPartitions().size(); ++k) {
            printf(" %d", tree_level[i]->getPartitions()[k]);
          }
          printf("\n");
          tmp_tree.push(tree_level[i]->getParent());
        }
      }
    }

    assert(tmp_tree.size() > 0);
    for (int i = 0; i < tree_level.size(); ++i) {
      int j;
      for (j = 0;
           j < tmp_tree.size() && tree_level[i]->getParent() != tmp_tree[j];
           ++j)
        ;
      if (j == tmp_tree.size()) {
        tmp_tree.push(tree_level[i]);
      }
    }
    tmp_tree.copyTo(tree_level);

    ++level;
  }

  printf("c\tLevel %d\n", level);
  printf("c\t\t");
  for (int k = 0; k < tree_level[0]->getPartitions().size(); ++k) {
    printf(" %d", tree_level[0]->getPartitions()[k]);
  }
  printf("\n");
}

/*
 * Computes the order in which partitions should be added to the SAT solver.
 * 'out_parts' is an output parameter where the order is returned. Partitions
 * with 0 soft clauses are discarded.
 * WARNING: assumes that the 'split' method has been invoked already.
 */
void PartMSU3::sortPartitions(vec<int> &out_parts) {
  assert(nPartitions() > 0);
  out_parts.clear();

  // Initialize algorithm's structures
  int best_part, old_best_part;
  double best_sum, best_max;
  vec<int> unit_parts;
  vec<double> weights(nPartitions(), 0);
  best_part = -1, old_best_part = -2;
  best_max = best_sum = 0;

  // Sort partitions
  while (old_best_part != best_part) {
    old_best_part = best_part;
    if (best_max == 0) {
      // Select one of the remaining partitions with highest rank
      double best_rank = -1;
      for (int i = 0; i < nPartitions(); ++i) {
        if (weights[i] >= 0 &&
            nPartitionSoft(i) >
                0) // partitions with negative weight are already sorted
        {
          int nvars = (nPartitionVars(i) > 0) ? nPartitionVars(i) : 1;
          double rank = (nPartitionSoft(i) + nPartitionHard(i)) / nvars;
          if (rank > best_rank) {
            best_rank = rank;
            best_part = i;
          }
        }
      }

      // Compute maximum weight of chosen partition, if one was chosen
      if (best_rank >= 0) {
        for (int i = 0; i < adjacentPartitionWeights(best_part).size(); ++i) {
          if (weights[adjacentPartitions(best_part)[i]] >= 0) {
            double weight = adjacentPartitionWeights(best_part)[i];
            best_max = (weight > best_max) ? weight : best_max;
          }
        }
      }
    } else {
      // Update weights
      for (int i = 0; i < adjacentPartitions(best_part).size(); ++i) {
        int adj_part = adjacentPartitions(best_part)[i];
        if (weights[adj_part] >= 0) {
          weights[adj_part] += adjacentPartitionWeights(best_part)[i];
        }
      }

      // Select next partition
      for (int i = 0; i < weights.size(); ++i) {
        if (weights[i] > weights[best_part]) {
          best_part = i;
        }
      }
      best_max = weights[best_part];
    }

    // Add selected partition if it contains soft clauses
    if (nPartitionSoft(best_part) > 0 && old_best_part != best_part) {
      if (nPartitionSoft(best_part) == 1) {
        unit_parts.push(best_part); // separate unit partitions
      } else {
        out_parts.push(best_part);
      }
    }
    weights[best_part] = -1;
  }

  // All unit partitions are to be solved last
  for (int i = 0; i < unit_parts.size(); ++i) {
    out_parts.push(unit_parts[i]);
  }
}

void PartMSU3::PartMSU3_sequential() {
  // nbInitialVariables = nVars();
  lbool res = l_True;
  vec<Lit> assumptions;
  vec<Lit> joinObjFunction;
  vec<Lit> currentObjFunction;
  vec<Lit> encodingAssumptions;

  Encoder *encoder = new Encoder();
  encoder->setIncremental(incremental_strategy);

  // Initialize partitions
  int part_index = 0;
  vec<int> parts;
  bool add_unit_parts = false;
  vec<int> unit_parts;
  if (nPartitions() == 0) {
    split(UNFOLDING_MODE, graph_type);
  }
  printConfiguration();

  if (merge_strategy == _PART_SEQUENTIAL_SORTED_) {
    sortPartitions(parts);
  }

  // Build solver
  initRelaxation();
  solver = rebuildSolver();

  activeSoft.growTo(maxsat_formula->nSoft(), false);
  for (int i = 0; i < maxsat_formula->nSoft(); i++)
    coreMapping[getAssumptionLit(i)] = i;

  for (;;) {

    res = searchSATSolver(solver, assumptions);
    if (res == l_True) {
      nbSatisfiable++;
      uint64_t newCost = computeCostModel(solver->model);
      if (nbSatisfiable == 1 || newCost < ubCost) {
        saveModel(solver->model);
        printf("o %" PRIu64 "\n", newCost);
        ubCost = newCost;
      }

      if (merge_strategy == _PART_SEQUENTIAL_) {
        // Select next partition
        while (nPartitionSoft(part_index) < 2 && part_index < nPartitions()) {
          if (nPartitionSoft(part_index) == 1) {
            unit_parts.push(part_index);
          }
          part_index++;
        }

        // Check if there are partitions left to be added
        if (part_index < nPartitions()) {
          for (int i = 0; i < nPartitionSoft(part_index); i++) {
            assumptions.push(~objFunction[communitySoft(part_index)[i]]);
          }
          part_index++;
        } else if (!add_unit_parts) {
          for (int i = 0; i < unit_parts.size(); ++i) {
            assumptions.push(~objFunction[communitySoft(unit_parts[i])[0]]);
          }
          add_unit_parts =
              true; // to know that unit partitions are being considered
        } else {
          printAnswer(_OPTIMUM_);
          exit(_OPTIMUM_);
        }
      } else {
        assert(merge_strategy == _PART_SEQUENTIAL_SORTED_);

        // Check if there are partitions left to be added
        if (part_index < parts.size()) {
          if (nPartitionSoft(parts[part_index]) == 1) {
            while (part_index < parts.size()) {
              assumptions.push(
                  ~objFunction[communitySoft(parts[part_index++])[0]]);
            }
          } else {
            for (int i = 0; i < nPartitionSoft(parts[part_index]); i++) {
              assumptions.push(
                  ~objFunction[communitySoft(parts[part_index])[i]]);
            }
            part_index++;
          }
        } else {
          printAnswer(_OPTIMUM_);
          exit(_OPTIMUM_);
        }
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

      if (solver->conflict.size() == 0) {
        printAnswer(_UNSATISFIABLE_);
        exit(_UNSATISFIABLE_);
      }

      joinObjFunction.clear();
      for (int i = 0; i < solver->conflict.size(); i++) {
        if (coreMapping.find(solver->conflict[i]) != coreMapping.end()) {
          assert(!activeSoft[coreMapping[solver->conflict[i]]]);
          activeSoft[coreMapping[solver->conflict[i]]] = true;
          joinObjFunction.push(
              getRelaxationLit(coreMapping[solver->conflict[i]]));
        }
      }

      currentObjFunction.clear();
      assumptions.clear();
      if (merge_strategy == _PART_SEQUENTIAL_) {
        for (int i = 0; i < part_index; ++i) {
          if (nPartitionSoft(i) > 1) {
            for (int j = 0; j < nPartitionSoft(i); ++j) {
              int soft_index = communitySoft(i)[j];
              if (activeSoft[soft_index]) {
                currentObjFunction.push(getRelaxationLit(soft_index));
              } else {
                assumptions.push(~getAssumptionLit(soft_index));
              }
            }
          }
        }
        if (add_unit_parts) {
          for (int i = 0; i < unit_parts.size(); ++i) {
            int soft_index = communitySoft(unit_parts[i])[0];
            if (activeSoft[soft_index]) {
              currentObjFunction.push(getRelaxationLit(soft_index));
            } else {
              assumptions.push(~getAssumptionLit(soft_index));
            }
          }
        }
      } else {
        assert(merge_strategy == _PART_SEQUENTIAL_SORTED_);

        for (int i = 0; i < part_index; ++i) {
          for (int j = 0; j < nPartitionSoft(parts[i]); ++j) {
            int soft_index = communitySoft(parts[i])[j];
            if (activeSoft[soft_index]) {
              currentObjFunction.push(getRelaxationLit(soft_index));
            } else {
              assumptions.push(~getAssumptionLit(soft_index));
            }
          }
        }
      }

      if (verbosity > 0)
        printf("c Relaxed soft clauses %d / %d\n", currentObjFunction.size(),
               objFunction.size());

      if (!encoder->hasCardEncoding()) {
        if (lbCost != (unsigned)currentObjFunction.size()) {
          encoder->buildCardinality(solver, currentObjFunction, lbCost);
          encoder->incUpdateCardinality(solver, currentObjFunction, lbCost,
                                        encodingAssumptions);
        }
      } else {
        // // Incremental construction of the encoding.
        // if (joinObjFunction.size() > 0)
        //     encoder.joinEncoding(solver, joinObjFunction, lbCost);

        // The right-hand side is constrained using assumptions.
        // NOTE: 'encodingAsssumptions' is modified in 'incrementalUpdate'.
        encoder->incUpdateCardinality(solver, joinObjFunction,
                                      currentObjFunction, lbCost,
                                      encodingAssumptions);
      }

      addVector(assumptions, encodingAssumptions);
    }
  }
}

void PartMSU3::PartMSU3_binary() {

  int nrelaxed = 0;
  // nbInitialVariables = nVars();
  lbool res = l_True;
  vec<Lit> assumptions;
  vec<Lit> joinObjFunction;
  vec<Lit> currentObjFunction;

  // Initialize partitions, compute guide tree and create encoders
  TreeNode *current_node = NULL;
  std::deque<TreeNode *> guide_tree;

  if (nPartitions() == 0) {
    split(UNFOLDING_MODE, graph_type);
  }
  printConfiguration();

  // printf("c Computing guide tree\n");
  computeGuideTree(guide_tree);
  for (std::deque<TreeNode *>::iterator it = guide_tree.begin();
       it != guide_tree.end(); ++it) {
    (*it)->setEncoder(new Encoder(incremental_strategy, encoding));
    (*it)->setEncodingAssumptions(new vec<Lit>());
  }

  // Build solver
  initRelaxation();
  solver = rebuildSolver();
  // printf("solver vars %d\n",solver->nVars());

  activeSoft.growTo(maxsat_formula->nSoft(), false);
  for (int i = 0; i < maxsat_formula->nSoft(); i++)
    coreMapping[getAssumptionLit(i)] = i;

  for (;;) {
    res = searchSATSolver(solver, assumptions);
    if (res == l_True) {
      nbSatisfiable++;
      uint64_t newCost = computeCostModel(solver->model);
      if (nbSatisfiable == 1 || newCost < ubCost) {
        saveModel(solver->model);
        printf("o %" PRIu64 "\n", newCost);
        ubCost = newCost;
      }

      if (nbSatisfiable == 1) {
        // assert(part_index == 0);
        current_node = guide_tree.front();
        guide_tree.pop_front();
        for (int i = 0; i < current_node->getPartitions().size(); ++i) {
          int comm = current_node->getPartitions()[i];
          for (int j = 0; j < nPartitionSoft(comm); ++j) {
            assumptions.push(~getAssumptionLit(communitySoft(comm)[j]));
          }
        }
      } else if (current_node->hasParent()) // no parent -> current_node is root
      {
        TreeNode *parent = current_node->getParent();
        if (parent->hasEncoder()) {
          guide_tree.push_back(parent);

          // Merge partitions by joining their encoders and lower bounds
          // vec<Lit> currentObjFunction;
          currentObjFunction.clear();
          for (int i = 0; i < parent->getPartitions().size(); ++i) {
            int comm = parent->getPartitions()[i];
            for (int j = 0; j < nPartitionSoft(comm); ++j) {
              if (activeSoft[communitySoft(comm)[j]]) {
                currentObjFunction.push(objFunction[communitySoft(comm)[j]]);
              }
            }
          }
          int64_t new_lb =
              parent->getLowerBound() + current_node->getLowerBound();
          if (parent->getEncoder()->hasCardEncoding() &&
              current_node->getEncoder()->hasCardEncoding()) {
            parent->getEncoder()->addCardinality(
                solver, *current_node->getEncoder(), new_lb);
            parent->getEncodingAssumptions()->clear();
            parent->getEncoder()->incUpdateCardinality(
                solver, currentObjFunction, new_lb,
                *(parent->getEncodingAssumptions()));
            delete current_node->getEncoder();
            delete current_node->getEncodingAssumptions();
          } else if (parent->getEncoder()->hasCardEncoding() ||
                     current_node->getEncoder()->hasCardEncoding()) {
            int encoded_start_index, encoded_end_index;
            int join_start_index, join_end_index;
            if (current_node->getPartitions()[0] ==
                parent->getPartitions()[0]) {
              encoded_start_index = join_end_index =
                  current_node->getPartitions().size();
              encoded_end_index = parent->getPartitions().size();
              join_start_index = 0;
            } else {
              encoded_start_index = 0;
              encoded_end_index = join_start_index =
                  parent->getPartitions().size() -
                  current_node->getPartitions().size();
              join_end_index = parent->getPartitions().size();
            }

            // Choose encoder that will be maintained
            if (current_node->getEncoder()->hasCardEncoding()) {
              delete parent->getEncoder();
              delete parent->getEncodingAssumptions();
              parent->setEncoder(current_node->getEncoder());
              parent->setEncodingAssumptions(
                  current_node->getEncodingAssumptions());
              swap(encoded_start_index, join_start_index);
              swap(encoded_end_index, join_end_index);
            } else {
              delete current_node->getEncoder();
              delete current_node->getEncodingAssumptions();
            }

            // Retrieve relaxed soft clauses in source partition and join
            // vec<Lit> joinObjFunction;
            joinObjFunction.clear();
            for (int i = join_start_index; i < join_end_index; ++i) {
              int comm = parent->getPartitions()[i];
              for (int j = 0; j < nPartitionSoft(comm); ++j) {
                if (activeSoft[communitySoft(comm)[j]]) {
                  Lit relax_lit = objFunction[communitySoft(comm)[j]];
                  joinObjFunction.push(relax_lit);
                }
              }
            }
            if (joinObjFunction.size() > 0) {
              parent->getEncodingAssumptions()->clear();
              parent->getEncoder()->joinEncoding(solver, joinObjFunction,
                                                 new_lb);
              parent->getEncoder()->incUpdateCardinality(
                  solver,
                  // joinObjFunction,
                  currentObjFunction, new_lb,
                  *(parent->getEncodingAssumptions()));
            }
          } else {
            delete current_node->getEncoder();
            delete current_node->getEncodingAssumptions();
          }
        } else {
          parent->setEncoder(current_node->getEncoder());
          parent->setEncodingAssumptions(
              current_node->getEncodingAssumptions());
        }
        parent->incrementLowerBound(current_node->getLowerBound());
        delete current_node;

        // Replace current partition
        assumptions.clear();
        current_node = guide_tree.front();
        guide_tree.pop_front();
        for (int i = 0; i < current_node->getPartitions().size(); ++i) {
          int comm = current_node->getPartitions()[i];
          for (int j = 0; j < communitySoft(comm).size(); ++j) {
            if (!activeSoft[communitySoft(comm)[j]]) {
              assumptions.push(~objFunction[communitySoft(comm)[j]]);
            }
          }
        }
        addVector(assumptions, *(current_node->getEncodingAssumptions()));
      } else {
        assert(guide_tree.empty());
        printAnswer(_OPTIMUM_);
        exit(_OPTIMUM_);
      }
    }

    if (res == l_False) {
      current_node->incrementLowerBound();
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

      if (solver->conflict.size() == 0) {
        printAnswer(_UNSATISFIABLE_);
        exit(_UNSATISFIABLE_);
      }

      joinObjFunction.clear();
      for (int i = 0; i < solver->conflict.size(); i++) {
        if (coreMapping.find(solver->conflict[i]) != coreMapping.end()) {
          assert(!activeSoft[coreMapping[solver->conflict[i]]]);
          activeSoft[coreMapping[solver->conflict[i]]] = true;
          joinObjFunction.push(
              getRelaxationLit(coreMapping[solver->conflict[i]]));
          nrelaxed++;
        }
      }

      currentObjFunction.clear();
      assumptions.clear();
      for (int i = 0; i < current_node->getPartitions().size(); ++i) {
        int comm = current_node->getPartitions()[i];
        for (int j = 0; j < nPartitionSoft(comm); ++j) {
          int soft_index = communitySoft(comm)[j];
          if (activeSoft[soft_index]) {
            currentObjFunction.push(getRelaxationLit(soft_index));
          } else {
            assumptions.push(~getAssumptionLit(soft_index));
          }
        }
      }

      if (verbosity > 0)
        printf("c Relaxed soft clauses %d / %d\n", nrelaxed,
               objFunction.size());

      if (!current_node->getEncoder()->hasCardEncoding()) {
        if (current_node->getLowerBound() != currentObjFunction.size()) {
          current_node->getEncoder()->buildCardinality(
              solver, currentObjFunction, current_node->getLowerBound());
          current_node->getEncoder()->incUpdateCardinality(
              solver, currentObjFunction, current_node->getLowerBound(),
              *(current_node->getEncodingAssumptions()));
        }
      } else {
        // Incremental construction of the encoding.
        if (joinObjFunction.size() > 0) {
          current_node->getEncoder()->joinEncoding(
              solver, joinObjFunction, current_node->getLowerBound());
        }

        // The right-hand side is constrained using assumptions.
        // NOTE: 'encodingAsssumptions' is modified in 'incrementalUpdate'.
        current_node->getEncodingAssumptions()->clear();
        current_node->getEncoder()->incUpdateCardinality(
            solver,
            // joinObjFunction,
            currentObjFunction, current_node->getLowerBound(),
            *(current_node->getEncodingAssumptions()));
      }
      addVector(assumptions, *(current_node->getEncodingAssumptions()));
    }
  }
}

void PartMSU3::search() {
  if (maxsat_formula->getProblemType() == _WEIGHTED_) {
    printf("Error: Currently algorithm MSU3 does not support weighted MaxSAT "
           "instances.\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }

  if (incremental_strategy == _INCREMENTAL_ITERATIVE_) {
    if (encoding != _CARD_TOTALIZER_) {
      printf("Error: Currently iterative encoding in PartMSU3 only "
             "supports the Totalizer encoding.\n");
      printf("s UNKNOWN\n");
      exit(_ERROR_);
    }

    switch (merge_strategy) {
    case _PART_SEQUENTIAL_:
      PartMSU3_sequential();
      break;
    case _PART_SEQUENTIAL_SORTED_:
      PartMSU3_sequential();
      break;
    case _PART_BINARY_:
      PartMSU3_binary();
      break;
    default:
      printf("Error: No partition merging strategy.\n");
      printf("s UNKNOWN\n");
      exit(_ERROR_);
    }
  } else {
    printf("Error: No incremental strategy.\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }
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
Solver *PartMSU3::rebuildSolver() {
  Solver *S = newSATSolver();

  for (int i = 0; i < maxsat_formula->nVars(); i++)
    newSATVariable(S);

  for (int i = 0; i < maxsat_formula->nHard(); i++)
    S->addClause(getHardClause(i).clause);

  vec<Lit> clause;
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    clause.clear();
    getSoftClause(i).clause.copyTo(clause);
    for (int j = 0; j < getSoftClause(i).relaxation_vars.size(); j++)
      clause.push(getSoftClause(i).relaxation_vars[j]);

    S->addClause(clause);
  }

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
void PartMSU3::initRelaxation() {
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    Lit l = maxsat_formula->newLiteral();
    getSoftClause(i).relaxation_vars.push(l);
    getSoftClause(i).assumption_var = l;
    objFunction.push(l);
  }
}

void PartMSU3::print_PartMSU3_configuration() {
  printf("c |  Algorithm: %23s                                             "
         "                      |\n",
         "PartMSU3");
  switch (merge_strategy) {
  case _PART_SEQUENTIAL_:
    printf("c |  Partition Strategy: %14s                                    "
           "                         |\n",
           "Sequential");
    break;
  case _PART_SEQUENTIAL_SORTED_:
    printf("c |  Partition Strategy: %14s                                    "
           "                               |\n",
           "Seq-Sorted");
    break;
  case _PART_BINARY_:
    printf("c |  Partition Strategy: %14s                                    "
           "                               |\n",
           "Binary");
    break;
  }
  switch (graph_type) {
  case VIG_GRAPH:
    printf("c |  Graph Type: %22s                                            "
           "                       |\n",
           "VIG");
    break;
  case CVIG_GRAPH:
    printf("c |  Graph Type: %22s                                            "
           "                       |\n",
           "CVIG");
    break;
  case RES_GRAPH:
    printf("c |  Graph Type: %22s                                            "
           "                       |\n",
           "Resolution");
    break;
  }

  printf("c |  Number of partitions: %12d                                      "
         "                             |\n",
         nPartitions());
  printf("c |  Soft partition ratio: %12.5f                                    "
         "                               |\n",
         (float)nPartitions() / maxsat_formula->nSoft());
}
