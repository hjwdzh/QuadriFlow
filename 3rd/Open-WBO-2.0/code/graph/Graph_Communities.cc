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

#include <iostream>
#include <stdlib.h>

#include "Graph.h"
#include "Graph_Communities.h"

#include "mtl/Vec.h"

#include <map>

using namespace openwbo;

#define PRECISION 0.000001

Graph_Communities::Graph_Communities() {
  _nCommunities = 0;
  _modularity = 0.0;
  _g = NULL;
}

Graph_Communities::~Graph_Communities() {}

int Graph_Communities::findCommunities(int mode, Graph *g) {
  // mode indicates the method used to identify communities...
  // Currentely, just the unfolding method is implemented.

  // Clear data from previous run
  _g = g;
  if (_vertexCommunity.size() > 0)
    _vertexCommunity.clear();
  _vertexCommunity.growTo(_g->nVertexes());
  for (int i = 0; i < g->nVertexes(); i++)
    _vertexCommunity[i] = i;

  resetInternalData();

  bool improvement = true;
  int level = 0;
  Graph *g_old = NULL;

  do {
    improvement = iterate();
    _modularity = modularity();

    ++level;

    g_old = _g;                // Save ptr to current graph
    _g = nextIterationGraph(); // Generate next iteration graph

    if (level > 1)
      delete g_old; // Delete previous graph, but never delete the original
                    // graph!!
    g_old = NULL;

    // Update community of original vertexes
    for (int i = 0; i < g->nVertexes(); i++)
      _vertexCommunity[i] = _renumber[_vertexToComm[_vertexCommunity[i]]];

    resetInternalData();

    if (level == 1) // do at least one more computation if partition is provided
      improvement = true;
  } while (improvement);

  // if (g_old != NULL) delete g_old;  // This can never happen!!!!

  return _nCommunities;
}

/// Internal

bool Graph_Communities::iterate() {
  double new_mod = modularity();
  double cur_mod = new_mod;
  bool better = false;

  // Generates a random order of vertexes
  int random_order[_g->nVertexes()];
  for (int i = 0; i < _g->nVertexes(); i++)
    random_order[i] = i;

  for (int i = 0; i < _g->nVertexes() - 1; i++) {
    int rand_pos = rand() % (_g->nVertexes() - i) + i;
    int tmp = random_order[i];
    random_order[i] = random_order[rand_pos];
    random_order[rand_pos] = tmp;
  }

  // Cycle to improve modularity
  do {
    cur_mod = new_mod;

    // For each vertex, tries to mode it to an adjacent community such that
    // modularity is increased
    for (int i = 0; i < _g->nVertexes(); i++) {
      int vertex = random_order[i];
      int comm = _vertexToComm[vertex];
      double factor = _g->weightedDegree(vertex) / _g->totalWeight();

      // computation of all neighboring communities with edges to vertex
      computeAdjCommunities(vertex);

      // remove vertex from its community
      remove(vertex, comm, _adjWeight[comm]);

      // determine the best adjacent community to mode the vertex to
      int best_comm = comm;
      double best_variation = 0.0;
      for (int i = 0; i < _adjComm.size(); i++) {
        // Calculate the modularity variation
        double variation =
            _adjWeight[_adjComm[i]] - (_total[_adjComm[i]] * factor);
        if (variation > best_variation) {
          best_comm = _adjComm[i];
          best_variation = variation;
        }
      }

      // insert vertex in the best adjacent community
      insert(vertex, best_comm, _adjWeight[best_comm]);

      if (best_comm != comm)
        better = true;
    }

    new_mod = modularity();

  } while (new_mod - cur_mod > PRECISION);

  return better;
}

void Graph_Communities::computeAdjCommunities(int vertex) {
  vec<int> &edges = _g->vertexEdges(vertex);
  vec<double> &weights = _g->vertexWeights(vertex);

  // Reset internal vectors
  for (int i = 0; i < _adjComm.size(); i++) {
    _adjWeight[_adjComm[i]] = 0.0;
    _adjMarked[_adjComm[i]] = false;
  }
  _adjComm.clear();

  // Consider current community
  _adjComm.push(_vertexToComm[vertex]);
  _adjMarked[_vertexToComm[vertex]] = true;

  // Mark adjacent communities and calculate weights
  for (int i = 0; i < edges.size(); i++) {
    int u = edges[i];
    int comm = _vertexToComm[u];

    if (u != vertex) {
      if (!_adjMarked[comm]) {
        _adjMarked[comm] = true;
        _adjComm.push(comm);
      }
      _adjWeight[comm] += weights[i];
    }
  }
}

Graph *Graph_Communities::nextIterationGraph() {
  // Compute the new number of communities
  for (int i = 0; i < _g->nVertexes(); i++)
    _renumber[i] = 0;
  for (int i = 0; i < _g->nVertexes(); i++)
    _renumber[_vertexToComm[i]]++;

  _nCommunities = 0;
  for (int i = 0; i < _g->nVertexes(); i++)
    if (_renumber[i] != 0)
      _renumber[i] = _nCommunities++;

  // Compute which vertexes belong to each community
  _communities.clear();
  _communities.growTo(_nCommunities);
  for (int i = 0; i < _g->nVertexes(); i++)
    _communities[_renumber[_vertexToComm[i]]].push(i);

  // Compute new weighted graph with colapsed communities
  Graph *g2 = new Graph(_nCommunities);

  for (int comm = 0; comm < _nCommunities; comm++) {
    map<int, double> m;
    map<int, double>::iterator it;
    int comm_size = _communities[comm].size();

    for (int u = 0; u < comm_size; u++) {
      vec<int> &edges = _g->vertexEdges(_communities[comm][u]);
      vec<double> &weights = _g->vertexWeights(_communities[comm][u]);

      for (int i = 0; i < edges.size(); i++) {
        int v = edges[i];
        int new_id = _renumber[_vertexToComm[v]];

        it = m.find(new_id);
        if (it == m.end())
          m.insert(make_pair(new_id, weights[i]));
        else
          it->second += weights[i];
      }
    }

    for (it = m.begin(); it != m.end(); it++)
      g2->addEdge(comm, it->first, it->second);
  }

  g2->mergeDuplicatedEdges();

  return g2;
}

void Graph_Communities::resetInternalData() {
  _vertexToComm.clear();
  _inside.clear();
  _total.clear();
  _renumber.clear();
  _adjComm.clear();
  _adjWeight.clear();
  _adjMarked.clear();

  _vertexToComm.growTo(_g->nVertexes());
  _inside.growTo(_g->nVertexes());
  _total.growTo(_g->nVertexes());
  _renumber.growTo(_g->nVertexes());
  _adjWeight.growTo(_g->nVertexes());
  _adjMarked.growTo(_g->nVertexes());

  for (int i = 0; i < _g->nVertexes(); i++) {
    _vertexToComm[i] = i;
    _inside[i] = _g->nSelfLoops(i);
    _total[i] = _g->weightedDegree(i);
    _renumber[i] = -1;
    _adjWeight[i] = 0;
    _adjMarked[i] = false;
  }
}

double Graph_Communities::modularity() {
  double mod = 0.;
  double tw2 = _g->totalWeight() * _g->totalWeight();

  for (int i = 0; i < _g->nVertexes(); i++) {
    if (_total[i] > 0) {
      mod += _inside[i] / _g->totalWeight() - (_total[i] * _total[i]) / tw2;
    }
  }
  return mod;
}

void Graph_Communities::remove(int node, int comm, double dnodecomm) {
  assert(node >= 0 && node < _g->nVertexes());

  _total[comm] -= _g->weightedDegree(node);
  _inside[comm] -= 2 * dnodecomm + _g->nSelfLoops(node);
  _vertexToComm[node] = -1;
}

void Graph_Communities::insert(int node, int comm, double dnodecomm) {
  assert(node >= 0 && node < _g->nVertexes());

  _total[comm] += _g->weightedDegree(node);
  _inside[comm] += 2 * dnodecomm + _g->nSelfLoops(node);
  _vertexToComm[node] = comm;
}
