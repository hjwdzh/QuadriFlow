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

using namespace openwbo;

Graph::Graph(int nVert) {
  _nSCC = 0;
  _nVert = nVert;
  _edges.growTo(_nVert);
  _weights.growTo(_nVert);
  _incomingEdges.growTo(_nVert);
  _totalWeights.growTo(_nVert);
  _nSelfLoops.growTo(_nVert);

  _nMarked = 0;
  _totalWeight = 0.0;

  for (int i = 0; i < _nVert; i++) {
    _incomingEdges[i] = 0;
    _totalWeights[i] = 0.0;
    _nSelfLoops[i] = 0.0;
  }
}

Graph::~Graph() {}

void Graph::addEdge(int u, int v, double w) {
  if (u == v)
    _nSelfLoops[u] += w;

  _edges[u].push(v);
  _weights[u].push(w);
  _incomingEdges[v]++;

  _totalWeights[u] += w;
  _totalWeight += w;
}

int Graph::nEdges() {
  int e = 0;
  mergeDuplicatedEdges();
  for (int i = 0; i < _nVert; i++)
    e += _edges[i].size();
  return e;
}

// Remove duplicated edges. O(V+E)

void Graph::mergeDuplicatedEdges() {
  // vmm - Vectors on stack?!?! Is this a good idea for large graphs?!?!
  vec<int> destVertexes;
  vec<double> w;
  int n = 0;

  w.growTo(_nVert);
  destVertexes.growTo(_nVert);
  for (int i = 0; i < _nVert; i++)
    w[i] = 0;

  for (int i = 0; i < _nVert; i++) {
    if (_edges[i].size() > 0) {
      for (int j = 0; j < _edges[i].size(); j++) {
        int v = _edges[i][j];
        if (w[v] == 0)
          destVertexes[n++] = v; // New destination vertex
        else
          _incomingEdges[v]--; // Duplicated edge
        w[v] += _weights[i][j];
      }

      for (int j = 0; j < n; j++) {
        _edges[i][j] = destVertexes[j];
        _weights[i][j] = w[destVertexes[j]];
        w[destVertexes[j]] = 0;
      }

      _edges[i].shrink_(_edges[i].size() - n);
      _weights[i].shrink_(_weights[i].size() - n);
      n = 0;
    }
  }
}
