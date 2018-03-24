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

#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include "Graph.h"
#include "NodeHeap.h"

using namespace openwbo;

#define INF DBL_MAX

// Depth-first search and topological sort

void Graph::visitedVertexes(int u, vec<int> &reachedVertexes) {
  DFSVisit(u, reachedVertexes);
  for (int i = 0; i < reachedVertexes.size(); i++) {
    _marks[reachedVertexes[i]] = WHITE;
  }
}

void Graph::DFSVisitIter(int u, vec<int> &reachedVertexes) {
  vec<int> *l = new vec<int>;
  int n = 0;

  l->push(u);

  while (l->size() > 0) {
    n++;
    u = l->last();
    l->pop();
    if (_marks[u] == WHITE) {
      _marks[u] = BLACK;

      for (int i = 0; i < _edges[u].size(); i++) {
        if (_marks[_edges[u][i]] == WHITE) {
          l->push(_edges[u][i]);
        }
      }
      reachedVertexes.push(u);
    }
  }

  printf("Reached %d vertexes.\n", n);
  delete l;
}

void Graph::DFSVisit(int u, vec<int> &reachedVertexes) {
  if (_marks[u] == WHITE) {
    _marks[u] = BLACK;

    for (int i = 0; i < _edges[u].size(); i++) {
      if (_marks[_edges[u][i]] == WHITE) {
        DFSVisit(_edges[u][i], reachedVertexes);
      }
    }
    reachedVertexes.push(u);
  }
}

void Graph::topologicalSort(vec<int> &vertexes) {
  for (int i = 0; i < _nVert; i++) {
    if (_marks[i] == WHITE && _edges[i].size()) {
      DFSVisit(i, vertexes);
    }
  }
  for (int i = 0; i < vertexes.size(); i++)
    _marks[vertexes[i]] = WHITE;
}

int Graph::connectedComponents() {
  int n = 0;
  vec<int> vertexes;

  for (int i = 0; i < _nVert; i++) {
    if (_marks[i] == WHITE && _edges[i].size()) {
      n++;
      DFSVisitIter(i, vertexes);
    }
  }
  return n;
}

// Strong Connected Components (SCC)

// void Graph::findAllScc() {
//   _time = 1;
//   _stack.clear();

//   clearSCC();

//   for(int i = 0; i < _nVert; i++) {
//     if (_marks[i] == WHITE) {
//       tarjan(i);
//     }
//   }

//   for(int m = 0; m < _nMarked; m++) {
//     int i = _markedVertexes[m];
//     _marks[i] = WHITE;
//     _index[i] = _lowlink[i] = 0;
//   }
//   _nMarked = 0;

// }

// void Graph::tarjan(int u) {
//   _marks[u] = GRAY;
//   _markedVertexes[_nMarked++] = u;
//   _index[u] = _time;
//   _lowlink[u] = _time++;
//   _stack.push_front(u);

//   for (list<int>::iterator iter = _edges[u]->begin(); iter !=
//   _edges[u]->end(); iter++) {
//     int v = *iter;
//     if (_marks[v] == WHITE) {
//       tarjan(v);
//       if (_lowlink[v] < _lowlink[u]) _lowlink[u] = _lowlink[v];
//     }
//     else if (_marks[v] == GRAY) {
//       if (_lowlink[v] < _lowlink[u]) _lowlink[u] = _lowlink[v];
//     }
//   }

//   if (_index[u] == _lowlink[u]) {
//     int v = -1;
//     list<int>* scc = new list<int>;
//     do {
//       v = _stack.front();
//       _stack.pop_front();
//       scc->push_front(v);
//       _marks[v] = BLACK;
//     } while (v != u);

//     if (scc->size() == 1) {
//       scc->clear();
//       delete scc;
//     }
//     else {
//       _sccs[_nSCC++] = scc;
//     }
//   }
// }

// void Graph::clearSCC() {
//   for (int i = 0; i < _nSCC; i++) {
//     _sccs[i]->clear();
//     delete _sccs[i];
//   }
//   _nSCC = 0;
// }

// Single Source Shortest Path (Dijkstra's algorithm)

// void Graph::dijkstra(int s, double* d, int* p, list<int>* pred, list<int>*
// stack) {
//   NodeHeap<double>* h = new NodeHeap<double>(_nVert, INF, true);

//   for (int i = 0; i < _nVert; i++) {
//     d[i] = INF;
//     p[i] = 0;
//     pred[i].clear();
//   }
//   d[s] = 0;
//   p[s] = 1;
//   h->changeValue(s, 0);

//   while (h->size()) {
//     int u = h->pop();
//     list<int>::iterator v;
//     list<double>::iterator w;
//     stack->push_front(u);

//     for (v = _edges[u]->begin(), w = _weights[u]->begin(); v !=
//     _edges[u]->end(); v++, w++) {
//       if (d[*v] > d[u]+*w+PRECISION) {
// 	d[*v] = d[u]+*w;
// 	h->changeValue(*v, d[*v]);
// 	pred[*v].clear();
// 	p[*v] = 0;
//       }
//       if (d[*v] < d[u]+*w+PRECISION && d[*v] > d[u]+*w-PRECISION) {
// 	pred[*v].push_front(u);
// 	p[*v] += p[u];
//       }
//     }
//   }
//   delete h;
// }
