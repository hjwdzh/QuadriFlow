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

#ifndef __GRAPH__
#define __GRAPH__

#include <stdint.h>
#include <string.h>

#include "mtl/Vec.h"

using namespace std;
using NSPACE::vec;

namespace openwbo {

enum color_ { WHITE, GRAY, BLACK };

class Graph {
public:
  // Constructor/Destructor:
  //
  Graph(int nVert);
  ~Graph();

  void addEdge(int u, int v, double w = 1.0);
  int nEdges();
  void mergeDuplicatedEdges();

  // Stats
  inline int nVertexes() { return _nVert; }
  inline vec<int> &vertexEdges(int u) { return _edges[u]; }
  inline vec<double> &vertexWeights(int u) { return _weights[u]; }
  inline int nNeighbors(int u) { return _edges[u].size(); }
  inline int nIncomingEdges(int u) { return _incomingEdges[u]; }
  inline double nSelfLoops(int u) { return _nSelfLoops[u]; }

  inline double weightedDegree(int u) { return _totalWeights[u]; }
  inline double totalWeight() { return _totalWeight; }

  inline double density() {
    return ((double)nEdges()) / (((double)_nVert * _nVert));
  }
  inline double avgDegree() { return ((double)nEdges()) / ((double)_nVert); }
  inline double avgWeightedDegree() {
    return ((double)_totalWeight) / ((double)_nVert);
  }

  // Graph Operations

  /* void findAllScc(); */
  /* inline int          nSCC    ()  { return _nSCC; } */
  /* inline list<int>*   getSCC  (int i)  { return _sccs[i];  } */

  void topologicalSort(vec<int> &vertexes);

  void visitedVertexes(int u, vec<int> &reachedVertexes);

  /* void dijkstra(int s, double* d, int* p, list<int>* pred, list<int>* stack);
   */

  int connectedComponents();

  // Labels, colors and output

protected:
  /* void tarjan(int u); */

  void DFSVisit(int u, vec<int> &reachedVertexes);
  void DFSVisitIter(int u, vec<int> &reachedVertexes);

  /* void clearSCC(); */

  /* void normalizeValues(double* v, int size, double nMin, double nMax); */

protected:
  int _nVert;
  vec<vec<int>> _edges;
  vec<vec<double>> _weights;
  vec<double> _totalWeights;
  double _totalWeight;
  vec<int> _incomingEdges;
  vec<double> _nSelfLoops;

  // utils
  vec<int> _marks;
  vec<int> _markedVertexes;
  int _nMarked;

  vec<int> _index;
  vec<int> _lowlink;
  vec<int> _stack;
  int _time;

  int _nSCC;
  vec<vec<int>> _sccs;
};

} // namespace openwbo

#endif
