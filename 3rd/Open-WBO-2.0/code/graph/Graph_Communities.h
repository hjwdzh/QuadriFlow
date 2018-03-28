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

#ifndef __GRAPH_COMMUNITIES__
#define __GRAPH_COMMUNITIES__

#include "Graph.h"
#include <string.h>

#include "mtl/Vec.h"

using namespace std;

namespace openwbo {

enum splitMode_ { RAND_MODE, UNFOLDING_MODE, LABEL_PROP_MODE };

class Graph_Communities {
public:
  // Constructor/Destructor:
  //
  Graph_Communities();
  ~Graph_Communities();

  int findCommunities(int mode, Graph *g);

  // Valid after findCommunities is called.
  inline int nCommunities() { return _nCommunities; }
  inline int vertexCommunity(int u) { return _vertexCommunity[u]; }
  inline double getModularity() { return _modularity; }

  inline const vec<int> &adjCommunities(int c) { return _g->vertexEdges(c); }
  inline const vec<double> &adjCommunityWeights(int c) {
    return _g->vertexWeights(c);
  }

protected:
  // Unfolding method
  Graph *nextIterationGraph();
  bool iterate();
  void computeAdjCommunities(int node);

  void resetInternalData();

  double modularity();

  void remove(int node, int comm, double dnodecomm);
  void insert(int node, int comm, double dnodecomm);

  // Label propagation method

protected:
  int _nCommunities;
  double _modularity;
  vec<int> _vertexCommunity;

  // Unfolding method
  Graph *_g; // Current working graph

  // Unfolding method - Utils
  vec<int> _vertexToComm;     // mapping of vertexes to communities
  vec<vec<int>> _communities; // vertexes of each community
  vec<double> _inside;        // total weight of edges inside each community
  vec<double> _total; // total weight of edges of vertexes in each community
                      // (inside and outside edges)

  vec<double> _adjWeight;
  vec<unsigned int> _adjComm;
  vec<bool> _adjMarked;

  vec<int> _renumber;

  // Label propagation method
};

} // namespace openwbo

#endif
