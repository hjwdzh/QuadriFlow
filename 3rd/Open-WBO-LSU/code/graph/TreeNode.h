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

#ifndef TREE_NODE_H_
#define TREE_NODE_H_

#include "../Encoder.h"
#include <deque>

using NSPACE::vec;
using NSPACE::Lit;

namespace openwbo {

class TreeNode {
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
  inline bool hasEncodingAssumptions() { return encoding_assumptions != NULL; }

  inline void incrementLowerBound(int64_t inc = 1) { lb += inc; }
};
} // namespace openwbo

#endif
