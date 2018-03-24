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

#include <stdlib.h>

#ifndef __HEAP__
#define __HEAP__

using namespace std;

template <class T> class NodeHeap {
public:
  // Constructor/Destructor:
  //
  NodeHeap(int nodes, T inf, bool min = true);
  ~NodeHeap();

  inline int top() { return _heap[0]; }
  inline int size() { return _size; }
  inline T value(int node) { return _values[node]; }

  void changeValue(int node, T value);
  int pop();

protected:
  void siftDown(int i);
  void siftUp(int i);
  void exchange(int i, int j);

  inline int parent(int i) { return (i - 1) / 2; }
  inline int left(int i) { return 2 * i + 1; }
  inline int right(int i) { return 2 * (i + 1); }

protected:
  bool _type;
  int _size;
  int *_heap;
  int *_location;
  T *_values;
  T _inf;
};

template <class T> NodeHeap<T>::NodeHeap(int nodes, T inf, bool min) {
  _size = 0;
  _inf = inf;
  _type = min;
  _heap = new int[nodes];
  _location = new int[nodes];
  _values = new T[nodes];
  for (int i = 0; i < nodes; i++) {
    _values[i] = _inf;
    _location[i] = -1;
  }
}

template <class T> NodeHeap<T>::~NodeHeap() {
  delete[] _heap;
  delete[] _values;
  delete[] _location;
}

template <class T> int NodeHeap<T>::pop() {
  if (_size == 0)
    return -1;

  int root = _heap[0];
  _heap[0] = _heap[--_size];
  _location[_heap[0]] = 0;
  _location[root] = -1;
  siftDown(0);
  return root;
}

template <class T> void NodeHeap<T>::changeValue(int node, T value) {
  T prev = _values[node];
  _values[node] = value;
  if (_location[node] == -1) {
    // Insert node
    _heap[_size] = node;
    _location[node] = _size;
    _size++;
    siftUp(_size - 1);
  } else if ((_type && prev > value) || (!_type && prev < value))
    siftUp(_location[node]);
  else if ((_type && prev < value) || (!_type && prev > value))
    siftDown(_location[node]);
}

template <class T> void NodeHeap<T>::siftDown(int i) {
  int l = left(i), r = right(i);
  int select = i;

  if (l < _size) {
    if ((_type && _values[_heap[l]] < _values[_heap[select]]) ||
        (!_type && _values[_heap[l]] > _values[_heap[select]]))
      select = l;
  }
  if (r < _size) {
    if ((_type && _values[_heap[r]] < _values[_heap[select]]) ||
        (!_type && _values[_heap[r]] > _values[_heap[select]]))
      select = r;
  }
  if (select != i) {
    exchange(i, select);
    siftDown(select);
  }
}

template <class T> void NodeHeap<T>::siftUp(int i) {
  int j = parent(i);
  if (j >= 0) {
    if ((_type && _values[_heap[j]] > _values[_heap[i]]) ||
        (!_type && _values[_heap[j]] < _values[_heap[i]])) {
      exchange(i, j);
      siftUp(j);
    }
  }
}

template <class T> void NodeHeap<T>::exchange(int i, int j) {
  int aux = _heap[i];
  _heap[i] = _heap[j];
  _heap[j] = aux;
  _location[_heap[i]] = i;
  _location[_heap[j]] = j;
}

#endif
