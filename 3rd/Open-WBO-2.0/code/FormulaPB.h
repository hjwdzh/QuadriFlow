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

#ifndef FormulaPB_h
#define FormulaPB_h

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

using NSPACE::vec;
using NSPACE::Lit;

namespace openwbo {

// Cardinality constraint of the form atMostK
class Card {

public:
  Card(vec<Lit> &lits, int64_t rhs, bool sign = false) {
    lits.copyTo(_lits);
    _rhs = rhs;
    if (sign) {
      int s = 0;
      for (int i = 0; i < _lits.size(); i++) {
        s += 1;
        _lits[i] = ~_lits[i];
      }
      _rhs = s - _rhs;
    }
  }

  Card() { _rhs = 0; }
  ~Card() {}

  void print() {
    printf("Card: ");

    for (int i = 0; i < _lits.size(); i++) {
      if (sign(_lits[i]))
        printf("~");
      printf("%d ", var(_lits[i]) + 1);
    }
    printf(" <= %d\n", (int)_rhs);
  }

  vec<Lit> _lits;
  int64_t _rhs;
};

// PB constraint. The constraint sign is encoded in the structure.
class PB {

public:
  PB(vec<Lit> &lits, vec<uint64_t> &coeffs, int64_t rhs, bool s = false) {
    lits.copyTo(_lits);
    coeffs.copyTo(_coeffs);
    _rhs = rhs;
    _sign = s;
  }

  PB() {
    _rhs = 0;
    _sign = false;
  }
  ~PB() {}

  void addProduct(Lit l, int64_t c) {
    _coeffs.push();
    _lits.push();
    if (c >= 0) {
      _coeffs[_coeffs.size() - 1] = c;
      _lits[_lits.size() - 1] = l;
    } else {
      _coeffs[_coeffs.size() - 1] = -c;
      _lits[_lits.size() - 1] = ~l;
      _rhs += -c;
    }
  }

  void addRHS(int64_t rhs) { _rhs += rhs; }

  void changeSign() {
    int s = 0;
    for (int i = 0; i < _coeffs.size(); i++) {
      s += _coeffs[i];
      _lits[i] = ~(_lits[i]);
    }
    _rhs = s - _rhs;
    _sign = !(_sign);
  }

  bool isClause() {
    // Assume _sign == false...
    bool sign = _sign;
    if (_sign)
      changeSign();
    if (_rhs != 1) {
      if (_sign != sign)
        changeSign();
      return false;
    }
    for (int i = 0; i < _coeffs.size(); i++) {
      if (_coeffs[i] != 1) {
        if (_sign != sign)
          changeSign();
        return false;
      }
    }
    return true;
  }

  bool isCardinality() {
    // Assume _sign == false...
    bool sign = _sign;
    if (_sign)
      changeSign();
    for (int i = 0; i < _coeffs.size(); i++) {
      if (_coeffs[i] != 1) {
        if (_sign != sign)
          changeSign();
        return false;
      }
    }
    return true;
  }

  void print() {
    // Assume _sign == false...
    if (isClause())
      printf("Clause: ");
    else if (isCardinality())
      printf("Card: ");
    else
      printf("PB: ");

    for (int i = 0; i < _coeffs.size(); i++) {
      printf("%d ", (int)_coeffs[i]);
      if (sign(_lits[i]))
        printf("~");
      printf("%d ", var(_lits[i]) + 1);
    }
    printf(" >= %d\n", (int)_rhs);
  }

  vec<uint64_t> _coeffs;
  vec<Lit> _lits;
  int64_t _rhs;
  bool _sign; // atLeast: false; atMost: true
};

class PBObjFunction {

public:
  PBObjFunction(vec<Lit> &lits, vec<uint64_t> &coeffs, int64_t c = 0) {
    lits.copyTo(_lits);
    coeffs.copyTo(_coeffs);
    _const = c;
  }

  PBObjFunction() { _const = 0; }
  ~PBObjFunction() {}

  void addProduct(Lit l, int64_t c) {
    _coeffs.push();
    _lits.push();
    if (c >= 0) {
      _coeffs[_coeffs.size() - 1] = c;
      _lits[_lits.size() - 1] = l;
    } else {
      _coeffs[_coeffs.size() - 1] = -c;
      _lits[_coeffs.size() - 1] = ~l;
      _const += c;
    }
  }

  vec<uint64_t> _coeffs;
  vec<Lit> _lits;
  int64_t _const;
};

} // namespace openwbo

#endif
