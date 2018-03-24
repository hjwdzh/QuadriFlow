/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
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

#include "Encodings.h"

using namespace openwbo;

// Creates an unit clause in the SAT solver
void Encodings::addUnitClause(Solver *S, Lit a, Lit blocking) {
  assert(clause.size() == 0);
  assert(a != lit_Undef);
  assert(var(a) < S->nVars());
  clause.push(a);
  if (blocking != lit_Undef)
    clause.push(blocking);
  S->addClause(clause);
  clause.clear();
}

// Creates a binary clause in the SAT solver
void Encodings::addBinaryClause(Solver *S, Lit a, Lit b, Lit blocking) {
  assert(clause.size() == 0);
  assert(a != lit_Undef && b != lit_Undef);
  assert(var(a) < S->nVars() && var(b) < S->nVars());
  clause.push(a);
  clause.push(b);
  if (blocking != lit_Undef)
    clause.push(blocking);
  S->addClause(clause);
  clause.clear();
}

// Creates a ternary clause in the SAT solver
void Encodings::addTernaryClause(Solver *S, Lit a, Lit b, Lit c, Lit blocking) {
  assert(clause.size() == 0);
  assert(a != lit_Undef && b != lit_Undef && c != lit_Undef);
  assert(var(a) < S->nVars() && var(b) < S->nVars() && var(c) < S->nVars());
  clause.push(a);
  clause.push(b);
  clause.push(c);
  if (blocking != lit_Undef)
    clause.push(blocking);
  S->addClause(clause);
  clause.clear();
}

// Creates a quaternary clause in the SAT solver
void Encodings::addQuaternaryClause(Solver *S, Lit a, Lit b, Lit c, Lit d,
                                    Lit blocking) {
  assert(clause.size() == 0);
  assert(a != lit_Undef && b != lit_Undef && c != lit_Undef && d != lit_Undef);
  assert(var(a) < S->nVars() && var(b) < S->nVars() && var(c) < S->nVars() &&
         var(d) < S->nVars());
  clause.push(a);
  clause.push(b);
  clause.push(c);
  clause.push(d);
  if (blocking != lit_Undef)
    clause.push(blocking);
  S->addClause(clause);
  clause.clear();
}
