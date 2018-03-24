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

#include "Enc_Ladder.h"

using namespace openwbo;

/*_________________________________________________________________________________________________
  |
  |  encode : (S : Solver *) (lits : vec<Lit>&)  ->  [void]
  |
  |  Description:
  |
  |    Encodes that exactly one literal from 'lits' is assigned value true.
  |    Uses the Ladder/Regular encoding for translating the AMO constraint into
  |    CNF.
  |
  |  For further details see:
  |    * Carlos Ansótegui, Felip Manyà: Mapping Problems with Finite-Domain
  |      Variables into Problems with Boolean Variables. SAT 2004
  |    * Ian Gent and Peter Nightingale. A New Encoding of All Different into
  |      SAT. ModRef 2004
  |
  |  Pre-conditions:
  |    * Assumes that 'lits' is not empty.
  |
  |  Post-conditions:
  |    * 'S' is updated with the clauses that encode the AMO constraint.
  |
  |________________________________________________________________________________________________@*/
void Ladder::encode(Solver *S, vec<Lit> &lits) {

  assert(lits.size() != 0);

  if (lits.size() == 1) {
    addUnitClause(S, lits[0]);
  } else {

    vec<Lit> seq_auxiliary;

    for (int i = 0; i < lits.size() - 1; i++) {
      seq_auxiliary.push(mkLit(S->nVars(), false));
      newSATVariable(S);
    }

    for (int i = 0; i < lits.size(); i++) {
      if (i == 0) {
        // With the clause below it becomes EO encoding.
        // addBinaryClause(S, lits[i], ~seq_auxiliary[i]);
        addBinaryClause(S, ~lits[i], seq_auxiliary[i]);
      } else if (i == lits.size() - 1) {
        addBinaryClause(S, lits[i], seq_auxiliary[i - 1]);
        addBinaryClause(S, ~lits[i], ~seq_auxiliary[i - 1]);
      } else {
        addBinaryClause(S, ~seq_auxiliary[i - 1], seq_auxiliary[i]);
        addTernaryClause(S, lits[i], ~seq_auxiliary[i], seq_auxiliary[i - 1]);
        addBinaryClause(S, ~lits[i], seq_auxiliary[i]);
        addBinaryClause(S, ~lits[i], ~seq_auxiliary[i - 1]);
      }
    }
  }
}