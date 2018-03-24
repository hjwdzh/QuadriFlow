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

#ifndef MaxTypes_h
#define MaxTypes_h

enum { _FORMAT_MAXSAT_ = 0, _FORMAT_PB_ };
enum { _VERBOSITY_MINIMAL_ = 0, _VERBOSITY_SOME_ };
enum { _UNWEIGHTED_ = 0, _WEIGHTED_ };
enum { _WEIGHT_NONE_ = 0, _WEIGHT_NORMAL_, _WEIGHT_DIVERSIFY_ };
enum {
  _ALGORITHM_WBO_ = 0,
  _ALGORITHM_LINEAR_SU_,
  _ALGORITHM_MSU3_,
  _ALGORITHM_PART_MSU3_,
  _ALGORITHM_OLL_,
  _ALGORITHM_BEST_
};
enum {
  _SATISFIABLE_ = 10,
  _UNSATISFIABLE_ = 20,
  _OPTIMUM_ = 30,
  _UNKNOWN_ = 40,
  _ERROR_ = 50
};
enum {
  _INCREMENTAL_NONE_ = 0,
  _INCREMENTAL_BLOCKING_,
  _INCREMENTAL_WEAKENING_,
  _INCREMENTAL_ITERATIVE_
};
enum { _CARD_CNETWORKS_ = 0, _CARD_TOTALIZER_, _CARD_MTOTALIZER_ };
enum { _AMO_LADDER_ = 0 };
enum { _PB_SWC_ = 0, _PB_GTE_ };
enum { _PART_SEQUENTIAL_ = 0, _PART_SEQUENTIAL_SORTED_, _PART_BINARY_ };

#endif
