/*!
 * \author Vasco Manquinho - vmm@sat.inesc-id.pt
 *
 * @section LICENSE
 *
 * Open-WBO, Copyright (c) 2013-2015, Ruben Martins, Vasco Manquinho, Ines Lynce
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

#ifndef __PB_PARSER__
#define __PB_PARSER__

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>

#include "MaxSATFormula.h"

using NSPACE::vec;
using std::cout;
using std::endl;

#define MAX_WORD_LENGTH 256

#ifndef _PB_MIN_
#define _PB_MIN_ 1
#define _PB_MAX_ 0
#endif

/*! Definition of possible constraint signs. */
enum pb_Sign { _PB_GREATER_OR_EQUAL_ = 0x1, _PB_LESS_OR_EQUAL_, _PB_EQUAL_ };

namespace openwbo {

/*! Generic parser class in open-wbo. All other parsers inherit from this class.
 */
class ParserPB {

public:
  //-------------------------------------------------------------------------
  // Constructor/destructor.
  //-------------------------------------------------------------------------

  ParserPB();
  virtual ~ParserPB();

  //-------------------------------------------------------------------------
  // Interface contract:
  //-------------------------------------------------------------------------

  virtual int parse(char *fileName);

  void parsePBFormula(char *fileName, MaxSATFormula *max) {
    maxsat_formula = max;
    parse(fileName);
  }

protected:
  // OPB instance parsing.
  virtual int parseLine();
  virtual int parseCostFunction();
  virtual int parseConstraint();
  virtual int parseProduct(int64_t *coeff, char *varName, int *varNameSize);
  virtual int getVariableID(char *varName, int varNameSize);

  inline char peek_char() { return *_fileStr; }
  inline char get_char() { return *_fileStr++; }
  inline void unget_char() { _fileStr--; }

  inline void skip_spaces() {
    while (get_char() == ' ')
      ;
    unget_char();
  }

  inline void readUntilEndOfLine() {
    static char c;
    while ((c = get_char()) != '\n' && c != '\0')
      ;
  }

  inline void parseNumber(int64_t *coeff) {
    static char word[MAX_WORD_LENGTH];
    int i = 0, c = peek_char();
    int64_t conv;

    *coeff = 1;
    while ((c == '-') || (c == '+')) {
      if (c == '-')
        *coeff *= -1;
      get_char();
      skip_spaces();
      c = peek_char();
    }
    while (isdigit(word[i] = get_char())) {
      i++;
    }
    unget_char();
    word[i] = '\0';
    assert(i > 0);

    std::istringstream ss(word);
    ss >> conv;

    // sscanf(word, "%d", &i);
    *coeff = (*coeff) * conv;
  }

  inline void parseWord(char *varName, int *varNameSize) {
    int i = 0;
    while (isgraph(varName[i] = get_char())) {
      i++;
    }
    unget_char();
    varName[i] = '\0';
    *varNameSize = i;
  }

protected:
  struct eqstr {
    bool operator()(const char *s1, const char *s2) const {
      return strcmp(s1, s2) == 0;
    }
  };

  char *_fileStr;
  int _fd;

  vec<int64_t> _coefficients;
  vec<int> _constraintVariables;

  int64_t _highestCoeffSum;

  MaxSATFormula *maxsat_formula;
};

} // namespace openwbo

#endif // __PB_PARSER__

/*****************************************************************************/
