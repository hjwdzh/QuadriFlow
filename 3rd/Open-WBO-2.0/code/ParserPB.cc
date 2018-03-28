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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <unistd.h>

#include "ParserPB.h"

using namespace openwbo;

//-------------------------------------------------------------------------
// Constructor/destructor.
//-------------------------------------------------------------------------

ParserPB::ParserPB() : _highestCoeffSum(0) {}

ParserPB::~ParserPB() {}

//! Parse an input file and loads it into the corresponding data structure.

int ParserPB::parse(char *fileName) {
  _highestCoeffSum = 0;

  if ((_fd = open(fileName, O_RDONLY)) < 0) {
    printf("c Error: Unable to open input stream for file %s\n", fileName);
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }
  cout << "c Instance file " << fileName << endl;

  struct stat statbuf;
  if (fstat(_fd, &statbuf) < 0) {
    cout << "c Error: Unable to get size of file " << fileName << endl;
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }
  cout << "c File size is " << statbuf.st_size << " bytes." << endl;

  if ((_fileStr = (char *)mmap(0, statbuf.st_size, PROT_READ, MAP_SHARED, _fd,
                               0)) == (caddr_t)-1) {
    cout << "c Error: Unable to put in memory file " << fileName << endl;
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }
  char *startMap = _fileStr;

  // vmm - Check why I need to do this!!!!
  if (strlen(_fileStr) != (unsigned)statbuf.st_size)
    _fileStr[statbuf.st_size] = '\0';

  int line = 0;
  while (peek_char() != '\0') {
    int error = parseLine();

    if (error != 0) {
      cout << "c Error: Parse Error " << error << " in line " << ++line << endl;
      printf("s UNKNOWN\n");
      exit(_ERROR_);
    }
    line++;
  }

  // Clear memory map of input file.
  munmap(startMap, statbuf.st_size);
  close(_fd);

  // cout << "c Highest Coefficient sum: " << _highestCoeffSum << endl;

  return 0;
}

//-------------------------------------------------------------------------
// PROTECTED
//-------------------------------------------------------------------------

//! Parse an input file line

int ParserPB::parseLine() {
  // static int nLine = 1;
  skip_spaces();
  char c = peek_char();

  // int size = strlen(_fileStr);
  // if (size < 100) printf("%s\n", _fileStr);
  // printf("c Parsing line #%d %d\n", nLine++, strlen(_fileStr));

  if (c == '*' || c == '\0' || c == 10 || c == '\n') {
    // Line is empty or end of line...
    readUntilEndOfLine();
    return 0;
  } else if (c == 'm') {
    return parseCostFunction();
  } else {
    // Line must represent a pseudo-boolean constraint
    return parseConstraint();
  }
}

//! Parse an input file line corresponding to the cost function.
/*!
  \return Returns _PB_PARSER_NO_ERROR_ if the line corresponds to the
  cost function and was correctly parsed.
  Otherwise, returns the respective code error.
*/

int ParserPB::parseCostFunction() {
  // int objective = _PB_MIN_;
  static char word[MAX_WORD_LENGTH];
  int i;

  // printf("c Parsing objective function...\n");

  parseWord(word, &i);

  // if (strncmp("max:", word, 4) == 0)
  // objective = _PB_MAX_;
  // Currently only supports min functions
  if (strncmp("min:", word, 4) != 0) {
    // Not a valid cost function
    cout << "c Error: Invalid objective function " << endl;
    cout << "s UNKNOWN" << endl;
    exit(_ERROR_);
  }

  int64_t coeff;
  char varName[MAX_WORD_LENGTH], c;
  int varNameSize;
  // int64_t coeffSum = 0;
  PBObjFunction *of = new PBObjFunction();

  skip_spaces();
  c = peek_char();
  if (c == '\n' || c == '\0' || c == 10 || c == 13) {
    while (c != '\0' && (c == 10 || c == 13 || c == '\n')) {
      get_char();
      c = peek_char();
    }
  }

  do {
    parseProduct(&coeff, varName, &varNameSize);

    int varID = getVariableID(varName, varNameSize);

    // cout << "c CF Product: " << coeff << " " << varName << " "
    // 	 << varNameSize << " " << varID << endl;
    of->addProduct(mkLit(varID), coeff);

    skip_spaces();
    c = peek_char();
    if (c == ';') {
      get_char();
      c = peek_char();
    }
  } while (c != '\0' && c != 10 && c != 13 && c != '\n');

  while (c != '\0' && (c == 10 || c == 13 || c == '\n')) {
    get_char();
    c = peek_char();
  }

  maxsat_formula->addObjFunction(of);

  delete of;

  return 0;
}

//! Parse a product that consists in a coefficient and a literal (variable
// and sign).
/*!
  \param coeff Reference to the coefficient
  \param varName Reference to the string containing the variable name
  \param varNameSize Reference to the number of characters in the variable name
  \param sign Reference to the literal sign
  \return Returns _PB_PARSER_NO_ERROR_ if the product was correctly parsed.
  Otherwise, returns the respective code error.
*/

int ParserPB::parseProduct(int64_t *coeff, char *varName, int *varNameSize) {

  skip_spaces();
  parseNumber(coeff);
  skip_spaces();
  if (peek_char() == '*')
    get_char(); // To allow for '*' between coefficient and variable name
  skip_spaces();

  parseWord(varName, varNameSize);
  if (varName[(*varNameSize) - 1] == ';') {
    // Removes possible ; from variable name
    (*varNameSize)--;
    varName[*varNameSize] = '\0';
  }

  return 0;
}

//! Parse a line corresponding to a pseudo-Boolean constraint
/*!
  \return Returns _PB_PARSER_NO_ERROR_ if the constraint was correctly parsed.
  Otherwise, returns the respective code error.
*/

int ParserPB::parseConstraint() {
  int64_t coeff;
  char varName[MAX_WORD_LENGTH], c;
  int varNameSize;
  PB *p = new PB();

  // printf("c Parsing Constraint...\n");

  // Read all products
  do {
    parseProduct(&coeff, varName, &varNameSize);
    int varID = getVariableID(varName, varNameSize);

    p->addProduct(mkLit(varID), coeff);

    // cout << coeff << " " << varName << " " << " [ " << varID << " ] " <<
    // endl;

    skip_spaces();
    c = peek_char();

    if (c == '\0' || c == 10 || c == 13 || c == '\n') {
      // At the end of the line and no sign was found!!!
      printf("c Error: end of constraint line without sign\n");
      printf("s UNKNOWN\n");
      exit(_ERROR_);
    }
  } while (c != '<' && c != '>' && c != '=');

  // printf("p sign %d\n",p->_sign);

  // Read constraint sign
  pb_Sign ctrSign = _PB_GREATER_OR_EQUAL_;
  if (c == '=')
    ctrSign = _PB_EQUAL_;
  else if (c == '<') {
    ctrSign = _PB_LESS_OR_EQUAL_;
    p->_sign = true;
  }

  get_char();
  c = peek_char();

  if (ctrSign != _PB_EQUAL_ && c != '=') {
    printf("c Error: invalid constraint sign.\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  } else if (ctrSign != _PB_EQUAL_)
    get_char();

  skip_spaces();

  // Read constraint rhs
  // int64_t rhs;
  parseNumber(&coeff);
  p->addRHS(coeff);
  if (ctrSign == _PB_LESS_OR_EQUAL_) {
    p->changeSign();
  }

  readUntilEndOfLine();

  if (ctrSign == _PB_EQUAL_) {
    PB *p2 = new PB(p->_lits, p->_coeffs, p->_rhs, true);
    assert(p->_sign == false);
    maxsat_formula->addPBConstraint(p);
    maxsat_formula->addPBConstraint(p2);
    delete p2;
  } else
    maxsat_formula->addPBConstraint(p);

  delete p;

  return 0;
}

//! Get the variable identifier corresponding to a given name. If the
// variable does not exist, a new identifier is created.

int ParserPB::getVariableID(char *varName, int varNameSize) {
  int id = maxsat_formula->varID(varName);
  if (id == var_Undef)
    id = maxsat_formula->newVarName(varName);
  return id;
}

/*****************************************************************************/
