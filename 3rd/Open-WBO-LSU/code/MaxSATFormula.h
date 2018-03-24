/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
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

#ifndef MaxSATFormula_h
#define MaxSATFormula_h

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "FormulaPB.h"
#include "MaxTypes.h"

#include <map>
#include <string>

using NSPACE::vec;
using NSPACE::Lit;
using NSPACE::lit_Undef;
using NSPACE::mkLit;

namespace openwbo {

typedef std::map<std::string, int> nameMap;
typedef std::map<int, std::string> indexMap;

class Soft {

public:
  /*! The soft class is used to model the soft clauses in a MaxSAT formula. */
  Soft(const vec<Lit> &soft, uint64_t soft_weight, Lit assump_var,
       const vec<Lit> &relax) {
    soft.copyTo(clause);
    weight = soft_weight;
    assumption_var = assump_var;
    relax.copyTo(relaxation_vars);
  }

  Soft() {}
  ~Soft() {
    clause.clear();
    relaxation_vars.clear();
  }

  vec<Lit> clause;    //!< Soft clause
  uint64_t weight;    //!< Weight of the soft clause
  Lit assumption_var; //!< Assumption variable used for retrieving the core
  vec<Lit> relaxation_vars; //!< Relaxation variables that will be added to the
                            //! soft clause
};

class Hard {
  /*! The hard class is used to model the hard clauses in a MaxSAT formula. */
public:
  Hard(const vec<Lit> &hard) { hard.copyTo(clause); }

  Hard() {}
  ~Hard() { clause.clear(); }

  vec<Lit> clause; //!< Hard clause
};

class MaxSATFormula {
  /*! This class contains the MaxSAT formula and methods for adding soft and
   * hard clauses. */
public:
  MaxSATFormula()
      : hard_weight(UINT64_MAX), problem_type(_UNWEIGHTED_), n_vars(0),
        n_soft(0), n_hard(0), n_initial_vars(0), sum_soft_weight(0),
        max_soft_weight(0) {
    objective_function = NULL;
    format = _FORMAT_MAXSAT_;
  }

  ~MaxSATFormula() {
    for (int i = 0; i < nSoft(); i++) {
      soft_clauses[i].clause.clear();
      soft_clauses[i].relaxation_vars.clear();
    }
    soft_clauses.clear();

    for (int i = 0; i < nHard(); i++)
      hard_clauses[i].clause.clear();
    hard_clauses.clear();
  }

  MaxSATFormula *copyMaxSATFormula();

  /*! Add a new hard clause. */
  void addHardClause(vec<Lit> &lits);

  /*! Add a new soft clause. */
  void addSoftClause(uint64_t weight, vec<Lit> &lits);

  /*! Add a new soft clause with predefined relaxation variables. */
  void addSoftClause(uint64_t weight, vec<Lit> &lits, vec<Lit> &vars);

  int nVars();   // Number of variables.
  int nSoft();   // Number of soft clauses.
  int nHard();   // Number of hard clauses.
  void newVar(); // New variable.

  Lit newLiteral(bool sign = false); // Make a new literal.

  void setProblemType(int type); // Set problem type.
  int getProblemType();          // Get problem type.

  void updateSumWeights(uint64_t weight); // Update initial 'ubCost'.
  uint64_t getSumWeights() { return sum_soft_weight; }

  void setMaximumWeight(uint64_t weight); // Set initial 'currentWeight'.
  uint64_t getMaximumWeight();            // Get 'currentWeight'.

  void setHardWeight(uint64_t weight); // Set initial 'hardWeight'.
  uint64_t getHardWeight() { return hard_weight; }

  /*! Return number of initial variables. */
  int nInitialVars();

  /*! Set initial number of variables. */
  void setInitialVars(int vars);

  /*! Return i-soft clause. */
  Soft &getSoftClause(int pos);

  /*! Return i-hard clause. */
  Hard &getHardClause(int pos);

  /*! Add a new cardinality constraint. */
  void addCardinalityConstraint(Card *card);

  /*! Return i-card constraint. */
  Card *getCardinalityConstraint(int pos) {
    return cardinality_constraints[pos];
  }

  /*! Add a new PB constraint. */
  void addPBConstraint(PB *pb);

  /*! Return i-PB constraint. */
  PB *getPBConstraint(int pos) { return pb_constraints[pos]; }

  int newVarName(char *varName);
  int varID(char *varName);

  void addObjFunction(PBObjFunction *of) {
    objective_function = new PBObjFunction(of->_lits, of->_coeffs, of->_const);
  }

  PBObjFunction *getObjFunction() { return objective_function; }

  int nCard() { return cardinality_constraints.size(); }

  int nPB() { return pb_constraints.size(); }

  void convertPBtoMaxSAT();

  void setFormat(int form) { format = form; }

  int getFormat() { return format; }

  indexMap &getIndexToName() { return _indexToName; }

protected:
  // MaxSAT database
  //
  vec<Soft> soft_clauses; //<! Stores the soft clauses of the MaxSAT formula.
  vec<Hard> hard_clauses; //<! Stores the hard clauses of the MaxSAT formula.

  // PB database
  //
  PBObjFunction *objective_function;   //<! Objective function for PB.
  vec<Card *> cardinality_constraints; //<! Stores the cardinality constraints.
  vec<PB *> pb_constraints;            //<! Stores the PB constraints.

  // Properties of the MaxSAT formula
  //
  uint64_t hard_weight; //<! Weight of the hard clauses.
  int problem_type;     //<! Stores the type of the MaxSAT problem.
  int n_vars;           //<! Number of variables used in the SAT solver.
  int n_soft;           //<! Number of soft clauses.
  int n_hard;           //<! Number of hard clauses.
  int n_initial_vars;   //<! Number of variables of the initial MaxSAT formula.
  uint64_t sum_soft_weight; //<! Sum of weights of soft clauses.
  uint64_t max_soft_weight; //<! Maximum weight of soft clauses.

  // Utils for PB formulas
  //
  nameMap _nameToIndex;  //<! Map from variable name to variable id.
  indexMap _indexToName; //<! Map from variable id to variable name.

  // Format
  //
  int format;
};

} // namespace openwbo

#endif
