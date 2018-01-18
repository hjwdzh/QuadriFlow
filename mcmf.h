// Cristinel Ababei, January 2009, Fargo ND
// E-mail: cristinel.ababei@ndsu.edu
//
// This is a C++ implementation of CS2 min-cost-max-flow scaling algorithm. 
//
// This is intended to be one of the cleanest and simplest to use minimum-cost 
// max-flow (MCMF) implementation using C++.  If you have a C++ application in 
// which you need to use a MCMF algo, then this may be your most elegant bet.
// See main() function for an example of how to use it.
// 
// This is an adapted (i.e., ported to C++) version of the faimous CS2 algo;
// CS2 is the second version of scaling algorithm for minimum-cost max-flow 
// problems.  For a detailed description of this famous algo, see:
// A.V. Goldberg, "An Efficient Implementation of a Scaling Minimum-Cost 
// Flow Algorithm", Journal of Algorithms, Vol. 22, pp. 1-29, 1997.

#ifndef _MCMF_H_
#define _MCMF_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//
// MCMF_CS2
//
////////////////////////////////////////////////////////////////////////////////

#define MAX_64 (0x7fffffffffffffffLL)
#define MAX_32 (0x7fffffff)
#define PRICE_MAX MAX_64

#define MAXLINE       100 // max line length in the input file
#define ARC_FIELDS      5 // no of fields in arc line
#define NODE_FIELDS     2 // no of fields in node line
#define P_FIELDS        3 // no of fields in problem line
#define PROBLEM_TYPE "min" //  name of problem type

#define UNFEASIBLE          2
#define ALLOCATION_FAULT    5
#define PRICE_OFL           6

// parameters
#define UPDT_FREQ      0.4
#define UPDT_FREQ_S    30
#define SCALE_DEFAULT  12.0
// PRICE_OUT_START may not be less than 1
#define PRICE_OUT_START  1
#define CUT_OFF_POWER    0.44
#define CUT_OFF_COEF     1.5
#define CUT_OFF_POWER2   0.75
#define CUT_OFF_COEF2    1
#define CUT_OFF_GAP      0.8
#define CUT_OFF_MIN      12
#define CUT_OFF_INCREASE 4

#define TIME_FOR_PRICE_IN1    2
#define TIME_FOR_PRICE_IN2    4
#define TIME_FOR_PRICE_IN3    6

#define MAX_CYCLES_CANCELLED 0
#define START_CYCLE_CANCEL   100

#define MAX( x, y ) ( ( (x) > (y) ) ?  x : y )
#define MIN( x, y ) ( ( (x) < (y) ) ? x : y )
#define ABS( x ) ( (x) >= 0 ) ? (x) : -(x)
#define N_NODE( i ) ( ( (i) == NULL ) ? -1 : ( (i) - _nodes + _node_min ) )
#define N_ARC( a ) ( ( (a) == NULL )? -1 : (a) - _arcs )

class MCMF_CS2
{
public:
	typedef long long int excess_t;
	typedef long long int price_t;

	class NODE;

	class ARC {
	public:
		long _rez_capacity; // residual capacity;
		price_t _cost; // cost of arc;
		NODE *_head; // head node;
		ARC *_sister; // opposite arc;
	public:
		ARC() {}
		~ARC() {}

		void set_rez_capacity(long rez_capacity) { _rez_capacity = rez_capacity; }
		void dec_rez_capacity(long delta) { _rez_capacity -= delta; }
		void inc_rez_capacity(long delta) { _rez_capacity += delta; }
		void set_cost(price_t cost) { _cost = cost; }
		void multiply_cost(price_t mult) { _cost *= mult; }
		void set_head(NODE *head) { _head = head; }
		void set_sister(ARC *sister) { _sister = sister; }
		long rez_capacity() { return _rez_capacity; }
		price_t cost() { return _cost; }
		NODE *head() { return _head; }
		ARC *sister() { return _sister; }
	};

	class NODE {
	public:
		excess_t _excess; // excess of the node;
		price_t _price; // distance from a sink;
		ARC *_first; // first outgoing arc;
		ARC *_current; // current outgoing arc;
		ARC *_suspended;
		NODE *_q_next; // next node in push queue
		NODE *_b_next; // next node in bucket-list;
		NODE *_b_prev; // previous node in bucket-list;
		long _rank; // bucket number;
		long _inp; // auxilary field;
	public:
		NODE() {}
		~NODE() {}

		void set_excess(excess_t excess) { _excess = excess; }
		void dec_excess(long delta) { _excess -= delta; }
		void inc_excess(long delta) { _excess += delta; }
		void set_price(price_t price) { _price = price; }
		void dec_price(long delta) { _price -= delta; }
		void set_first(ARC *first) { _first = first; }
		void set_current(ARC *current) { _current = current; }
		void inc_current() { _current++; }
		void set_suspended(ARC *suspended) { _suspended = suspended; }
		void set_q_next(NODE *q_next) { _q_next = q_next; }
		void set_b_next(NODE *b_next) { _b_next = b_next; }
		void set_b_prev(NODE *b_prev) { _b_prev = b_prev; }
		void set_rank(long rank) { _rank = rank; }
		void set_inp(long inp) { _inp = inp; }
		excess_t excess() { return _excess; }
		price_t price() { return _price; }
		ARC *first() { return _first; }
		void dec_first() { _first--; }
		void inc_first() { _first++; }
		ARC *current() { return _current; }
		ARC *suspended() { return _suspended; }
		NODE *q_next() { return _q_next; }
		NODE *b_next() { return _b_next; }
		NODE *b_prev() { return _b_prev; }
		long rank() { return _rank; }
		long inp() { return _inp; }
	};

	class BUCKET {
	private:
		// 1st node with positive excess or simply 1st node in the buket;
		NODE *_p_first;
	public:
		BUCKET(NODE *p_first) : _p_first(p_first) {}
		BUCKET() {}
		~BUCKET() {}

		void set_p_first(NODE *p_first) { _p_first = p_first; }
		NODE *p_first() { return _p_first; }
	};

public:
	long _n; // number of nodes
	long _m; // number of arcs

	long *_cap; // array containig capacities
	NODE *_nodes; // array of nodes
	NODE *_sentinel_node; // next after last
	NODE *_excq_first; // first node in push-queue
	NODE *_excq_last; // last node in push-queue
	ARC *_arcs; // array of arcs
	ARC *_sentinel_arc; // next after last

	BUCKET *_buckets; // array of buckets
	BUCKET *_l_bucket; // last bucket
	long _linf; // number of l_bucket + 1
	int _time_for_price_in;

	price_t _epsilon; // quality measure
	price_t _dn; // cost multiplier = number of nodes + 1
	price_t _price_min; // lowest bound for prices
	price_t _mmc; // multiplied maximal cost
	double _f_scale; // scale factor
	double _cut_off_factor; // multiplier to produce cut_on and cut_off from n and epsilon
	double _cut_on; // the bound for returning suspended arcs
	double _cut_off; // the bound for suspending arcs
	excess_t _total_excess; // total excess

	// if = 1 - signal to start price-in ASAP - 
	// maybe there is infeasibility because of susoended arcs 
	int _flag_price;
	// if = 1 - update failed some sources are unreachable: either the 
	// problem is unfeasible or you have to return suspended arcs
	int _flag_updt;
	// maximal number of cycles cancelled during price refine 
	int _snc_max;

	// dummy variables;
	ARC _d_arc; // dummy arc - for technical reasons
	NODE _d_node; // dummy node - for technical reasons
	NODE *_dummy_node; // the address of d_node
	NODE *_dnode;

	long _n_rel; // number of relabels from last price update
	long _n_ref; // current number of refines
	long _n_src; // current number of nodes with excess
	long _n_push;
	long _n_relabel;
	long _n_discharge;
	long _n_refine;
	long _n_update;
	long _n_scan;
	long _n_prscan;
	long _n_prscan1;
	long _n_prscan2;
	long _n_bad_pricein;
	long _n_bad_relabel;
	long _n_prefine;

	bool _no_zero_cycles; // finds an optimal flow with no zero-cost cycles
	bool _check_solution; // check feasibility/optimality. HIGH OVERHEAD!
	bool _comp_duals; // compute prices?
	bool _cost_restart; // to be able to restart after a cost function change
	bool _print_ans;
	long long int *_node_balance;

	// sketch variables used during reading in arcs;
	long _node_min; // minimal no of nodes
	long _node_max; // maximal no of nodes
	long *_arc_first; // internal array for holding
	// - node degree
	// - position of the first outgoing arc
	long *_arc_tail; // internal array: tails of the arcs
	long _pos_current;
	ARC *_arc_current;
	ARC *_arc_new;
	ARC *_arc_tmp;
	price_t _max_cost; // maximum cost
	excess_t _total_p; // total supply
	excess_t _total_n; // total demand
	// pointers to the node structure
	NODE *_i_node;
	NODE *_j_node;

public:
	MCMF_CS2(long num_nodes, long num_arcs) {
		_n = num_nodes;
		_m = num_arcs;

		_flag_price = 0;
		_flag_updt = 0;
		_n_push = 0;
		_n_relabel = 0;
		_n_discharge = 0;
		_n_refine = 0;
		_n_update = 0;
		_n_scan = 0;
		_n_prscan = 0;
		_n_prscan1 = 0;
		_n_prscan2 = 0;
		_n_bad_pricein = 0;
		_n_bad_relabel = 0;
		_n_prefine = 0;
		_no_zero_cycles = false;
		_check_solution = false;
		_comp_duals = false;
		_cost_restart = false;
		_print_ans = true;
		// allocate arrays and prepare for "receiving" arcs;
		// will also reset _pos_current, etc.;
		allocate_arrays();
	}
	~MCMF_CS2() {}

	void err_end(int cc);
	void allocate_arrays();
	void deallocate_arrays();
	void set_arc(long tail_node_id, long head_node_id,
		long low_bound, long up_bound, price_t cost);
	void set_supply_demand_of_node(long id, long excess);
	void pre_processing();
	void cs2_initialize();
	void up_node_scan(NODE *i);
	void price_update();
	int relabel(NODE *i);
	void discharge(NODE *i);
	int price_in();
	void refine();
	int price_refine();
	void compute_prices();
	void price_out();
	int update_epsilon();
	int check_feas();
	int check_cs();
	int check_eps_opt();
	void init_solution();
	void cs_cost_reinit();
	void cs2_cost_restart(double *objective_cost);
	void print_solution();
	void print_graph();
	void finishup(double *objective_cost);
	void cs2(double *objective_cost);
	int run_cs2();

	// shared utils;
	void increase_flow(NODE *i, NODE *j, ARC *a, long df) {
		i->dec_excess(df);
		j->inc_excess(df);
		a->dec_rez_capacity(df);
		a->sister()->inc_rez_capacity(df);
	}
	bool time_for_update() {
		return (_n_rel > _n * UPDT_FREQ + _n_src * UPDT_FREQ_S);
	}
	// utils for excess queue;
	void reset_excess_q() {
		for (; _excq_first != NULL; _excq_first = _excq_last) {
			_excq_last = _excq_first->q_next();
			_excq_first->set_q_next(_sentinel_node);
		}
	}
	bool out_of_excess_q(NODE *i) { return (i->q_next() == _sentinel_node); }
	bool empty_excess_q() { return (_excq_first == NULL); }
	bool nonempty_excess_q() { return (_excq_first != NULL); }
	void insert_to_excess_q(NODE *i) {
		if (nonempty_excess_q()) {
			_excq_last->set_q_next(i);
		}
		else {
			_excq_first = i;
		}
		i->set_q_next(NULL);
		_excq_last = i;
	}
	void insert_to_front_excess_q(NODE *i) {
		if (empty_excess_q()) {
			_excq_last = i;
		}
		i->set_q_next(_excq_first);
		_excq_first = i;
	}
	void remove_from_excess_q(NODE *i) {
		i = _excq_first;
		_excq_first = i->q_next();
		i->set_q_next(_sentinel_node);
	}
	// utils for excess queue as a stack;
	bool empty_stackq() { return empty_excess_q(); }
	bool nonempty_stackq() { return nonempty_excess_q(); }
	void reset_stackq() { reset_excess_q(); }
	void stackq_push(NODE *i) {
		i->set_q_next(_excq_first);
		_excq_first = i;
	}
	void stackq_pop(NODE *i) {
		remove_from_excess_q(i);
	}
	// utils for buckets;
	void reset_bucket(BUCKET *b) { b->set_p_first(_dnode); }
	bool nonempty_bucket(BUCKET *b) { return ((b->p_first()) != _dnode); }
	void insert_to_bucket(NODE *i, BUCKET *b) {
		i->set_b_next(b->p_first());
		b->p_first()->set_b_prev(i);
		b->set_p_first(i);
	}
	void get_from_bucket(NODE *i, BUCKET *b) {
		i = b->p_first();
		b->set_p_first(i->b_next());
	}
	void remove_from_bucket(NODE *i, BUCKET *b) {
		if (i == b->p_first()) {
			b->set_p_first(i->b_next());
		}
		else {
			i->b_prev()->set_b_next(i->b_next());
			i->b_next()->set_b_prev(i->b_prev());
		}
	}
	// misc utils;
	void update_cut_off() {
		if (_n_bad_pricein + _n_bad_relabel == 0) {
			_cut_off_factor = CUT_OFF_COEF2 * pow((double)_n, CUT_OFF_POWER2);
			_cut_off_factor = MAX(_cut_off_factor, CUT_OFF_MIN);
			_cut_off = _cut_off_factor * _epsilon;
			_cut_on = _cut_off * CUT_OFF_GAP;
		}
		else {
			_cut_off_factor *= CUT_OFF_INCREASE;
			_cut_off = _cut_off_factor * _epsilon;
			_cut_on = _cut_off * CUT_OFF_GAP;
		}
	}
	void exchange(ARC *a, ARC *b) {
		if (a != b) {
			ARC *sa = a->sister();
			ARC *sb = b->sister();
			long d_cap;

			_d_arc.set_rez_capacity(a->rez_capacity());
			_d_arc.set_cost(a->cost());
			_d_arc.set_head(a->head());

			a->set_rez_capacity(b->rez_capacity());
			a->set_cost(b->cost());
			a->set_head(b->head());

			b->set_rez_capacity(_d_arc.rez_capacity());
			b->set_cost(_d_arc.cost());
			b->set_head(_d_arc.head());

			if (a != sb) {
				b->set_sister(sa);
				a->set_sister(sb);
				sa->set_sister(b);
				sb->set_sister(a);
			}

			d_cap = _cap[a - _arcs];
			_cap[a - _arcs] = _cap[b - _arcs];
			_cap[b - _arcs] = d_cap;
		}
	}
};

#endif