#ifndef FLOW_H_
#define FLOW_H_
#include <vector>
#include <map>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>


using namespace boost;

typedef int EdgeWeightType;

typedef adjacency_list_traits < vecS, vecS, directedS > Traits;
typedef adjacency_list < vecS, vecS, directedS,
	property < vertex_name_t, std::string,
	property < vertex_index_t, long,
	property < vertex_color_t, boost::default_color_type,
	property < vertex_distance_t, long,
	property < vertex_predecessor_t, Traits::edge_descriptor > > > > >,

	property < edge_capacity_t, EdgeWeightType,
	property < edge_residual_capacity_t, EdgeWeightType,
	property < edge_reverse_t, Traits::edge_descriptor > > > > Graph;

void AddEdge(Traits::vertex_descriptor &v1,
	Traits::vertex_descriptor &v2,
	property_map < Graph, edge_reverse_t >::type &rev,
	const double capacity,
	Graph &g);

void AddDirectEdge(Traits::vertex_descriptor &v1,
	Traits::vertex_descriptor &v2,
	property_map < Graph, edge_reverse_t >::type &rev,
	const int capacity, const int inv_capacity,
	Graph &g);

class MaxFlowHelper
{
public:
	MaxFlowHelper()
	{
		rev = get(edge_reverse, g);
		num = 0;
	}
	int num;
	void resize(int n) {
		vertex_descriptors.resize(n);
		for (int i = 0; i < n; ++i) {
			vertex_descriptors[i] = add_vertex(g);
		}
		num = n;
	}
	int compute()
	{
		EdgeWeightType flow = boykov_kolmogorov_max_flow(g, vertex_descriptors.front(), vertex_descriptors.back());
		return flow;
	}
	void AddEdge(int x, int y, int c, int rc) {
		AddDirectEdge(vertex_descriptors[x], vertex_descriptors[y], rev, c, rc, g);
	}
	void Apply(std::unordered_map<int64_t, std::pair<int, int> >& edge_to_variable, std::vector<Vector2i>& edge_diff)
	{
		property_map<Graph, edge_capacity_t>::type
			capacity = get(edge_capacity, g);
		property_map<Graph, edge_residual_capacity_t>::type
			residual_capacity = get(edge_residual_capacity, g);


		graph_traits<Graph>::vertex_iterator u_iter, u_end;
		graph_traits<Graph>::out_edge_iterator ei, e_end;
		for (tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter)
			for (tie(ei, e_end) = out_edges(*u_iter, g); ei != e_end; ++ei)
				if (capacity[*ei] > 0) {
					int flow = (capacity[*ei] - residual_capacity[*ei]);
					if (flow > 0) {
						int64_t key = (int64_t)*u_iter * num + target(*ei, g);
						auto q = edge_to_variable[key];
						edge_diff[q.first / 2][q.first % 2] += q.second * flow;
					}
				}
	}
	Graph g;
	property_map < Graph, edge_reverse_t >::type rev;
	std::vector<Traits::vertex_descriptor> vertex_descriptors;

};
int flow(std::vector<std::map<int, std::pair<int, int> > >& graph);


inline void AddEdge(Traits::vertex_descriptor &v1, Traits::vertex_descriptor &v2, property_map < Graph, edge_reverse_t >::type &rev, const double capacity, Graph &g)
{
	Traits::edge_descriptor e1 = add_edge(v1, v2, g).first;
	Traits::edge_descriptor e2 = add_edge(v2, v1, g).first;
	put(edge_capacity, g, e1, capacity);
	put(edge_capacity, g, e2, capacity);

	rev[e1] = e2;
	rev[e2] = e1;
}

inline void AddDirectEdge(Traits::vertex_descriptor &v1, Traits::vertex_descriptor &v2, property_map < Graph, edge_reverse_t >::type &rev, const int capacity,
	const int inv_capacity, Graph &g)
{
	Traits::edge_descriptor e1 = add_edge(v1, v2, g).first;
	Traits::edge_descriptor e2 = add_edge(v2, v1, g).first;
	put(edge_capacity, g, e1, capacity);
	put(edge_capacity, g, e2, inv_capacity);

	rev[e1] = e2;
	rev[e2] = e1;
}
#endif
