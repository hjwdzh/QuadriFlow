#ifndef FLOW_H_
#define FLOW_H_
#include <vector>
#include <map>
#include <list>
#include <Eigen/Core>

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

class AdvanceMaxFlowHelper
{
public:
    struct FlowInfo
    {
        int id;
        int capacity, flow;
    };
    AdvanceMaxFlowHelper()
    {
        num = 0;
    }
    int num;
    std::vector<FlowInfo*> variable_to_edge;
    void resize(int n, int m) {
        graph.resize(n);
        variable_to_edge.resize(m, 0);
        num = n;
    }
    void AddEdge(int x, int y, int c, int rc, int v) {
        FlowInfo flow;
        flow.id = y;
        flow.capacity = c;
        flow.flow = 0;
        graph[x].push_back(flow);
        flow.id = x;
        flow.capacity = rc;
        flow.flow = 0;
        graph[y].push_back(flow);
        if (v != -1)
            variable_to_edge[v] = &graph[y].back();
    }
    bool CheckFlip(std::vector<Vector2i>& edge_diff, std::vector<Vector3i>& face_edgeIds, std::vector<Vector3i>& face_edgeOrients, int fid)
    {
        if (fid < 0)
            return true;
        Vector2i diff[3];
        for (int j = 0; j < 3; ++j) {
            int eid = face_edgeIds[fid][j];
            Vector2i d = edge_diff[eid];
            if (variable_to_edge[eid * 2])
                d[0] += variable_to_edge[eid * 2]->flow;
            if (variable_to_edge[eid * 2 + 1])
                d[1] += variable_to_edge[eid * 2 + 1]->flow;
            diff[j] = rshift90(d, face_edgeOrients[fid][j]);
        }
        // check manifold is meaningless when it is still a singularity
//        if (diff[0] + diff[1] + diff[2] != Vector2i::Zero())
//            return true;
        for (int j = 0; j < 3; ++j) {
            int nj = (j + 2) % 3;
            if (-diff[j][0] * diff[nj][1] + diff[j][1] * diff[nj][0] < 0)
                return false;
        }
        return true;
    }
    void ApplyFlow(int v1, int v2, int flow) {
        for (auto& it : graph[v1]) {
            if (it.id == v2) {
                it.flow += flow;
                break;
            }
        }
    }
    int compute(std::vector<Vector2i>& edge_diff, std::vector<Vector3i>& face_edgeIds, std::vector<Vector3i>& face_edgeOrients, bool checkFlip) {
        int total_flow = 0;
        printf("%d %d\n", num, face_edgeIds.size());
        int count = 0;
        while (true) {
            count += 1;
            std::vector<int> vhash(num, 0);
            std::vector<std::pair<int, int> > q;
            q.push_back(std::make_pair(0, -1));
            vhash[0] = 1;
            int q_front = 0;
            bool found = false;
            while (q_front < q.size()) {
                int vert = q[q_front].first;
                int prev_loc = q[q_front].second;
                int prev_vert = -1;
                if (prev_loc != -1)
                    prev_vert = q[prev_loc].first;
                
                if (prev_vert != -1) {
                    ApplyFlow(prev_vert, vert, 1);
                    ApplyFlow(vert, prev_vert, -1);
                }
                for (auto& l : graph[vert]) {
                    if (vhash[l.id] || l.capacity <= l.flow)
                        continue;
                    l.flow += 1;
                    ApplyFlow(l.id, vert, -1);
                    if (vert != 0 && checkFlip && !CheckFlip(edge_diff, face_edgeIds, face_edgeOrients, (vert - 1) / 2)) {
                        l.flow -= 1;
                        ApplyFlow(l.id, vert, 1);
                        continue;
                    }
                    l.flow -= 1;
                    ApplyFlow(l.id, vert, 1);
                    q.push_back(std::make_pair(l.id, q_front));
                    vhash[l.id] = 1;
                    if (l.id == num - 1) {
                        found = true;
                        break;
                    }
                }
                if (prev_vert != -1) {
                    ApplyFlow(prev_vert, vert, -1);
                    ApplyFlow(vert, prev_vert, 1);
                }
                if (found)
                    break;
                q_front += 1;
            }
            if (q_front == q.size())
                break;
            int loc = q.size() - 1;
            /*
            std::vector<int> constraints_id;
            std::vector<int> locs;
            locs.push_back(loc);
            while (q[loc].second != -1) {
                constraints_id.push_back(q[loc].first);
                loc = q[loc].second;
                locs.push_back(loc);
            }
            constraints_id.push_back(q[loc].first);
             */
            while (q[loc].second != -1) {
                int current_v = q[loc].first;
                loc = q[loc].second;
                int prev_v = q[loc].first;
                ApplyFlow(prev_v, current_v, 1);
                ApplyFlow(current_v, prev_v, -1);
            }
            total_flow += 1;
        }
        return total_flow;
    }
    void Apply(std::unordered_map<int64_t, std::pair<int, int> >& edge_to_variable, std::vector<Vector2i>& edge_diff)
    {
        for (int i = 0; i < variable_to_edge.size(); ++i) {
            FlowInfo* info = variable_to_edge[i];
            if (info != 0)
                edge_diff[i/2][i%2] += info->flow;
        }
/*        for (int i = 0; i < graph.size(); ++i) {
            for (auto& l : graph[i]) {
                if (l.flow > 0) {
                    int64_t key = (int64_t)i * num + l.id;
                    auto q = edge_to_variable[key];
                    edge_diff[q.first / 2][q.first % 2] += q.second * l.flow;
                }
            }
        }*/
    }
    std::vector<std::list<FlowInfo> > graph;
};
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
