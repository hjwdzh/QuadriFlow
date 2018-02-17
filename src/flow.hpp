#ifndef FLOW_H_
#define FLOW_H_
#include <Eigen/Core>
#include <list>
#include <map>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>

using namespace boost;
using namespace Eigen;

typedef int EdgeWeightType;

// clang-format off
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
// clang-format on

void AddEdge(Traits::vertex_descriptor& v1, Traits::vertex_descriptor& v2,
             property_map<Graph, edge_reverse_t>::type& rev, const double capacity, Graph& g);

void AddDirectEdge(Traits::vertex_descriptor& v1, Traits::vertex_descriptor& v2,
                   property_map<Graph, edge_reverse_t>::type& rev, const int capacity,
                   const int inv_capacity, Graph& g, Traits::edge_descriptor& e1,
                   Traits::edge_descriptor& e2);

class MaxFlowHelper {
   public:
    MaxFlowHelper() {}
    virtual ~MaxFlowHelper(){};
    virtual void resize(int n, int m) = 0;
    virtual int compute() = 0;
    virtual void AddEdge(int x, int y, int c, int rc, int v) = 0;
    virtual void Apply(std::vector<Vector2i>& edge_diff) = 0;
};

class BoykovMaxFlowHelper : public MaxFlowHelper {
   public:
    BoykovMaxFlowHelper() {
        rev = get(edge_reverse, g);
        num = 0;
    }
    int num;
    void resize(int n, int m) {
        vertex_descriptors.resize(n);
        for (int i = 0; i < n; ++i) {
            vertex_descriptors[i] = add_vertex(g);
        }
        num = n;
    }
    int compute() {
        EdgeWeightType flow =
            boykov_kolmogorov_max_flow(g, vertex_descriptors.front(), vertex_descriptors.back());
        return flow;
    }
    void AddEdge(int x, int y, int c, int rc, int v) {
        Traits::edge_descriptor e1, e2;
        AddDirectEdge(vertex_descriptors[x], vertex_descriptors[y], rev, c, rc, g, e1, e2);
        if (v != -1) {
            edge_to_variables[e1] = std::make_pair(v, -1);
            edge_to_variables[e2] = std::make_pair(v, 1);
        }
    }
    void Apply(std::vector<Vector2i>& edge_diff) {
        property_map<Graph, edge_capacity_t>::type capacity = get(edge_capacity, g);
        property_map<Graph, edge_residual_capacity_t>::type residual_capacity =
            get(edge_residual_capacity, g);

        graph_traits<Graph>::vertex_iterator u_iter, u_end;
        graph_traits<Graph>::out_edge_iterator ei, e_end;
        for (tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter)
            for (tie(ei, e_end) = out_edges(*u_iter, g); ei != e_end; ++ei)
                if (capacity[*ei] > 0) {
                    int flow = (capacity[*ei] - residual_capacity[*ei]);
                    if (flow > 0) {
                        auto it = edge_to_variables.find(*ei);
                        if (it != edge_to_variables.end()) {
                            edge_diff[it->second.first / 2][it->second.first % 2] +=
                                it->second.second * flow;
                        }
                        /*
                        int64_t key = (int64_t)*u_iter * num + target(*ei, g);
                        auto q = edge_to_variable[key];
                        edge_diff[q.first / 2][q.first % 2] += q.second * flow;
                        if (abs(edge_diff[q.first/2][q.first%2]) > 2) {
                            printf("Edge %d %d\n", *u_iter, target(*ei, g));
                            printf("Apply error 1: %d %d %d\n",
                                   edge_diff[q.first / 2][q.first % 2],
                                   q.second, flow);
                            printf("capacity %d\n", capacity[*ei]);
                            exit(0);
                        }
                         */
                    }
                }
    }
    Graph g;
    property_map<Graph, edge_reverse_t>::type rev;
    std::vector<Traits::vertex_descriptor> vertex_descriptors;
    std::map<Traits::edge_descriptor, std::pair<int, int>> edge_to_variables;
};

int flow(std::vector<std::map<int, std::pair<int, int>>>& graph);

inline void AddEdge(Traits::vertex_descriptor& v1, Traits::vertex_descriptor& v2,
                    property_map<Graph, edge_reverse_t>::type& rev, const double capacity,
                    Graph& g) {
    Traits::edge_descriptor e1 = add_edge(v1, v2, g).first;
    Traits::edge_descriptor e2 = add_edge(v2, v1, g).first;
    put(edge_capacity, g, e1, capacity);
    put(edge_capacity, g, e2, capacity);

    rev[e1] = e2;
    rev[e2] = e1;
}

inline void AddDirectEdge(Traits::vertex_descriptor& v1, Traits::vertex_descriptor& v2,
                          property_map<Graph, edge_reverse_t>::type& rev, const int capacity,
                          const int inv_capacity, Graph& g, Traits::edge_descriptor& e1,
                          Traits::edge_descriptor& e2) {
    e1 = add_edge(v1, v2, g).first;
    e2 = add_edge(v2, v1, g).first;
    put(edge_capacity, g, e1, capacity);
    put(edge_capacity, g, e2, inv_capacity);

    rev[e1] = e2;
    rev[e2] = e1;
}

class ECMaxFlowHelper : public MaxFlowHelper {
   public:
    struct FlowInfo {
        int id;
        int capacity, flow;
        int v, d;
        FlowInfo* rev;
    };
    struct SearchInfo {
        SearchInfo(int _id, int _prev_id, FlowInfo* _info)
            : id(_id), prev_id(_prev_id), info(_info) {}
        int id;
        int prev_id;
        FlowInfo* info;
    };
    ECMaxFlowHelper() { num = 0; }
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
        flow.v = v;
        flow.d = -1;
        graph[x].push_back(flow);
        auto& f1 = graph[x].back();
        flow.id = x;
        flow.capacity = rc;
        flow.flow = 0;
        flow.v = v;
        flow.d = 1;
        graph[y].push_back(flow);
        auto& f2 = graph[y].back();
        f2.rev = &f1;
        f1.rev = &f2;
    }

    void ApplyFlow(int v1, int v2, int flow) {
        for (auto& it : graph[v1]) {
            if (it.id == v2) {
                it.flow += flow;
                break;
            }
        }
    }
    int compute() {
        int total_flow = 0;
        int count = 0;
        while (true) {
            count += 1;
            std::vector<int> vhash(num, 0);
            std::vector<SearchInfo> q;
            q.push_back(SearchInfo(0, -1, 0));
            vhash[0] = 1;
            int q_front = 0;
            bool found = false;
            while (q_front < q.size()) {
                int vert = q[q_front].id;
                for (auto& l : graph[vert]) {
                    if (vhash[l.id] || l.capacity <= l.flow) continue;
                    q.push_back(SearchInfo(l.id, q_front, &l));
                    vhash[l.id] = 1;
                    if (l.id == num - 1) {
                        found = true;
                        break;
                    }
                }
                if (found) break;
                q_front += 1;
            }
            if (q_front == q.size()) break;
            int loc = q.size() - 1;
            while (q[loc].prev_id != -1) {
                q[loc].info->flow += 1;
                q[loc].info->rev->flow -= 1;
                loc = q[loc].prev_id;
                //                int prev_v = q[loc].id;
                //                ApplyFlow(prev_v, current_v, 1);
                //                ApplyFlow(current_v, prev_v, -1);
            }
            total_flow += 1;
        }
        return total_flow;
    }
    void Apply(std::vector<Vector2i>& edge_diff) {
        for (int i = 0; i < graph.size(); ++i) {
            for (auto& flow : graph[i]) {
                if (flow.flow > 0 && flow.v != -1) {
                    if (flow.flow > 0) {
                        edge_diff[flow.v / 2][flow.v % 2] += flow.d * flow.flow;
                        if (abs(edge_diff[flow.v / 2][flow.v % 2]) > 2) {
                        }
                    }
                }
            }
        }
    }
    std::vector<std::list<FlowInfo>> graph;
};
#endif
