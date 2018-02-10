#include "parametrizer.hpp"
#include <unordered_map>
#include <vector>

double Parametrizer::QuadEnergy(std::vector<int>& loop_vertices, std::vector<Vector4i>& res_quads, int level)
{
    if (loop_vertices.size() == 4) {
        double energy = 0;
        for (int j = 0; j < 4; ++j) {
            int v0 = loop_vertices[j];
            int v2 = loop_vertices[(j + 1) % 4];
            int v1 = loop_vertices[(j + 3) % 4];
            Vector3d pt1 = (O_compact[v1] - O_compact[v0]).normalized();
            Vector3d pt2 = (O_compact[v2] - O_compact[v0]).normalized();
            Vector3d n = pt1.cross(pt2);
            double sina = n.norm();
            if (n.dot(N_compact[v0]) < 0)
                sina = -sina;
            double cosa = pt1.dot(pt2);
            double angle = atan2(sina, cosa) / 3.141592654 * 180.0;
            if (angle < 0)
                angle = 360 + angle;
            energy += angle * angle;
        }
        res_quads.push_back(Vector4i(loop_vertices[0],loop_vertices[3],loop_vertices[2],loop_vertices[1]));
        return energy;
    }
    double max_energy = 1e30;
    for (int seg1 = 2; seg1 < loop_vertices.size(); seg1 += 2) {
        for (int seg2 = seg1 + 1; seg2 < loop_vertices.size(); seg2 += 2) {
            std::vector<Vector4i> quads[4];
            std::vector<int> vertices = {loop_vertices[0], loop_vertices[1], loop_vertices[seg1], loop_vertices[seg2]};
            double energy = 0;
            energy += QuadEnergy(vertices, quads[0], level + 1);
            if (seg1 > 2) {
                std::vector<int> vertices(loop_vertices.begin() + 1, loop_vertices.begin() + seg1);
                vertices.push_back(loop_vertices[seg1]);
                energy += QuadEnergy(vertices, quads[1], level + 1);
            }
            if (seg2 != seg1 + 1) {
                std::vector<int> vertices(loop_vertices.begin() + seg1, loop_vertices.begin() + seg2);
                vertices.push_back(loop_vertices[seg2]);
                energy += QuadEnergy(vertices, quads[2], level + 2);
            }
            if (seg2 + 1 != loop_vertices.size()) {
                std::vector<int> vertices(loop_vertices.begin() + seg2, loop_vertices.end());
                vertices.push_back(loop_vertices[0]);
                energy += QuadEnergy(vertices, quads[3], level + 1);
            }
            if (max_energy > energy) {
                max_energy = energy;
                res_quads.clear();
                for (int i = 0; i < 4; ++i) {
                    for (auto& v : quads[i]) {
                        res_quads.push_back(v);
                    }
                }
            }
        }
    }
    return max_energy;
}

void Parametrizer::FixHoles() {
    std::vector<int> detected_boundary(E2E_compact.size(), 0);
    for (int i = 0; i < E2E_compact.size(); ++i) {
        if (detected_boundary[i] != 0 || E2E_compact[i] != -1)
            continue;
        std::vector<int> loop_edges;
        int current_e = i;
        
        while (detected_boundary[current_e] == 0) {
            detected_boundary[current_e] = 1;
            loop_edges.push_back(current_e);
            current_e = current_e / 4 * 4 + (current_e + 1) % 4;
            while (E2E_compact[current_e] != -1) {
                current_e = E2E_compact[current_e];
                current_e = current_e / 4 * 4 + (current_e + 1) % 4;
            }
        }
        std::vector<int> loop_vertices(loop_edges.size());
        for (int j = 0; j < loop_edges.size(); ++j) {
            loop_vertices[j] = F_compact[loop_edges[j]/4][loop_edges[j]%4];
        }
        std::vector<std::vector<int> > loop_vertices_array, loop_edges_array;
        std::unordered_map<int, int> map_loops;
        for (int i = 0; i < loop_vertices.size(); ++i) {
            if (map_loops.count(loop_vertices[i])) {
                int j = map_loops[loop_vertices[i]];
                loop_vertices_array.push_back(std::vector<int>());
                loop_edges_array.push_back(std::vector<int>());
                if (i - j > 3 && (i - j) % 2 == 0) {
                    for (int k = j; k < i; ++k) {
                        if (map_loops.count(loop_vertices[k])) {
                            loop_vertices_array.back().push_back(loop_vertices[k]);
                            loop_edges_array.back().push_back(loop_edges[k]);
                            map_loops.erase(loop_vertices[k]);
                        }
                    }
                }
            }
            map_loops[loop_vertices[i]] = i;
        }
        if (map_loops.size() > 2 && map_loops.size() % 2 == 0) {
            loop_vertices_array.push_back(std::vector<int>());
            loop_edges_array.push_back(std::vector<int>());
            for (int k = 0; k < loop_vertices.size(); ++k) {
                if (map_loops.count(loop_vertices[k])) {
                    if (map_loops.count(loop_vertices[k])) {
                        loop_vertices_array.back().push_back(loop_vertices[k]);
                        loop_edges_array.back().push_back(loop_edges[k]);
                        map_loops.erase(loop_vertices[k]);
                    }
                }
            }
        }
        for (int i = 0; i < loop_vertices_array.size(); ++i) {
            auto& loop_vertices = loop_vertices_array[i];
            std::vector<Vector4i> quads;
#ifdef LOG_OUTPUT
            printf("Compute energy for loop: %d\n", (int)loop_vertices.size());
#endif
            QuadEnergy(loop_vertices, quads, 0);
            for (auto& p : quads) {
                F_compact.push_back(p);
            }
        }
    }
}

void Parametrizer::FixFlipHierarchy() {
    Hierarchy fh;
    fh.DownsampleEdgeGraph(face_edgeOrients, face_edgeIds, edge_diff, -1);
    fh.FixFlip();
    fh.UpdateGraphValue(face_edgeOrients, face_edgeIds, edge_diff);
}


