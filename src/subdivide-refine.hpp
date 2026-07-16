#pragma once

#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <unordered_map>
#include <Eigen/Core>
#include "config.hpp"

#ifdef WITH_OMP
#include <omp.h>
#endif

namespace qflow {

using namespace Eigen;

// ---------------------------------------------------------------------------
// Closest point on a triangle to query point p (Eberly's method)
// ---------------------------------------------------------------------------
inline Vector3d ClosestPointOnTriangle(const Vector3d& p,
                                        const Vector3d& a,
                                        const Vector3d& b,
                                        const Vector3d& c) {
    Vector3d ab = b - a, ac = c - a, ap = p - a;
    double d1 = ab.dot(ap), d2 = ac.dot(ap);
    if (d1 <= 0.0 && d2 <= 0.0) return a;

    Vector3d bp = p - b;
    double d3 = ab.dot(bp), d4 = ac.dot(bp);
    if (d3 >= 0.0 && d4 <= d3) return b;

    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
        return a + (d1 / (d1 - d3)) * ab;

    Vector3d cp = p - c;
    double d5 = ab.dot(cp), d6 = ac.dot(cp);
    if (d6 >= 0.0 && d5 <= d6) return c;

    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
        return a + (d2 / (d2 - d6)) * ac;

    double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + w * (c - b);
    }
    double inv = 1.0 / (va + vb + vc);
    return a + (vb * inv) * ab + (vc * inv) * ac;
}

// ---------------------------------------------------------------------------
// Grid-based spatial index for fast closest-point projection onto a mesh.
// V columns are 3D vertices, F columns are triangle vertex indices.
// ---------------------------------------------------------------------------
struct MeshProjector {
    const MatrixXd& V;
    const MatrixXi& F;
    Vector3d bmin;
    double cell;
    int gx, gy, gz;
    std::vector<std::vector<int>> grid;

    MeshProjector(const MatrixXd& V_, const MatrixXi& F_) : V(V_), F(F_) {
        bmin = V.rowwise().minCoeff();
        Vector3d bmax = V.rowwise().maxCoeff();
        Vector3d extent = bmax - bmin + Vector3d(1e-6, 1e-6, 1e-6);

        // Target ~8 triangles per cell
        double vol = extent.prod();
        cell = std::cbrt(vol / std::max(1, (int)F.cols()) * 8.0);
        if (cell < 1e-10) cell = 1.0;

        gx = std::max(1, (int)std::ceil(extent[0] / cell));
        gy = std::max(1, (int)std::ceil(extent[1] / cell));
        gz = std::max(1, (int)std::ceil(extent[2] / cell));

        // Cap to avoid excessive memory (~64 MB)
        while ((long long)gx * gy * gz > 8000000LL) {
            gx = std::max(1, gx - 1);
            gy = std::max(1, gy - 1);
            gz = std::max(1, gz - 1);
            cell = std::max({extent[0] / gx, extent[1] / gy, extent[2] / gz});
        }

        grid.resize((long long)gx * gy * gz);

        // Insert each triangle into every cell its bounding box overlaps
        for (int i = 0; i < F.cols(); ++i) {
            Vector3d a = V.col(F(0, i)), b = V.col(F(1, i)), c = V.col(F(2, i));
            Vector3d tmin = (a.cwiseMin(b).cwiseMin(c) - bmin).cwiseMax(Vector3d::Zero());
            Vector3d tmax =  a.cwiseMax(b).cwiseMax(c) - bmin;
            int x0 = (int)(tmin[0] / cell), x1 = std::min(gx-1, (int)(tmax[0] / cell));
            int y0 = (int)(tmin[1] / cell), y1 = std::min(gy-1, (int)(tmax[1] / cell));
            int z0 = (int)(tmin[2] / cell), z1 = std::min(gz-1, (int)(tmax[2] / cell));
            for (int x = x0; x <= x1; ++x)
                for (int y = y0; y <= y1; ++y)
                    for (int z = z0; z <= z1; ++z)
                        grid[((long long)x * gy + y) * gz + z].push_back(i);
        }
    }

    // Return the closest point on the mesh to p.
    Vector3d project(const Vector3d& p) const {
        Vector3d lp = p - bmin;
        int cx = std::min(gx-1, std::max(0, (int)(lp[0] / cell)));
        int cy = std::min(gy-1, std::max(0, (int)(lp[1] / cell)));
        int cz = std::min(gz-1, std::max(0, (int)(lp[2] / cell)));

        double best_d2 = 1e30;
        Vector3d best_p = p;

        int max_r = std::max({cx, gx-1-cx, cy, gy-1-cy, cz, gz-1-cz});
        for (int r = 0; r <= max_r; ++r) {
            // Once we've found a hit, the closest possible triangle in the
            // next shell is at least (r-1)*cell away; early-exit if that's
            // already farther than our best.
            if (r > 1 && best_d2 <= ((r - 1) * cell) * ((r - 1) * cell)) break;

            for (int dx = -r; dx <= r; ++dx)
            for (int dy = -r; dy <= r; ++dy)
            for (int dz = -r; dz <= r; ++dz) {
                if (std::abs(dx) != r && std::abs(dy) != r && std::abs(dz) != r) continue;
                int nx = cx+dx, ny = cy+dy, nz = cz+dz;
                if (nx < 0 || nx >= gx || ny < 0 || ny >= gy || nz < 0 || nz >= gz) continue;
                for (int fi : grid[((long long)nx * gy + ny) * gz + nz]) {
                    Vector3d cp = ClosestPointOnTriangle(p,
                        V.col(F(0, fi)), V.col(F(1, fi)), V.col(F(2, fi)));
                    double d2 = (cp - p).squaredNorm();
                    if (d2 < best_d2) { best_d2 = d2; best_p = cp; }
                }
            }
        }
        return best_p;
    }
};

// ---------------------------------------------------------------------------
// One level of simple midpoint subdivision on a pure quad mesh.
// Each quad (v0,v1,v2,v3) → 4 quads via edge midpoints + face centre.
// Edge midpoints are shared between adjacent quads.
// ---------------------------------------------------------------------------
inline void SubdivideQuadMesh(std::vector<Vector3d>& verts,
                               std::vector<Vector4i>& faces) {
    std::unordered_map<int64_t, int> edge_map;
    edge_map.reserve(faces.size() * 5);

    auto midpoint = [&](int a, int b) -> int {
        if (a > b) std::swap(a, b);
        int64_t key = ((int64_t)a << 32) | (uint32_t)b;
        auto it = edge_map.find(key);
        if (it != edge_map.end()) return it->second;
        int idx = (int)verts.size();
        verts.push_back((verts[a] + verts[b]) * 0.5);
        edge_map[key] = idx;
        return idx;
    };

    int nf = (int)faces.size();
    std::vector<Vector4i> new_faces;
    new_faces.reserve(nf * 4);

    for (int i = 0; i < nf; ++i) {
        int v0 = faces[i][0], v1 = faces[i][1],
            v2 = faces[i][2], v3 = faces[i][3];
        // Face centre
        int fc = (int)verts.size();
        verts.push_back((verts[v0] + verts[v1] + verts[v2] + verts[v3]) * 0.25);
        // Edge midpoints (shared with adjacent faces)
        int e01 = midpoint(v0, v1), e12 = midpoint(v1, v2),
            e23 = midpoint(v2, v3), e30 = midpoint(v3, v0);
        // Four child quads, winding consistent with the parent
        new_faces.push_back(Vector4i(v0, e01, fc, e30));
        new_faces.push_back(Vector4i(e01, v1, e12, fc));
        new_faces.push_back(Vector4i(fc, e12, v2, e23));
        new_faces.push_back(Vector4i(e30, fc, e23, v3));
    }
    faces = std::move(new_faces);
}

// ---------------------------------------------------------------------------
// Subdivide a quad mesh `levels` times and project every vertex back onto
// the original triangle mesh (origV, origF), which must be in the same
// coordinate space as verts.
// ---------------------------------------------------------------------------
inline void SubdivideAndRefit(std::vector<Vector3d>& verts,
                               std::vector<Vector4i>& faces,
                               const MatrixXd& origV,
                               const MatrixXi& origF,
                               int levels) {
    if (levels <= 0) return;

    // Subdivide
    for (int l = 0; l < levels; ++l) {
        printf("  Subdivide level %d/%d: %d faces", l + 1, levels, (int)faces.size());
        fflush(stdout);
        SubdivideQuadMesh(verts, faces);
        printf(" -> %d faces, %d vertices\n", (int)faces.size(), (int)verts.size());
        fflush(stdout);
    }

    // Build spatial index once, project all vertices in parallel
    printf("  Building spatial index (%d triangles)...\n", (int)origF.cols());
    fflush(stdout);
    MeshProjector proj(origV, origF);

    printf("  Projecting %d vertices onto original mesh...\n", (int)verts.size());
    fflush(stdout);
#ifdef WITH_OMP
#pragma omp parallel for schedule(dynamic, 512)
#endif
    for (int i = 0; i < (int)verts.size(); ++i)
        verts[i] = proj.project(verts[i]);

    printf("  Refit complete.\n");
    fflush(stdout);
}

} // namespace qflow
