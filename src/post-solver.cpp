//
//  post-solver.cpp
//  parametrize
//
//  Created by Jingwei on 2/5/18.
//

#include "post-solver.hpp"

/*
 * TODO: Optimize O_quad, and possibly N_quad
 * Input:
 *  O_quad[i], initialized i-th vertex position of the quad mesh
 *  N_quad[i], initialized i-th vertex normal of the quad mesh
 *  Q_quad[i], initialized i-th vertex orientation of the quad mesh, guaranteed to be orthogonal to N_quad[i]
 *  F_quad[i], 4 vertex index of the i-th quad face
 *
 *  Concept: i-th directed edge is the (i%4)-th edge of the (i/4)-th face of the quad mesh
 *  V2E_quad[i], one directed edge from i-th vertex of the quad mesh
 *  E2E_quad[i], the reverse directed edge's index of the i-th directed edge of the quad mesh
 *
 *  V.col(i), i-th vertex position of the triangle mesh
 *  N.col(i), i-th vertex normal of the triangle mesh
 *  Q.col(i), i-th vertex orientation of the triangle mesh, guaranteed to be orthogonal to N.col(i)
 *  O.col(i), "quad position" associated with the i-th vertex in the triangle mesh (see InstantMesh position field)
 *  F.col(i), i-th triangle of the triangle mesh
 *
 *  V2E[i], one directed edge from the i-th vertex of the triangle mesh
 *  E2E[i], the reverse directed edge's index of the i-th directed edge of the triangle mesh
 *
 *  j = disajoint_tree.Index(i)
 *      the j-th vertex of the quad mesh is corresponding to the i-th vertex of the triangle mesh
 *      the relation is one-to-multiple
 *      O_quad can be viewed as an average of corresponding O
 *      N_quad can be viewed as an average of corresponding N
 *      Q_quad can be viewed as aggregation of corresponding Q
 *          Method that aggregates qi to qj with weights wi and wj:
 *              value = compat_orientation_extrinsic_4(qj, nj, qi, ni)
 *              result = (value.first * wj + value.second * wi).normalized()
 *
 * Output:
 *  Optimized O_quad, (possibly N_quad)
 */
void optimize_quad_positions(std::vector<Vector3d>& O_quad,
                             std::vector<Vector3d>& N_quad,
                             std::vector<Vector3d>& Q_quad,
                             std::vector<Vector4i>& F_quad,
                             VectorXi& V2E_quad,
                             VectorXi& E2E_quad,
                             MatrixXd& V,
                             MatrixXd& N,
                             MatrixXd& Q,
                             MatrixXd& O,
                             MatrixXi& F,
                             VectorXi& V2E,
                             VectorXi& E2E,
                             DisajointTree& disajoint_tree)
{
    //Information for the quad mesh
    printf("Quad mesh info:\n");
    printf("Number of vertices with normals and orientations: %d = %d = %d\n",
           O_quad.size(), N_quad.size(), Q_quad.size());
    printf("Number of faces: %d\n", F_quad.size());
    printf("Number of directed edges: %d\n", E2E_quad.size());
    //Information for the original mesh
    printf("Triangle mesh info:\n");
    printf("Number of vertices with normals, orientations and associated quad positions: %d = %d = %d = %d\n",
           V.cols(), N.cols(), Q.cols(), O.cols());
    printf("Number of faces: %d\n", F.cols());
    printf("Number of directed edges: %d\n", E2E.size());
    
    /* initial quad flips
     currently there are many flips, (82 flips in hand.obj)
     By uncommenting the ComputePosition() function call in Parametrizer.cpp, (current post linear solver)
     the flip number will be reduced to 14.
     As we discussed, we hope to implement this function to replace the previous solver
     */
    int flip_count = 0;
    for (int i = 0; i < F_quad.size(); ++i) {
        bool flipped = false;
        for (int j = 0; j < 4; ++j) {
            int v1 = F_quad[i][j];
            int v2 = F_quad[i][(j + 1) % 4];
            int v3 = F_quad[i][(j + 3) % 4];
            Vector3d face_norm = (O_quad[v2] - O_quad[v1]).cross(O_quad[v3] - O_quad[v1]);
            Vector3d vertex_norm = N_quad[v1];
            if (face_norm.dot(vertex_norm) < 0) {
                flipped = true;
            }
        }
        if (flipped) {
            flip_count++;
        }
    }
    printf("Flipped Quads: %d\n", flip_count);
}

