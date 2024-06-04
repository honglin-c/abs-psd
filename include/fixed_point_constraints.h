#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>

//Input:
//  q_size - total number of scalar generalized coordinates (3 times number of vertices in the mesh)
//  indices - indices (row ids in V) for fixed vertices
//  dim  - dimension of V (3 for 3D mesh)
//Output:
//  P - 3*mx3*n sparse matrix which projects out fixed vertices
void fixed_point_constraints(Eigen::SparseMatrix<double> &P, unsigned int q_size, unsigned int dim, std::vector<unsigned int> indices);