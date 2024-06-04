#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrix<double> &P, unsigned int q_size, unsigned int dim, const std::vector<unsigned int> indices) {
    P = Eigen::SparseMatrix<double>(q_size - dim * indices.size(), q_size);
    // calculate the index for non-fixed vertices
    Eigen::VectorXi nonfixed_id(q_size/dim);
    nonfixed_id.setZero();
    for(int i = 0; i < indices.size(); i++)
        nonfixed_id[indices[i]] = -1;
    int count = 0;
    for(int i = 0; i < q_size/dim; i++) {
        if(nonfixed_id[i] != -1)
            nonfixed_id[i] = count++;
    }

    // set the cooresponding dimxdim blocks to identity
    std::vector<Eigen::Triplet<double>> triplets;
    auto setIdentity = [](std::vector<Eigen::Triplet<double>>& triplets, int dim, int start_row, int start_col) {
        for(int i = 0; i < dim; i++)
            triplets.push_back(Eigen::Triplet<double>(start_row+i, start_col+i, 1));
    };
    for(int i = 0; i < q_size/dim; i++) {
        if(nonfixed_id[i] != -1) {
            int start_row = nonfixed_id[i] * dim;
            int start_col = i * dim;
            setIdentity(triplets, dim, start_row, start_col);
        }
    }
    P.setFromTriplets(triplets.begin(), triplets.end());  
}