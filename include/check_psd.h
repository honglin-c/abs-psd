#ifndef CHECK_PSD_H
#define CHECK_PSD_H

#include <Eigen/SparseCholesky>

bool check_psd(const Eigen::SparseMatrix<double> &H)
{ 
   	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt;
    llt.compute(H);
    if (llt.info() != Eigen::Success)
    {
      return false;
    }
    return true;
}


#endif