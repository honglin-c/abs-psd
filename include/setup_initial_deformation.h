#ifndef SETUP_INITIAL_DEFORMATION_H
#define SETUP_INITIAL_DEFORMATION_H

#include <Eigen/Dense>
#include <vector>
#include <string>

void setup_initial_deformation(const Eigen::MatrixXd V,
                               const Eigen::MatrixXi F,
                               const std::string pose_label,
                               const double deformation_magnitude,
                               const double deformation_ratio,
                               const double fixed_boundary_range,
                               Eigen::MatrixXd &U,
                               std::vector<unsigned int> &indices_fixed);
                          

#endif