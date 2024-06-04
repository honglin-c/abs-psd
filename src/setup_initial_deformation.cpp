#include "setup_initial_deformation.h"

void setup_initial_deformation(const Eigen::MatrixXd V,
                               const Eigen::MatrixXi F,
                               const std::string pose_label,
                               const double deformation_magnitude,
                               const double deformation_ratio,
                               const double fixed_boundary_range,
                               Eigen::MatrixXd &U,
                               std::vector<unsigned int> &indices_fixed)
{
  U = V;
  indices_fixed.clear();

  // find the index of the leftmost and rightmost and topmost and bottommost and frontmost and backmost vertex
  int leftmost_idx = 0;
  int rightmost_idx = 0;
  int topmost_idx = 0;
  int bottommost_idx = 0;
  int frontmost_idx = 0;
  int backmost_idx = 0;
  for (int i = 0; i < V.rows(); ++i) {
    if (V(i, 0) < V(leftmost_idx, 0)) {
      leftmost_idx = i;
    }
    if (V(i, 0) > V(rightmost_idx, 0)) {
      rightmost_idx = i;
    }
    if (V(i, 1) < V(bottommost_idx, 1)) {
      bottommost_idx = i;
    }
    if (V(i, 1) > V(topmost_idx, 1)) {
      topmost_idx = i;
    }
    if (V(i, 2) < V(backmost_idx, 2)) {
      backmost_idx = i;
    }
    if (V(i, 2) > V(frontmost_idx, 2)) {
      frontmost_idx = i;
    }
  }

  double x_range = V(rightmost_idx, 0) - V(leftmost_idx, 0);
  double y_range = V(topmost_idx, 1) - V(bottommost_idx, 1);
  double z_range = V(frontmost_idx, 2) - V(backmost_idx, 2);

  // find all the vertices on the left and right and top and bottom and front and back
  std::vector<int> leftmost_vertices;
  std::vector<int> rightmost_vertices;
  std::vector<int> topmost_vertices;
  std::vector<int> bottommost_vertices;
  std::vector<int> frontmost_vertices;
  std::vector<int> backmost_vertices;
  std::vector<int> x_center_vertices;
  std::vector<int> y_center_vertices;
  std::vector<int> z_center_vertices;

  double x_mean = V(leftmost_idx, 0) + 0.5 * x_range;
  double y_mean = V(bottommost_idx, 1) + 0.5 * y_range;
  double z_mean = V(backmost_idx, 2) + 0.5 * z_range;
  for (int i = 0; i < V.rows(); ++i) {
    if (V(i, 0) < V(leftmost_idx, 0) + fixed_boundary_range) {
      leftmost_vertices.push_back(i);
    }
    if (V(i, 0) > V(rightmost_idx, 0) - fixed_boundary_range) {
      rightmost_vertices.push_back(i);
    }
    if (V(i, 1) < V(bottommost_idx, 1) + fixed_boundary_range) {
      bottommost_vertices.push_back(i);
    }
    if (V(i, 1) > V(topmost_idx, 1) - fixed_boundary_range) {
      topmost_vertices.push_back(i);
    }
    if (V(i, 2) < V(backmost_idx, 2) + fixed_boundary_range) {
      backmost_vertices.push_back(i);
    }
    if (V(i, 2) > V(frontmost_idx, 2) - fixed_boundary_range) {
      frontmost_vertices.push_back(i);
    }
    if (V(i, 0) > x_mean - 0.5*fixed_boundary_range && V(i, 0) < x_mean + 0.5*fixed_boundary_range) {
      x_center_vertices.push_back(i);
    }
    if (V(i, 1) > y_mean - 0.5*fixed_boundary_range && V(i, 1) < y_mean + 0.5*fixed_boundary_range) {
      y_center_vertices.push_back(i);
    }
    if (V(i, 2) > z_mean - 0.5*fixed_boundary_range && V(i, 2) < z_mean + 0.5*fixed_boundary_range) {
      z_center_vertices.push_back(i);
    }
  }

  if (pose_label == "stretch") {
    for (int i = 0; i < rightmost_vertices.size(); ++i) {
      U(rightmost_vertices[i], 0) += deformation_magnitude;
    }
  }
  else if (pose_label == "stretch_percentage") {
    for (int i = 0; i < rightmost_vertices.size(); ++i) {
      U(rightmost_vertices[i], 0) += deformation_ratio * x_range;
    }
  }
  else if (pose_label == "stretch_top") {
    for (int i = 0; i < topmost_vertices.size(); ++i) {
      U(topmost_vertices[i], 1) += deformation_magnitude;
    }
  }
  else if (pose_label == "stretch_front") {
    for (int i = 0; i < frontmost_vertices.size(); ++i) {
      U(frontmost_vertices[i], 2) += deformation_magnitude;
    }
  }
  else if (pose_label == "stretch_longest_axis") {
    // stretch the vertices on the longest axis by deformation_magnitude
    // find the longest axis
    if (x_range > y_range && x_range > z_range) {
      // x is the longest axis
      for (int i = 0; i < rightmost_vertices.size(); ++i) {
        U(rightmost_vertices[i], 0) += deformation_ratio * x_range;
      }
    }
    else if (y_range > x_range && y_range > z_range) {
      // y is the longest axis
      for (int i = 0; i < topmost_vertices.size(); ++i) {
        U(topmost_vertices[i], 1) += deformation_ratio * y_range;
      }
    }
    else {
      // z is the longest axis
      for (int i = 0; i < frontmost_vertices.size(); ++i) {
        U(frontmost_vertices[i], 2) += deformation_ratio * z_range;
      }
    }
  }
  else if (pose_label == "shear") {
    for (int i = 0; i < rightmost_vertices.size(); ++i) {
      U(rightmost_vertices[i], 1) += deformation_magnitude;
    }
  }
  else if (pose_label == "shear_percentage") {
    for (int i = 0; i < rightmost_vertices.size(); ++i) {
      U(rightmost_vertices[i], 1) += deformation_ratio * y_range;
    }
  }
  else if (pose_label == "stretch_shear") {
    for (int i = 0; i < rightmost_vertices.size(); ++i) {
      U(rightmost_vertices[i], 0) += deformation_magnitude;
      U(rightmost_vertices[i], 1) += deformation_magnitude;
    }
  }
  else if (pose_label == "stretch_shear_longest_axis") {
    // stretch the vertices on the longest axis by deformation_magnitude
    // find the longest axis
    if (x_range > y_range && x_range > z_range) {
      // x is the longest axis
      for (int i = 0; i < rightmost_vertices.size(); ++i) {
        U(rightmost_vertices[i], 0) += deformation_ratio * x_range;
        U(rightmost_vertices[i], 1) += deformation_ratio * y_range * 0.5;
      }
    }
    else if (y_range > x_range && y_range > z_range) {
      // y is the longest axis
      for (int i = 0; i < topmost_vertices.size(); ++i) {
        U(topmost_vertices[i], 1) += deformation_ratio * y_range;
        U(topmost_vertices[i], 0) += deformation_ratio * x_range * 0.5;
      }
    }
    else {
      // z is the longest axis
      for (int i = 0; i < frontmost_vertices.size(); ++i) {
        U(frontmost_vertices[i], 2) += deformation_ratio * z_range;
        U(frontmost_vertices[i], 0) += deformation_ratio * x_range * 0.5;
      }
    }
  }
  else if (pose_label == "stretch_shear_percentage") {
    for (int i = 0; i < rightmost_vertices.size(); ++i) {
      U(rightmost_vertices[i], 0) += deformation_ratio * x_range;
      U(rightmost_vertices[i], 1) += deformation_ratio * y_range;
    }
  }
  else if (pose_label == "compress") {
    for (int i = 0; i < rightmost_vertices.size(); ++i) {
      U(rightmost_vertices[i], 0) -= deformation_magnitude;
    }
  }
  else if (pose_label == "compress_percentage") {
    for (int i = 0; i < rightmost_vertices.size(); ++i) {
      U(rightmost_vertices[i], 0) -= deformation_ratio * x_range;
    }
  }
  else if (pose_label == "compress_top") {
    for (int i = 0; i < topmost_vertices.size(); ++i) {
      U(topmost_vertices[i], 1) -= deformation_magnitude;
    }
  }
  else if (pose_label == "compress_longest_axis") {
    // compress the vertices on the longest axis by deformation_magnitude
    // find the longest axis
    if (x_range > y_range && x_range > z_range) {
      // x is the longest axis
      for (int i = 0; i < rightmost_vertices.size(); ++i) {
        U(rightmost_vertices[i], 0) -= deformation_ratio * x_range;
      }
    }
    else if (y_range > x_range && y_range > z_range) {
      // y is the longest axis
      for (int i = 0; i < topmost_vertices.size(); ++i) {
        U(topmost_vertices[i], 1) -= deformation_ratio * y_range;
      }
    }
    else {
      // z is the longest axis
      for (int i = 0; i < frontmost_vertices.size(); ++i) {
        U(frontmost_vertices[i], 2) -= deformation_ratio * z_range;
      }
    }
  }
  else if (pose_label == "stretch_twist") {
    // rotate the rightmost vertices around the x axis by 90 degrees around the center of the rightmost vertices
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    for (int i = 0; i < rightmost_vertices.size(); ++i) {
      center += U.row(rightmost_vertices[i]);
    }
    center /= rightmost_vertices.size();
    for (int i = 0; i < rightmost_vertices.size(); ++i) {
      Eigen::Vector3d p = U.row(rightmost_vertices[i]);
      p -= center;
      p = Eigen::AngleAxisd(0.5 * M_PI, Eigen::Vector3d::UnitX()) * p;
      p += center;
      U.row(rightmost_vertices[i]) = p;
      U(rightmost_vertices[i], 0) += deformation_ratio * x_range;
    }
  }
  else if (pose_label == "stretch_twist_top") {
    // rotate the rightmost vertices around the x axis by 90 degrees around the center of the rightmost vertices
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    for (int i = 0; i < topmost_vertices.size(); ++i) {
      center += U.row(topmost_vertices[i]);
    }
    center /= topmost_vertices.size();
    for (int i = 0; i < topmost_vertices.size(); ++i) {
      Eigen::Vector3d p = U.row(topmost_vertices[i]);
      p -= center;
      p = Eigen::AngleAxisd(0.5 * M_PI, Eigen::Vector3d::UnitY()) * p;
      p += center;
      U.row(topmost_vertices[i]) = p;
      U(topmost_vertices[i], 1) += deformation_ratio * x_range;
    }
  }
  else if (pose_label == "teaser") {
    // move the topmost and bottommost vertices 
    for (int i = 0; i < topmost_vertices.size(); ++i) {
      U(topmost_vertices[i], 0) += deformation_magnitude;
      U(topmost_vertices[i], 2) += deformation_magnitude * 0.3;
    }
    for (int i = 0; i < bottommost_vertices.size(); ++i) {
      U(bottommost_vertices[i], 0) += deformation_magnitude;
      U(bottommost_vertices[i], 2) += deformation_magnitude * 0.3;
    }
  }

  if (pose_label == "compress_top" || pose_label == "stretch_top" || pose_label == "stretch_twist_top") {
    indices_fixed.insert(indices_fixed.end(), topmost_vertices.begin(), topmost_vertices.end());
    indices_fixed.insert(indices_fixed.end(), bottommost_vertices.begin(), bottommost_vertices.end());
  }
  else if (pose_label == "stretch_front") {
    indices_fixed.insert(indices_fixed.end(), backmost_vertices.begin(), backmost_vertices.end());
    indices_fixed.insert(indices_fixed.end(), frontmost_vertices.begin(), frontmost_vertices.end());
  }
  else if (pose_label == "compress_longest_axis" || pose_label == "stretch_longest_axis" || pose_label == "stretch_shear_longest_axis") {
    // fix the two ends of the longest axis
    if (x_range > y_range && x_range > z_range) {
      // x is the longest axis
      indices_fixed.insert(indices_fixed.end(), rightmost_vertices.begin(), rightmost_vertices.end());
      indices_fixed.insert(indices_fixed.end(), leftmost_vertices.begin(), leftmost_vertices.end());
    }
    else if (y_range > x_range && y_range > z_range) {
      // y is the longest axis
      indices_fixed.insert(indices_fixed.end(), topmost_vertices.begin(), topmost_vertices.end());
      indices_fixed.insert(indices_fixed.end(), bottommost_vertices.begin(), bottommost_vertices.end());
    }
    else {
      // z is the longest axis
      indices_fixed.insert(indices_fixed.end(), backmost_vertices.begin(), backmost_vertices.end());
      indices_fixed.insert(indices_fixed.end(), frontmost_vertices.begin(), frontmost_vertices.end());
    }
  }
  else if (pose_label.substr(0,6) == "teaser") {
    indices_fixed.insert(indices_fixed.end(), topmost_vertices.begin(), topmost_vertices.end());
    indices_fixed.insert(indices_fixed.end(), bottommost_vertices.begin(), bottommost_vertices.end());
    indices_fixed.insert(indices_fixed.end(), backmost_vertices.begin(), backmost_vertices.end());
    indices_fixed.insert(indices_fixed.end(), frontmost_vertices.begin(), frontmost_vertices.end());
  }
  else {
    indices_fixed.insert(indices_fixed.end(), leftmost_vertices.begin(), leftmost_vertices.end());
    indices_fixed.insert(indices_fixed.end(), rightmost_vertices.begin(), rightmost_vertices.end());
  }
}