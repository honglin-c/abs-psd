#include <igl/readMESH.h>
#include <igl/writeOBJ.h>
#include <igl/writeMESH.h>

#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>

#include <igl/boundary_facets.h>
#include <igl/per_vertex_normals.h>
#include <igl/fit_rotations.h>
#include <igl/volume.h>

#include <igl/opengl/glfw/Viewer.h>
#include <thread>
#include <mutex>
#include <string>
#include <filesystem>

#include "cxxopts.hpp"
#include "fixed_point_constraints.h"
#include "setup_initial_deformation.h"
 
#include <chrono>
using namespace std::chrono;


int main(int argc, char** argv)
{
  // parse command line arguments
  cxxopts::Options options("example", "This program performs quasistatic simulation using Projected Newton with absolute eigenvalue projection by defaults");

  options.add_options()
    ("clamp", "use eigenvalue clamping strategy instread", cxxopts::value<bool>()->default_value("false"))
    ("p,epsilon", "Projection threshold for the eigenvalue clamping", cxxopts::value<std::string>()->default_value("0"))
    ("n,mesh_name", "Mesh name", cxxopts::value<std::string>()->default_value("bimba"))
    ("l,pose_label", "Pose label", cxxopts::value<std::string>()->default_value("stretch"))
    ("r,lambda_mu_ratio", "Ratio between lambda and mu (lambda / mu)", cxxopts::value<double>()->default_value("1e2"))
    ("g,deformation_magnitude", "The magnitude of the deformation", cxxopts::value<double>()->default_value("2.0"))
    ("t,deformation_ratio", "The ratio of the deformation", cxxopts::value<double>()->default_value("2.0"))
    ("b,fixed_boundary_range", "The range of fixed vertices on the boundary", cxxopts::value<double>()->default_value("0.1"))
    ("c,convergence_eps", "The convergence threshold", cxxopts::value<double>()->default_value("1e-5"))
    ("ym", "Young's modulus (need to set both YM and PR to enable this option, otherwise lambda_mu_ratio is used instead)", cxxopts::value<double>()->default_value("0"))
    ("pr", "Poisson's ratio (need to set both YM and PR to enable this option, otherwise lambda_mu_ratio is used instead)", cxxopts::value<double>()->default_value("0"))
    ("h,help", "show help")
  ;
  
  auto result = options.parse(argc, argv);
  if (result.count("help"))
  {
      std::cout << options.help() << "\n";
      return 0;
  }

  const bool clamp = result["clamp"].as<bool>();
  std::string eps_str = result["epsilon"].as<std::string>();
  double eps = std::stod(eps_str);
  std::string mesh_name = result["mesh_name"].as<std::string>();
  std::string pose_label = result["pose_label"].as<std::string>();
  const double lambda_mu_ratio = result["lambda_mu_ratio"].as<double>();
  const double deformation_magnitude = result["deformation_magnitude"].as<double>();
  const double deformation_ratio = result["deformation_ratio"].as<double>();
  const double fixed_boundary_range = result["fixed_boundary_range"].as<double>();
  const double convergence_eps = result["convergence_eps"].as<double>();
  const double YM = result["ym"].as<double>();
  const double PR = result["pr"].as<double>();

  if (!clamp) {
    // we use eps < 0 as a flag for abs projection, see lines 71-78 in our modified `TinyAD/include/TinyAD/Utils/HessianProjection.hh`
    eps_str = "abs";
    eps = -1;  
  }
  else {
    eps_str = "clamp";
  }
    
  // create output folder 
  if (!std::filesystem::exists("../results/"))
    std::filesystem::create_directory("../results/");

  const std::string output_folder = "../results/" + mesh_name + "_" + pose_label + "/";
  {
    if (!std::filesystem::exists(output_folder))
      std::filesystem::create_directory(output_folder);
    if (!std::filesystem::exists(output_folder + "obj/"))
      std::filesystem::create_directory(output_folder + "obj/");
    if (!std::filesystem::exists(output_folder + "obj/" + mesh_name + "_" + eps_str + "/"))
      std::filesystem::create_directory(output_folder + "obj/" + mesh_name + "_" + eps_str + "/");
    if (!std::filesystem::exists(output_folder + "hist/"))
      std::filesystem::create_directory(output_folder + "hist/");
    if (!std::filesystem::exists(output_folder + "iter/"))
      std::filesystem::create_directory(output_folder + "iter/");
  }

  const std::string output_tag = "mesh_" + mesh_name + "_" + eps_str;

  // set up viewer
  igl::opengl::glfw::Viewer viewer;

  const double MU = (YM > 0 && PR > 0) ? YM / (2 * (1 + PR)) : 3.3344e7;
  double LAMBDA = (YM > 0 && PR > 0) ? YM * PR / ((1 + PR) * (1 - 2 * PR)) : MU * lambda_mu_ratio;
 
  {
    TINYAD_DEBUG_OUT("Clamp: " << clamp);
    TINYAD_DEBUG_OUT("MU: " << MU);
    TINYAD_DEBUG_OUT("LAMBDA: " << LAMBDA);
    TINYAD_DEBUG_OUT("Projection threshold: " << eps);
    TINYAD_DEBUG_OUT("Mesh name: " << mesh_name);
    TINYAD_DEBUG_OUT("Pose label: " << pose_label);
    TINYAD_DEBUG_OUT("Lambda / mu ratio: " << lambda_mu_ratio);
    TINYAD_DEBUG_OUT("Deformation magnitude: " << deformation_magnitude);
    TINYAD_DEBUG_OUT("Deformation ratio: " << deformation_ratio);
    TINYAD_DEBUG_OUT("Fixed vertices boundary range: " << fixed_boundary_range);
    TINYAD_DEBUG_OUT("Convergence threshold: " << convergence_eps);
    TINYAD_DEBUG_OUT("Young's modulus: " << YM);
    TINYAD_DEBUG_OUT("Poisson's ratio: " << PR);
  }

  Eigen::MatrixXd V, U; // #V-by-3 3D vertex positions
  Eigen::MatrixXi F, FF; // #T-by-4 indices into V
  Eigen::VectorXi TriTag, TetTag;
  if (std::filesystem::exists(std::string(SOURCE_PATH) + "/data/" + mesh_name + ".mesh"))
    igl::readMESH(std::string(SOURCE_PATH) + "/data/" + mesh_name + ".mesh", V, F, FF);
  else if (std::filesystem::exists(std::string(SOURCE_PATH) + "/data/" + mesh_name + ".msh"))
    igl::readMSH(std::string(SOURCE_PATH) + "/data/" + mesh_name + ".msh", V, FF, F, TriTag, TetTag);
  else {
    std::cout << "Mesh " << mesh_name << " not found!" << std::endl;
    exit(1);
  }

  TINYAD_DEBUG_OUT("Read mesh with " << V.rows() << " vertices and " << F.rows() << " tetrahedrons.");

  // get boundary vertices
  igl::boundary_facets(F, FF);
  FF = FF.rowwise().reverse().eval();

  TINYAD_DEBUG_OUT("Boundary has " << FF.rows() << " faces.");

  // normalize the mesh V
  Eigen::RowVector3d mean = V.colwise().mean();
  V.rowwise() -= mean;
  double max_norm = V.rowwise().norm().maxCoeff();
  V /= max_norm;

  // set up mesh U
  U = V;

  // fixed point constraints
  Eigen::SparseMatrix<double> P;
  std::vector<unsigned int> indices_fixed;

  // set up fixed point constraints and initial deformation
  setup_initial_deformation(V, F, pose_label, deformation_magnitude, deformation_ratio, fixed_boundary_range, U, indices_fixed);
  fixed_point_constraints(P, 3*V.rows(), 3, indices_fixed);

  TINYAD_DEBUG_OUT("Finish setting up fixed point constraints.");

  bool redraw = false;
  std::mutex m;
  std::thread optimization_thread(
    [&]()
    {
      // Pre-compute triangle rest shapes in local coordinate systems
      std::vector<Eigen::Matrix3d> rest_shapes(F.rows());
      for (int f_idx = 0; f_idx < F.rows(); ++f_idx)
      {
        // Get 3D vertex positions
        Eigen::Vector3d ar = V.row(F(f_idx, 0));
        Eigen::Vector3d br = V.row(F(f_idx, 1));
        Eigen::Vector3d cr = V.row(F(f_idx, 2));
        Eigen::Vector3d dr = V.row(F(f_idx, 3));

        // Save 3-by-3 matrix with edge vectors as colums
        rest_shapes[f_idx] = TinyAD::col_mat(br - ar, cr - ar, dr - ar);
      };

      // Set up function with 3d vertex positions as variables.
      auto func = TinyAD::scalar_function<3>(TinyAD::range(V.rows()));

      // Add objective term per face. Each connecting 4 vertices.
      func.add_elements<4>(TinyAD::range(F.rows()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element) {
          // Evaluate element using either double or TinyAD::Double
          using T = TINYAD_SCALAR_TYPE(element);

          // Get variable 3d vertex positions
          Eigen::Index f_idx = element.handle;
          Eigen::Vector3<T> a = element.variables(F(f_idx, 0));
          Eigen::Vector3<T> b = element.variables(F(f_idx, 1));
          Eigen::Vector3<T> c = element.variables(F(f_idx, 2));
          Eigen::Vector3<T> d = element.variables(F(f_idx, 3));

          Eigen::Matrix3<T> M = TinyAD::col_mat(b - a, c - a, d - a);
          Eigen::Matrix3d Mr = rest_shapes[f_idx];
          Eigen::Matrix3<T> J = M * Mr.inverse();

          // Compute the stable Neo-Hookean energy
          double A = 0.5 * Mr.determinant();
          const double mu = MU;
          const double lambda = LAMBDA;
          auto Ic = (J.transpose() * J).trace();
          auto detF = J.determinant();
          // Here we use Eq. (6.9) from [Kim and Eberle 2022]
          // PS: Using Eq. 14 from [Smith et al. 2018] will also give similar results
          double alpha = 1.0 + mu / lambda;
          auto W = mu / 2.0 * (Ic - 3.0) + lambda / 2.0 * (detF - alpha) * (detF - alpha);
          return A * W;
      });

      // Assemble inital x vector from U matrix.
      // x_from_data(...) takes a lambda function that maps
      // each variable handle (vertex index) to its initial 2D value (Eigen::Vector2d).
      Eigen::VectorXd x = func.x_from_data([&] (int v_idx) {
          return U.row(v_idx);
      });

      TINYAD_DEBUG_OUT("Initial energy: " << func.eval(x));

      // Projected Newton
      TinyAD::LinearSolver solver;
      int max_iters = 200;
      Eigen::VectorXd d;
      Eigen::SparseMatrix<double> H;
      std::vector<double> hist;
      for (int i = 0; i < max_iters; ++i)
      {
        igl::writeOBJ(output_folder + "obj/" + mesh_name + "_" + eps_str + "/" + output_tag + "_iter_" + std::to_string(i) + ".obj", U, FF);
        igl::writeMESH(output_folder + "obj/" + mesh_name + "_" + eps_str + "/" + output_tag + "_iter_" + std::to_string(i) + ".mesh", U, F, FF);
      
        auto [f, g, H_proj] = func.eval_with_hessian_proj(x, eps);

        // record the energy
        hist.push_back(f);

        TINYAD_DEBUG_OUT("Energy in iteration " << i << ": " << f);

        // choose optimization method
        H = H_proj;
        // newton's iteration 
        Eigen::VectorXd Pg = P*g;
        Eigen::SparseMatrix<double> PHP = P*H*P.transpose();
        d = TinyAD::newton_direction(Pg, PHP, solver);
        d = P.transpose() * d;

        TINYAD_DEBUG_OUT("Newton decrement " << i << " =  " << TinyAD::newton_decrement(d, g));
        TINYAD_DEBUG_OUT("d norm " << i << " =  " << d.norm());
        TINYAD_DEBUG_OUT("g norm " << i << " =  " << g.norm());

        if (TinyAD::newton_decrement(d, g) < convergence_eps * LAMBDA)
          break;

        // line search
        x = TinyAD::line_search(x, d, f, g, func, 1.0, 0.8, 100, 1e-8);

      // Write final x vector to U matrix.
      // x_to_data(...) takes a lambda function that writes the final value
      // of each variable (Eigen::Vector2d) back to our U matrix.
        func.x_to_data(x, [&] (int v_idx, const Eigen::Vector3d& p) {
            U.row(v_idx) = p;
            });
        {
          std::lock_guard<std::mutex> lock(m);
          redraw = true; 
        }
      }

      TINYAD_DEBUG_OUT("Final energy: " << func.eval(x));
      hist.push_back(func.eval(x));

      {
        igl::writeOBJ(output_folder + "obj/" + mesh_name + "_" + eps_str + "/" + output_tag + "_iter_" + std::to_string(hist.size()-1) + ".obj", U, FF);
        igl::writeMESH(output_folder + "obj/" + mesh_name + "_" + eps_str + "/" + output_tag + "_iter_" + std::to_string(hist.size()-1) + ".mesh", U, F, FF);
        
        std::ofstream output_file_fixed(output_folder + "obj/" + mesh_name + "_" + eps_str + "/" + output_tag + "_fixed_vid.txt");
        std::ostream_iterator<int> output_iterator_fixed(output_file_fixed, "\n");
        std::copy(std::begin(indices_fixed), std::end(indices_fixed), output_iterator_fixed);

        std::ofstream output_file(output_folder + "hist/" + output_tag + ".txt");
        std::ostream_iterator<double> output_iterator(output_file, "\n");
        std::copy(std::begin(hist), std::end(hist), output_iterator);

        std::ofstream output_file_iter(output_folder + "iter/" + output_tag + ".txt");
        output_file_iter << (hist.size()-1) << std::endl;
      }
    });

  // Plot mesh
  viewer.core().is_animating = true;
  viewer.data().set_mesh(U, FF);
  viewer.core().align_camera_center(U);
  viewer.data().show_lines = false;

  viewer.callback_pre_draw = [&] (igl::opengl::glfw::Viewer& viewer)
  {
    if(redraw)
    {
      viewer.data().set_vertices(U);
      viewer.core().align_camera_center(U);
      {
        std::lock_guard<std::mutex> lock(m);
        redraw = false;
      }
    }
    return false;
  };

  viewer.launch();
  if(optimization_thread.joinable())
  {
    optimization_thread.join();
  }

  return 0;
}
