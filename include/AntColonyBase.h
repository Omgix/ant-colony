#ifndef ANT_COLONY_ANTCOLONYBASE_H
#define ANT_COLONY_ANTCOLONYBASE_H

#include <Eigen/Eigen>

#include <fstream>
#include <omp.h>
#include <string>
#include <deque>

// NOTE: Only accept EDGE_WEIGHT_TYPE of EUC_2D

class AntColonyBase
{
 public:
  explicit AntColonyBase(const char *filename, double alpha = 15, double beta = 20,
      double rho = 0.1, double colony_eff = 1.0, unsigned maxiter = 500);
  explicit AntColonyBase(const std::string &filename, double alpha = 15, double beta = 20,
                         double rho = 0.1, double colony_eff = 1.0, unsigned maxiter = 500);
  AntColonyBase(const AntColonyBase&) = delete;
  AntColonyBase &operator=(const AntColonyBase&) = delete;
  int calcTSP();
  int recalcTSP();
  std::deque<int> &get_path();
  std::deque<double> &get_mintour_each();
  std::deque<double> &get_mintour_global();
  void printAdj(std::ostream &os);
  double total_len();
 private:
  enum NODE_COORD_TYPE {
    NONE,
    COORD,
    NO_COORDS,
  };
  struct Point {
    double x;
    double y;
  };
  double _alpha;                 // Regulate the influence of the intensity of pheromone
  double _beta;                   // Regulate the influence of visibility of city
  double _rho;                     //Rate at which each pheromone disappears
  double _colony_eff;
  unsigned _maxiter;

  int _dim;
  Eigen::MatrixXd _adj_matrix;
  bool _caculated;
  std::deque<int> _path;
  std::deque<double> _mintour_each;
  std::deque<double> _mintour_global;
  std::string error_msg;
};

#endif //ANT_COLONY_ANTCOLONYBASE_H
