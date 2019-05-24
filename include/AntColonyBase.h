#ifndef ANT_COLONY_ANTCOLONYBASE_H
#define ANT_COLONY_ANTCOLONYBASE_H

#include <Eigen/Eigen>

#include <fstream>
#include <string>
#include <deque>

// NOTE: Only accept EDGE_WEIGHT_TYPE of EUC_2D

class AntColonyBase
{
 public:
  explicit AntColonyBase(const char *filename);
  explicit AntColonyBase(const std::string &filename);
  AntColonyBase(const AntColonyBase&) = delete;
  AntColonyBase &operator=(const AntColonyBase&) = delete;
  int calcTSP();
  int recalcTSP();
  std::deque<int> &get_path();
 private:
  struct Point {
    double x;
    double y;
  };
  int _dim;
  Eigen::MatrixXd _adj_matrix;
  bool _caculated;
  std::deque<int> _path;
};

#endif //ANT_COLONY_ANTCOLONYBASE_H
