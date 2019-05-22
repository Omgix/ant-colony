#ifndef ANT_COLONY_ANTCOLONYCLASS_H
#define ANT_COLONY_ANTCOLONYCLASS_H

#include <Eigen/Eigen>

#include <fstream>
#include <string>
#include <vector>

class AntColony
{
 public:
  explicit AntColony(const char *filename);
  explicit AntColony(const std::string &filename);
  ~AntColony();
  int calcTSP();
  int recalcTSP();
  std::vector<int> &get_path();
 private:
  int _n_nodes;
  Eigen::MatrixXd _adj_matrix;
  bool _caculated;
  std::vector<int> _path;
};

#endif //ANT_COLONY_ANTCOLONYCLASS_H
