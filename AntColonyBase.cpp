#include "AntColonyBase.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

AntColonyBase::AntColonyBase(const char *filename) {
  std::ifstream input (filename);
  bool node_coord_sec = false;
  std::deque<Point> points;
  _caculated = false;

  for (std::string line; getline(input, line);) {
    if (node_coord_sec) {

      std::istringstream coord (line);
      Point point;

      coord >> point.x;
      coord >> point.x;
      coord >> point.y;
      points.push_back(point);

      for (int i = 1; i < _dim; ++i) {
        getline(input, line);
        std::istringstream coord (line);
        coord >> point.x;
        coord >> point.x;
        coord >> point.y;
        points.push_back(point);
      }

      break;
    } else {

      if (line.find("DIMENSION") != std::string::npos) {
        size_t pos = line.find(":");
        std::istringstream str_n_nodes (line.substr(pos + 1));
        str_n_nodes >> _dim;
      } else if (line.find("NODE_COORD_SECTION") != std::string::npos) {
        node_coord_sec = true;
      }

    }
  }

  _adj_matrix = Eigen::MatrixXd(_dim, _dim);
  for (int i = 0; i < _dim; ++i)
    for (int j = 0; j < _dim; ++i) {
      double d_x = points[i].x - points[j].x;
      double d_y = points[i].y - points[j].y;
      _adj_matrix(i,j) = sqrt(d_x*d_x + d_y*d_y);
    }

}

AntColonyBase::AntColonyBase(const std::string &filename) : AntColonyBase(filename.c_str()) {}

int
AntColonyBase::calcTSP()
{
  if (!_caculated)
    recalcTSP();
}

int
AntColonyBase::recalcTSP()
{
  // Save path to the member _path, start at 0;
  // Return value:
  // 0 Calculate success
  // -1 Fail, max iterations reach.

  // Note: change member _caculated to true at last
  // TODO
}

std::deque<int> &
AntColonyBase::get_path() {
  return _path;
}