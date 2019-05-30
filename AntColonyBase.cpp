#include "AntColonyBase.h"

#include <cmath>
#include <chrono>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <iostream>

AntColonyBase::AntColonyBase(const char *filename, double alpha, double beta,
                             double rho, double colony_eff, unsigned maxiter):
                             _alpha(alpha), _beta(beta), _rho(rho), _colony_eff(colony_eff), _maxiter(maxiter)
{
  std::ifstream input (filename);
  NODE_COORD_TYPE type = NONE;
  _caculated = false;
  std::string line;

  // Read file to get dimension and points coord.
  while (getline(input, line)) {

    if (line.find("DIMENSION") != std::string::npos) {
      size_t pos = line.find(":");
      std::istringstream str_n_nodes (line.substr(pos + 1));
      str_n_nodes >> _dim;
    } else if (line.find("NODE_COORD_SECTION") != std::string::npos) {
      type = COORD;
      break;
    } else if (line.find("EDGE_WEIGHT_SECTION") != std::string::npos) {
      type = NO_COORDS;
      break;
    }
  }

  _mintour_each = std::deque<double>(_dim, 0);
  _mintour_global = std::deque<double>(_dim, 0);

  if (type == COORD) {

    std::deque<Point> points;
    Point point;

    for (int i = 0; i < _dim; ++i) {
      getline(input, line);
      std::istringstream coord(line);
      coord >> point.x;
      coord >> point.x;
      coord >> point.y;
      points.push_back(point);
    }

    _adj_matrix = Eigen::MatrixXd(_dim, _dim);

    for (int i = 0; i < _dim; ++i)
      for (int j = 0; j <= i; ++j) {
        double d_x = points[i].x - points[j].x;
        double d_y = points[i].y - points[j].y;
        _adj_matrix(j, i) = _adj_matrix(i, j) = sqrt(d_x * d_x + d_y * d_y);
      }

  } else if (type == NO_COORDS) {

    _adj_matrix = Eigen::MatrixXd(_dim, _dim);

    for (int i = 0; i < _dim; ++i) {
      getline(input, line);
      std::istringstream coord(line);
      double weight;
      for (int j = 0; j < _dim; ++j) {
        coord >> weight;
        _adj_matrix(i,j) = weight;
      }
    }

  } else {
    error_msg = "Unknown NODE_COORD_TYPE";
  }
}

AntColonyBase::AntColonyBase(const std::string &filename, double alpha, double beta,
                             double rho, double colony_eff, unsigned maxiter) :
                             AntColonyBase(filename.c_str(), alpha, beta, rho, colony_eff, maxiter) {}

int
AntColonyBase::calcTSP()
{
  if (!_caculated)
    return recalcTSP();
  else
    return 1;
}

int
AntColonyBase::recalcTSP()
{
  // Save path to the member _path, start at 0;
  // Return value:
  // 0 Calculate success
  // -1 Fail, max iterations reach.

  // Note: change member _caculated to true at last
  const int NMax = 500;               //max citys
  int m = ceil(_colony_eff * _dim);                       //number of ants
  const double Q = 999;              //flexible
  Eigen::MatrixXd Phe = Eigen::MatrixXd::Constant(_dim, _dim, 1);    // Pheromone
  int ant;                            //Ant's current location
  int i,j,k,p;                        //loop variables
  bool Passed[NMax];                  //Used to determine if the city has passed, can it be selected

  long seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_int_distribution<int> start_die (0, _dim - 1);
  double min_L = 1e32;
  Eigen::VectorXi min_path (_dim);

  for(k = 0;k < _maxiter;k++){

    Eigen::MatrixXd deltaPhe = Eigen::MatrixXd::Zero(_dim, _dim);
    double min_L_local = 1e32;

    //For each ant, perform a loop
    for(i = 0;i < m;i++){

      Eigen::MatrixXd deltaPheSingle = Eigen::MatrixXd::Zero(_dim, _dim);
      double LK = 0;
      Eigen::VectorXi path (_dim);
      int start = start_die(generator);
      ant = start;
      path(0) = start;
      for(j = 0;j < _dim;j++)
        Passed[j] = false;

      //Each ant chooses n-1 times to complete the tour
      for(j = 1;j < _dim;j++){
        Passed[ant] = true;
        // Calculate the probability to go to the next city.
        std::vector<double> probability;
        for(p = 0;p < _dim;p++)
          if(!Passed[p])
            probability.push_back(pow(Phe(j, p), _alpha) * pow(1/_adj_matrix(i,j), _beta));
          else
            probability.push_back(0.0);

        //for (double n: probability)
          //std::cout << n << ' ';
        //std::cout << std::endl;
        // Construct a generator with such probability;
        std::discrete_distribution<int>
            next_die (probability.begin(), probability.end());
        // Choose next city to visit;
        int next = next_die(generator);
        Passed[next] = true;
        path(j) = next;
        LK += _adj_matrix(ant, next);
        deltaPheSingle(ant, next) = Q;
        ant = next;
      }

      LK += _adj_matrix(ant, start);
      deltaPheSingle(ant, start) = Q;
      deltaPhe += deltaPheSingle / LK;

      if (LK < min_L) {
        min_L = LK;
        min_path = path;
      }
      if (LK < min_L_local)
        min_L_local = LK;
    }

    _mintour_global[k] = min_L;
    _mintour_each[k] = min_L_local;
    if (k == _maxiter - 1) {
      for (i = 0; i < _dim; ++i)
        if (min_path(i) == 0)
          break;

      for (j = 0; j < _dim; ++j)
        _path.push_back(min_path((i + j) % _dim));

    } else
      Phe = Phe * _rho + deltaPhe;//After each cycle, the pheromone disappears
  }

  _caculated = true;
  return 0;
}

std::deque<int> &
AntColonyBase::get_path() {
  return _path;
}

std::deque<double> &
AntColonyBase::get_mintour_each() {
  return _mintour_each;
}

std::deque<double> &
AntColonyBase::get_mintour_global() {
  return _mintour_global;
}

void
AntColonyBase::printAdj(std::ostream &os) {
  for (int i = 0; i < _dim; ++i) {
    for (int j = 0; j < _dim; ++j)
      os << _adj_matrix(i,j) << ' ';
    os << std::endl;
  }
}

double
AntColonyBase::total_len() {
  double result = 0.0;
  for (unsigned i = 0; i < _path.size() - 1; ++i)
    result += _adj_matrix(_path[i], _path[i + 1]);
  result += _adj_matrix(_path[_path.size() - 1], _path[0]);
  return result;
}