#include "AntColonyBase.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

AntColonyBase::AntColonyBase(const char *filename) {
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

AntColonyBase::AntColonyBase(const std::string &filename) : AntColonyBase(filename.c_str()) {}

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
  const int MMax = 9999;              //max number of ants
  const int NMax = 500;               //max citys
  int m = _dim;                       //number of ants
  const double Q = 999 ;              //flexible
  const int K = 1000;                 //ITER
  Eigen::MatrixXd Phe(_dim, _dim);    // Pheromone
  double LK;                          //total length
  Eigen::MatrixXi Path(m, _dim);   //Record the path of the ant to prevent duplicate paths. Recorded is the point
  int ant;                            //Ant's current location
  int i,j,k,p;                        //loop variables
  double _adj_matrixis = 0.1;         //Rate at which each pheromone disappears
  int sameNum,samePhe[NMax];          //Every time I go to find the side with the most pheromone, as in the initial situation,
                                      // when the amount of pheromone is the same,
  int bugNum,bugTry[NMax];            //Selection made in the event of an error
  double bugP = 0.90;                 //Selection made in the event of an error
  int start = 0;                      //Starting point, the city number is from 0 - n-1.
  double Max;                         //Used to select the most pheromone side
  bool Passed[NMax];                  //Used to determine if the city has passed, can it be selected

  for(i = 0;i < _dim;i++)
    for(j = 0; j < _dim;j++)
      Phe(i,j) = 1; //Initialize the pheromone concentration on each side
  for(i = 0;i< m;i++)
    Path(i,0) = start;  //The starting point of each ant is fixed
  for(k = 0;k < K;k++){
    for(i = 0;i < _dim;i++)
      for(j = 0; j < _dim;j++)
        Phe(i,j) *= _adj_matrixis ;//After each cycle, the pheromone disappears
    srand((unsigned)time(nullptr));
    for(i = 0;i < m;i++){//For each ant, perform a loop
      ant = start;
      for(j = 0;j < _dim;j++)
        Passed[j] = false;
      Passed[ant] = true;
      for(j = 1;j < _dim;j++){//Each ant chooses n-1 times
        Max = 0;
        sameNum  = 0 ;
        bugNum = 0;
        for(p = 0;p < _dim;p++)
          if(!Passed[p])
            Max = Max > Phe(ant,p) ? Max : Phe(ant,p) ;//Find the maximum value of the pheromone around the edge
        for(p = 0;p < _dim;p++)
          if(Max == Phe(ant,p))
            if(!Passed[p])
              samePhe[sameNum++] = p;//When the record pheromone takes the maximum value, the corresponding city number and quantity
        for(p = 0;p < _dim;p++)
          if(!Passed[p])
            bugTry[bugNum++] = p;
        if( (double)rand() /32765 < bugP)
          ant = samePhe[ rand() % sameNum ] ;//Select a side from it with a random number
        else
          ant = bugTry [ rand() % bugNum ] ;//In case of error, randomly select a side
        Passed[ant] = true;
        Path(i,j) = ant;
      }
    }
    //After completing the operation of each ant, perform the operation of adding pheromone,
    // using Ant-Circle System
    for(i = 0; i < m;i++){
      LK  = 0 ;
      for(j = 0; j < _dim-1;j++)
        LK += _adj_matrix(Path(i,j),Path(i,j+1));//Calculate the total distance of ants in a loop
      LK += _adj_matrix(Path(i,j),Path(i,0));
      for(j = 0; j < _dim-1;j++)
        Phe(Path(i,j),Path(i,j+1)) += Q/LK ;
      Phe(Path(i,j),Path(i,0)) += Q/LK ;
    }
  }


  double curr = 1e32; //Although the operation has been completed, we have to intuitively find the
  // shortest path from all existing paths.
  for(i = 0;i < m;i++){
    LK = 0;
    for(j = 0;j < _dim-1;j++)
      LK += _adj_matrix(Path(i,j),Path(i,j+1));//Calculate the total distance of ants in a loop
    LK += _adj_matrix(Path(i,j),Path(i,0));//Back to the initial point
    if(LK < curr){
      curr = LK;
      start = i;
    }
  }
  for(i = 0;i < _dim; i++)
    _path.push_back(Path(start, i));

  _caculated = true;
  return 0;
}

std::deque<int> &
AntColonyBase::get_path() {
  return _path;
}

void
AntColonyBase::printAdj(std::ostream &os) {
  for (int i = 0; i < _dim; ++i) {
    for (int j = 0; j < _dim; ++j)
      os << _adj_matrix(i,j) << ' ';
    os << std::endl;
  }
}