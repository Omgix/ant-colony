#include "AntColonyBase.h"

#include <iostream>

int main(int argc, char **argv) {
  AntColonyBase data(argv[1]);
  //data.printAdj(std::cout);
  //std::cout << std::endl;
  data.calcTSP();
  std::deque<int> &path = data.get_path();

  for (int node: path)
    std::cout << node + 1 << ' ';
  std::cout << std::endl;
  std::cout << data.total_len() << std::endl;
}