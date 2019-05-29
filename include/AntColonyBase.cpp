#include <iostream>
#include <fstream>

// NOTE: Only accept EDGE_WEIGHT_TYPE of EUC_2D


using namespace std;
//using 10 ants
const int MMax = 9999;//max number of ants
const int NMax = 500;//max citys
int m,n;// number of ants and city
const double Q = 999 ;//flixible
const int K = 1000;//ITER
int D[NMax][NMax];//legth
double Phe[NMax][NMax];//xinxisu
int LK;//totla length
int Path[MMax][NMax];//Record the path of the ant to prevent duplicate paths. Recorded is the point
int ant;//Ant's current location
int i,j,k,p;//recycle
double Dis = 0.1;//Rate at which each pheromone disappears
int sameNum,samePhe[NMax];//Every time I go to find the side with the most pheromone, as in the initial situation, when the amount of pheromone is the same,
int bugNum,bugTry[NMax];//Selection made in the event of an error
double bugP = 0.90;//Selection made in the event of an error
int start;//Starting point, the city number is from 0 - n-1.
double Max;//Used to select the most pheromone side
bool Passed[NMax];//Used to determine if the city has passed, can it be selected
int main()
{
  /*
    fstream f("data.txt",ios::in);
    f >> n;
    if( n > NMax)//本来想直接赋值，后来又决定从文件中读入。
        return 0;
    for(i = 0;i<n;i++)
        for(j = 0;j<n ;j++)
            if(j != i)
                f>>D[i][j];
    for(i = 0;i < n;i++)
        D[i][i] = 0;
    f >> start;
    if(start > n-1)return 0;//没必要的检测，如果文件未正确书写，发生意外的错误，概不负责。
    f.close();

    */
  n=3;
  int D[4][4]={
          {0,3,6,7},
          {5,0,2,3},
          {6,4,0,2},
          {3,7,5,0}
  };
  for(i = 0;i < n;i++)
    for(j = 0; j < n;j++)
      Phe[i][j] = 1;//Initialize the pheromone concentration on each side
  for(i = 0;i< m;i++)
    Path[i][0] = start;//The starting point of each ant is fixed
  m = 999;
  for(k = 0;k < K;k++){
    for(i = 0;i < n;i++)
      for(j = 0; j < n;j++)
        Phe[i][j] *= Dis ;//After each cycle, the pheromone disappears
    srand((unsigned)time(NULL));
    for(i = 0;i < m;i++){//For each ant, perform a loop
      ant = start;
      for(j = 0;j < n;j++)
        Passed[j] = false;
      Passed[ant] = true;
      for(j = 1;j < n;j++){//Each ant chooses n-1 times
        Max = 0;
        sameNum  = 0 ;
        bugNum = 0;
        for(p = 0;p < n;p++)
          if(!Passed[p])
            Max = Max > Phe[ant][p] ? Max : Phe[ant][p] ;//Find the maximum value of the pheromone around the edge
        for(p = 0;p < n;p++)
          if(Max == Phe[ant][p])
            if(!Passed[p])
              samePhe[sameNum++] = p;//When the record pheromone takes the maximum value, the corresponding city number and quantity
        for(p = 0;p < n;p++)
          if(!Passed[p])
            bugTry[bugNum++] = p;
        if( (double)rand() /32765 < bugP)
          ant = samePhe[ rand() % sameNum ] ;//Select a side from it with a random number
        else
          ant = bugTry [ rand() % bugNum ] ;//In case of error, randomly select a side
        Passed[ant] = true;
        Path[i][j] = ant;
      }
    }
    //After completing the operation of each ant, perform the operation of adding pheromone, using Ant-Circle System
    for(i = 0; i < m;i++){
      LK  = 0 ;
      for(j = 0; j < n-1;j++)
        LK += D[Path[i][j]][Path[i][j+1]];//计Calculate the total distance of ants in a loop
      LK += D[Path[i][j]][Path[i][0]];
      for(j = 0; j < n-1;j++)
        Phe[Path[i][j]][Path[i][j+1]] += Q/LK ;
      Phe[Path[i][j]][Path[i][0]] += Q/LK ;
    }
  }
  p = 0xfffff;//Although the operation has been completed, we have to intuitively find the shortest path from all existing paths.
  for(i = 0;i < m;i++){
    LK = 0;
    for(j = 0;j < n-1;j++)
      LK += D[Path[i][j]][Path[i][j+1]];//Calculate the total distance of ants in a loop
    LK += D[Path[i][j]][Path[i][0]];//Back to the initial point
    if(LK < p){
      p = LK;
      start = i;
    }
  }
  for(i = 0;i < n; i++)
    cout << Path[start][i]<<"->";
  cout << Path[start][0]<<endl;
  return 0;
}