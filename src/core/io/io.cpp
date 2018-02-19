#include "io.h"
#include "misc.h"
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

template<typename T>
void WriteMatrix(const std::string fname,T* mat,int n_row,int n_col,
  int n_acc,char sep,std::string header){
  std::ofstream ofs;
  ofs.open(fname,std::ofstream::out);
  ofs<<header<<std::endl;
  // The width 16: 1(sign bit)+1(positive)+1(dot)+9(n_acc)+4(exponent)
  ofs<<std::scientific<<std::setprecision(n_acc);
  for(int i=0;i<n_row;i++){
    for(int j=0;j<n_col-1;j++){
      ofs<<std::setw(16)<<mat[i*n_col+j]<<sep;
    }
    ofs<<std::setw(16)<<mat[i*n_col+n_col-1]<<std::endl;
  }
  ofs.close();
} // WriteMatrix

template<typename T>
void ReadMatrix(const std::string fname,std::vector<T>& mat,int& n_row,int& n_col,
  char sep,int skiprows){
  std::ifstream ifs;
  ifs.open(fname,std::ifstream::in);
  std::string line;
  std::string ele;
  n_row=0;
  int count=0;
  while(std::getline(ifs,line)){
    if(n_row>=skiprows){
      std::istringstream iss(line);
      while(std::getline(iss,ele,sep)){
        if(ele.size()==0){
          continue;
        }
        std::istringstream ele_iss(ele);
        T val;
        ele_iss>>val;
        mat.push_back(val);
      }
    }
    n_row++;
  }
  count=mat.size();
  n_row-=skiprows;
  n_col=count/n_row;
  ifs.close();
} // ReadMatrix

template void WriteMatrix(const std::string fname,Real* mat,int n_row,int n_col,
  int n_acc,char sep,std::string header);
template void ReadMatrix(const std::string fname,std::vector<Real>& mat,
  int& n_row,int& n_col,char sep,int skiprows);
