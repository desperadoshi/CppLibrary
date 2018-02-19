#ifndef IO_H
#define IO_H

#include <string>
#include <vector>

template<typename T>
void WriteMatrix(const std::string fname,T* mat,int n_row,int n_col,
  int n_acc=9,char sep=',',std::string header="");

template<typename T>
void ReadMatrix(const std::string fname,std::vector<T>& mat,int& n_row,int& n_col,
  char sep=',',int skiprows=0);

#endif
