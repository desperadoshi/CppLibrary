#ifndef MISC_H
#define MISC_H

#include <iostream>

typedef double Real;

#define error_here(IN_STR) \
  std::cerr<<#IN_STR<<" at line "<<__LINE__<<" in "<<__FILE__<<", compiled " \
  <<__DATE__<<" at "<<__TIME__<<std::endl;

#define not_implemented(IN_STR) \
  std::cout<<#IN_STR<<" is not implemented!"<<std::endl;

#endif
