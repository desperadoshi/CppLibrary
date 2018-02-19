#ifndef MISC_H
#define MISC_H

typedef double Real;

#define error_here() \
  std::cerr<<__FILE__<<", line "<<__LINE__<<", compiled "<<__DATE__<<" at " \
  <<__TIME__<<std::endl;

#endif