#ifndef _OPTIONS_H_
#define _OPTIONS_H_

#include <stdlib.h>
//#include <string>
//#include <cstring>
#include <stdio.h>

extern int PFMODE;//boltzman sampling or stochastic sampling and  partition function mode or dS mode, it is mode as defined and used for partition function of sfold
extern int NODANGLEMODE;//no dangling at all means d0
extern int D2MODE;//d2 mode
extern int DEFAULTMODE;//default mode

extern char seqfile[200];


#ifdef __cplusplus
extern "C" {
#endif
void help();
void parse_options(int argc, char** argv);

#ifdef __cplusplus
}
#endif



#endif
