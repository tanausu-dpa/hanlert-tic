#include <unistd.h>

// This lonley bit of C code serves one purpose: to get the process ID
// from the OS.
//
// It was a suggestion from Brian Vanderwende (NCAR CISL):
//
// "Unfortunately, there is no single method to get the PID in Fortran
// that works across compilers (it's easy in GCC but not in Intel or
// PGI), so the simplest approach would be to write a simple C
// interface that gets called in Fortran"
void fgetpid_(int *id)
{
  *id = (int)getpid();
}
