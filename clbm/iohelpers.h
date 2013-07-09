#ifndef IOHELPERS_H_
#define IOHELPERS_H_
#include <stdio.h>
#include <string.h>

void read_uint(FILE *, unsigned int *);
void read_n_uints(FILE *, unsigned int *, size_t);
void read_double(FILE *, double *);
void read_n_doubles(FILE *, double *, size_t);

void write_uint(FILE *, unsigned int);
void write_n_uints(FILE *, const unsigned int *, size_t);
void write_double(FILE *, double);
void write_n_doubles(FILE *, const double *, size_t);

#endif /* IOHELPERS_H_ */
