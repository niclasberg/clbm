#include "iohelpers.h"

void read_uint(FILE * handle, unsigned int * val)
{
	read_n_uints(handle, val, 1);
}

void read_n_uints(FILE * handle, unsigned int * val, size_t num)
{
	if(fread(val, sizeof(unsigned int), num, handle) != num)
		fprintf(stderr, "Unable to read the unsigned int(s)\n");
}

void read_double(FILE * handle, double * val)
{
	read_n_doubles(handle, val, 1);
}

void read_n_doubles(FILE * handle, double * val, size_t num)
{
	if(fread(val, sizeof(double), num, handle) != num)
		fprintf(stderr, "Unable to read the double(s)\n");
}

void write_uint(FILE * handle, unsigned int val)
{
	write_n_uints(handle, &val, 1);
}

void write_n_uints(FILE * handle, const unsigned int * val, size_t num)
{
	if(fwrite(val, sizeof(unsigned int), num, handle) != num)
		fprintf(stderr, "Unable to write the unsigned int(s)\n");
}

void write_double(FILE * handle, double val)
{
	write_n_doubles(handle, &val, 1);
}

void write_n_doubles(FILE * handle, const double * val, size_t num)
{
	if(fwrite(val, sizeof(double), num, handle) != num)
		fprintf(stderr, "Unable to write the double(s)\n");
}
