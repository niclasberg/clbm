#ifndef MICRO_BC_H_
#define MICRO_BC_H_
#include "data_types.h"

int bc_zou_he_east, bc_zou_he_north, bc_zou_he_west, bc_zou_he_south,
	bc_zou_he_north_east, bc_zou_he_north_west, bc_zou_he_south_west, bc_zou_he_south_east,
	bc_regularized_east, bc_regularized_north, bc_regularized_west, bc_regularized_south,
	bc_diffusive_east, bc_diffusive_north, bc_diffusive_west, bc_diffusive_south;

int INIT_DIFFUSIVE, INIT_ZOU_HE, INIT_REGULARIZED;

void bc_init(int, LbmState *, unsigned int, unsigned int);
void bc_destroy();
void micro_bc(Node *, unsigned int);
static void bc_zou_he_straight(double *, double, double, double, int);
static void bc_regularized_straight(double *, double, double, double, int);

static void zou_he_east(Node *);
static void zou_he_north(Node *);
static void zou_he_west(Node *);
static void zou_he_south(Node *);

static void zou_he_north_east(Node *);
static void zou_he_north_west(Node *);
static void zou_he_south_east(Node *);
static void zou_he_south_west(Node *);

static void diffusive_south(Node *);
static void diffusive_north(Node *);
static void diffusive_east(Node *);
static void diffusive_west(Node *);


#endif /* MICRO_BC_H_ */
