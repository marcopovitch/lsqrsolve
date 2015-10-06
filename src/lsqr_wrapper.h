#include <stdio.h>
#include <stdlib.h>

#ifdef __APPLE__
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

#include <assert.h>
#include <math.h>

#include <sparse/sparse.h>

#include "lsqr.h"

#ifndef __LSQR_WRAPPER_H__
#define __LSQR_WRAPPER_H__

extern const char *lsqr_msg[];

void sparseMATRIXxVECTOR(long int mode, dvec * x, dvec * y, void *data);
void MATRIXxVECTOR(long int mode, dvec * x, dvec * y, void *data);
#endif
