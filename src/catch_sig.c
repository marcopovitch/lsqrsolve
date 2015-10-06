#include <stdio.h>

#include "catch_sig.h"
#include "extern.h"

void emergency_halt()
{
    fprintf(stderr, "Caugth Ctrl-C signal !\n");
    fprintf(stderr, "lsqrsole will stop at the end of this iteration\n");
    please_stop_lsqr = 1;
}
