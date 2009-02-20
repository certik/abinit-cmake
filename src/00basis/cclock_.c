/*  Zachary Levine   7 Dec 1994

    Simple ANSI c routine which returns cpu time in seconds.
    Suitable for a call by a Fortran program.
    Usage (fortran): call cclock(cpu) for double precision cpu.
    Results are implementation-dependent with no particular guarantee
    of accuracy in the results.

*/

#include <time.h>

void cclock_ (cpu)

double* cpu;

{
    *cpu = ((double) clock()) / CLOCKS_PER_SEC;
}
