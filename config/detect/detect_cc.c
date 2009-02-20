#include <stdlib.h>
#include <stdio.h>

int main(void)
{
#if defined __GNUC__
#if defined __i386
#define arch "x86_32"
#elif defined __ia64
#define arch "IA64"
#else
#define arch "UNKNOWN"
#endif
#if defined __linux__
#define system "LINUX"
#elif defined __macosx__
#define system "MACOSX"
#else
#define system "UNKNOWN"
#endif
 printf("GCC %d.%d %s %s\n",__GNUC__,__GNUC_MINOR__,arch,system);
#endif

 return(0);
}
