#include <sys/time.h>
#include <sys/resource.h>

double etime_(tt)
float tt[2];
{
  int who;
  struct rusage used;
  who = 0;
  getrusage(who,&used);
  tt[0] = used.ru_utime.tv_sec+((used.ru_utime.tv_usec)/1000000.);
  tt[1] = used.ru_stime.tv_sec+((used.ru_stime.tv_usec)/1000000.);
  return(tt[0]+tt[1]);
}

double etime(tt)
float tt[2];
{
  int who;
  struct rusage used;
  who = 0;
  getrusage(who,&used);
  tt[0] = used.ru_utime.tv_sec+((used.ru_utime.tv_usec)/1000000.);
  tt[1] = used.ru_stime.tv_sec+((used.ru_stime.tv_usec)/1000000.);
  return(tt[0]+tt[1]);
}

double ETIME(tt)
float tt[2];
{
  int who;
  struct rusage used;
  who = 0;
  getrusage(who,&used);
  tt[0] = used.ru_utime.tv_sec+((used.ru_utime.tv_usec)/1000000.);
  tt[1] = used.ru_stime.tv_sec+((used.ru_stime.tv_usec)/1000000.);
  return(tt[0]+tt[1]);
}

