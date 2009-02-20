#include "ab_dtset_c.h"
#include <string.h>
#include <stdlib.h>

/* The way the fortran compiler set the names of module routines. */
#define FC_MOD_NAME(A) __ab_dtset__ ## A
/* g95 : #define FC_MOD_NAME(A) ab_dtset_MP_ ## A */
#define FC_MOD_CALL(A,...) FC_MOD_NAME(A)(__VA_ARGS__)

/* Fortran interface. */
void FC_MOD_NAME(ab_dtset_new)(int *dt, const char *filename, int *len);
void FC_MOD_NAME(ab_dtset_new_from_string)(int *dt, const char *string, int *len);
void FC_MOD_NAME(ab_dtset_free)(int *dt);
void FC_MOD_NAME(ab_dtset_get_ndtset)(int *dt, int *ndtset, uint *errno);
void FC_MOD_NAME(ab_dtset_get_integer)(int *dt, int *value, int *att, int *idtset, uint *errno);
void FC_MOD_NAME(ab_dtset_get_real)(int *dt, double *value, int *att, int *idtset, uint *errno);
void FC_MOD_NAME(ab_dtset_get_shape)(int *dt, int *dims, int *size, int *att, int *idtset, uint *errno);
void FC_MOD_NAME(ab_dtset_get_integer_array)(int *dt, int *values, int *size, int *att, int *idtset, uint *errno);
void FC_MOD_NAME(ab_dtset_get_real_array)(int *dt, double *values, int *size, int *att, int *idtset, uint *errno);

#include "dtset_c.static.h"

AbDtsetTypes ab_dtset_get_type_from_id(AbDtsetIds id)
{
  return ((id >= 0 && id < AB_DTSET_N_IDS)?ab_dtsets_types[id]:_OTHER);
}

AbDtsets* ab_dtset_new(const char *filename)
{
  int n;
  AbDtsets *dt, tmpDt;
  
  n = strlen(filename);
  
  FC_MOD_CALL(ab_dtset_new, &tmpDt, filename, &n);
  if (tmpDt > 0)
    {
      dt = malloc(sizeof(AbDtsets));
      *dt = tmpDt;
    }
  else
    dt = (AbDtsets*)0;
  return dt;
}
AbDtsets* ab_dtset_new_from_string(const char *string)
{
  int n;
  AbDtsets *dt, tmpDt;
  
  n = strlen(string);
  
  FC_MOD_CALL(ab_dtset_new_from_string, &tmpDt, string, &n);
  if (tmpDt > 0)
    {
      dt = malloc(sizeof(AbDtsets));
      *dt = tmpDt;
    }
  else
    dt = (AbDtsets*)0;
  return dt;
}

void ab_dtset_free(AbDtsets *dt)
{
  FC_MOD_CALL(ab_dtset_free, dt);
  free(dt);
}

AbError ab_dtset_get_ndtset(AbDtsets *dt, int *ndtset)
{
  AbError errno;

  FC_MOD_CALL(ab_dtset_get_ndtset, dt, ndtset, &errno);
  return errno;
}

AbError ab_dtset_get_integer(AbDtsets *dt, AbDtsetIds id, int idtset, int *value)
{
  AbError errno;

  if (AB_DTSET_TYPE(id) != _INT_SCALAR) return AB_ERROR_DTSET_ATT;
  FC_MOD_CALL(ab_dtset_get_integer, dt, value, (int*)&id, &idtset, &errno);
  return errno;
}

AbError ab_dtset_get_real(AbDtsets *dt, AbDtsetIds id, int idtset, double *value)
{
  AbError errno;

  if (AB_DTSET_TYPE(id) != _DOUBLE_SCALAR) return AB_ERROR_DTSET_ATT;
  FC_MOD_CALL(ab_dtset_get_real, dt, value, (int*)&id, &idtset, &errno);
  return errno;
}

AbError ab_dtset_get_shape(AbDtsets *dt, int *n, int dims[7], AbDtsetIds id, int idtset)
{
  AbError errno;

  FC_MOD_CALL(ab_dtset_get_shape, dt, dims, n, (int*)&id, &idtset, &errno);
  return errno;
}

AbError ab_dtset_get_integer_array(AbDtsets *dt, int *values, size_t n, AbDtsetIds id, int idtset)
{
  AbError errno;

  if (AB_DTSET_TYPE(id) != _INT_ARRAY) return AB_ERROR_DTSET_ATT;
  FC_MOD_CALL(ab_dtset_get_integer_array, dt, values, (int*)&n, (int*)&id, &idtset, &errno);
  return errno;
}

AbError ab_dtset_get_real_array(AbDtsets *dt, double *values, size_t n, AbDtsetIds id, int idtset)
{
  AbError errno;

  if (AB_DTSET_TYPE(id) != _DOUBLE_ARRAY) return AB_ERROR_DTSET_ATT;
  FC_MOD_CALL(ab_dtset_get_real_array, dt, values, (int*)&n, (int*)&id, &idtset, &errno);
  return errno;
}
