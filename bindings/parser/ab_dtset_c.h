#ifndef AB_DTSET_H
#define AB_DTSET_H

#include <stdlib.h>

/**
 * AbDtsetTypes:
 * @_INT_SCALAR: a 32 bits integer.
 * @_DOUBLE_SCALAR: a 64 bits float.
 * @_INT_ARRAY: an array of 32 bits integers.
 * @_DOUBLE_ARRAY: an array of 64 bits floats.
 *
 * The possible types of the attributes of datasets.
 */
typedef enum
  {
    _INT_SCALAR,
    _INT_ARRAY,
    _DOUBLE_SCALAR,
    _DOUBLE_ARRAY,
    _OTHER
  } AbDtsetTypes;

#include "dtset_c.h"

AbDtsetTypes ab_dtset_get_type_from_id(AbDtsetIds id);
/**
 * AB_DTSET_TYPE:
 * @A: an #AbDtsetIds id.
 *
 * Get the type of a given attribute of Dtset structure.
 *
 * Returns: a #AbDtsetTypes id.
 */
#define AB_DTSET_TYPE(A) ab_dtset_get_type_from_id(A)

/**
 * AbDtsets:
 *
 * An object to handle an array of ABINIT datasets, read from a file.
 */
typedef int AbDtsets;

/**
 * AbError:
 * @AB_NO_ERROR: no error.
 * @AB_ERROR_DTSET_OBJ: wrong dataset object.
 * @AB_ERROR_DTSET_ATT: wrong attribute in dataset.
 * @AB_ERROR_DTSET_ID: wrong dataset index.
 * @AB_ERROR_DTSET_SIZE: wrong size when accessing arrays.
 *
 * An error code.
 */
typedef enum
  {
    AB_NO_ERROR,
    AB_ERROR_DTSET_OBJ,
    AB_ERROR_DTSET_ATT,
    AB_ERROR_DTSET_ID,
    AB_ERROR_DTSET_SIZE
  } AbError;


/**
 * ab_dtset_new:
 * @filename: a string, NULL terminated.
 *
 * Parse the given file using ABINIT routines and allocate a
 * dtsets array. This array must be deallocated after use with
 * ab_dtset_free().
 *
 * Returns: an #AbDtsets object or NULL on failure.
 */
AbDtsets* ab_dtset_new(const char *filename);
/**
 * ab_dtset_new_from_string:
 * @string: a string, NULL terminated.
 *
 * Parse the given string using ABINIT routines and allocate a
 * dtsets array. This array must be deallocated after use with
 * ab_dtset_free().
 *
 * Returns: an #AbDtsets object or NULL on failure.
 */
AbDtsets* ab_dtset_new_from_string(const char *string);
/**
 * ab_dtset_free:
 * @dt: the dataset array to handle.
 *
 * Clean all allocated memory from the data set allocation.
 */
void ab_dtset_free(AbDtsets *dt);

/**
 * ab_dtset_get_ndtset:
 * @dt: the dataset array to handle.
 * @ndtset: a location to store the returned value.
 *
 * An array of datasets may contain more than one. Test it with this
 * routine. @ndtset will contains the number of allocated datasets (in
 * addition to the default one).
 *
 * Returns: #AB_NO_ERROR if @dt is valid and correctly parsed.
 */
AbError ab_dtset_get_ndtset(AbDtsets *dt, int *ndtset);
/**
 * ab_dtset_get_integer:
 * @dt: the dataset array to handle.
 * @id: an attribute id, see dtset_c.h.
 * @idtset: the number of the dtset to read, 0 is default value.
 * @value: a location to store the returned value.
 *
 * Use this method to get the value of an integer attribute. @idtset
 * must be in [0;n] where n is the returned value of
 * ab_dtset_get_ndtset(). If @id is unknown, return value is
 * 0. For real attributes, see ab_dtset_get_real().
 *
 * Returns: #AB_NO_ERROR if values are correctly read.
 */
AbError ab_dtset_get_integer(AbDtsets *dt, AbDtsetIds id, int idtset, int *value);
/**
 * ab_dtset_get_real:
 * @dt: the dataset array to handle.
 * @id: an attribute id, see dtset_c.h ;
 * @idtset: the number of the dtset to read, 0 is default.
 * @value: a location to store the returned value.
 *
 * Use this method to get the value of a double attribute. @idtset
 * must be in [0;n] where n is the return value of
 * ab_dtset_get_ndtset(). If @id is unknown, return value is
 * undefined. For integer attributes, see ab_dtset_get_integer().
 *
 * Returns: #AB_NO_ERROR if values are correctly read.
 */
AbError ab_dtset_get_real(AbDtsets *dt, AbDtsetIds id, int idtset, double *value);

/**
 * ab_dtset_get_shape:
 * @dt: the dataset array to handle.
 * @n: a location to store the number of dimensions.
 * @dims: an array with 7 integers ;
 * @id: an attribute id, see dtset_c.h ;
 * @idtset: the number of the dtset to read, 0 is default.
 *
 * This method is used to poll the size of an array attribute. The
 * shape of the attribute is stored in @dims. Only the @n first values
 * of @dims are relevant.
 *
 * Returns: #AB_NO_ERROR if values are correctly read.
 */
AbError ab_dtset_get_shape(AbDtsets *dt, int *n, int dims[7],
			   AbDtsetIds id, int idtset);
/**
 * ab_dtset_get_integer_array:
 * @dt: the dataset array to handle.
 * @values: an allocated array of @n values ;
 * @n: the size of the given array ;
 * @id: an attribute id, see dtset_c.h ;
 * @idtset: the number of the dtset to read, 0 is default.
 *
 * This method is used to read the values of an array. The array must
 * already be allocated. To know its size, use ab_dtset_get_shape().
 *
 * Returns: #AB_NO_ERROR if values are correctly read.
 */
AbError ab_dtset_get_integer_array(AbDtsets *dt, int *values, size_t n,
				   AbDtsetIds id, int idtset);
/**
 * ab_dtset_get_real_array:
 * @dt: the dataset array to handle.
 * @values: an allocated array of @n values ;
 * @n: the size of the given array ;
 * @id: an attribute id, see dtset_c.h ;
 * @idtset: the number of the dtset to read, 0 is default.
 *
 * This method is used to read the values of an array. The array must
 * already be allocated. To know its size, use ab_dtset_get_shape().
 *
 * Returns: #AB_NO_ERROR if values are correctly read.
 */
AbError ab_dtset_get_real_array(AbDtsets *dt, double *values, size_t n,
				AbDtsetIds id, int idtset);

#endif
