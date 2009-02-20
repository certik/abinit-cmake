/* errors.h */

/*
 * Copyright (c) 2008 ABINIT Group (MG)
 * All rights reserved.
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#if !defined HAVE_ERROR_HANDLERS 
/* Bulletproof macros for emergency cases! 
 * They might be useful if __FILE__ expands to the full path name and exceedss 
 * the max number of columns in a Fortran file. ISO doesn't define any standard!
 * */ 
#define ABI_COMMENT(MSG) call wrap_wrtout(std_out,MSG,"COLL") 
#define ABI_WARNING(MSG) call wrap_wrtout(std_out,MSG,"COLL")

/**#ifdef FC_PGI6**/
#if 0
#define ABI_ERROR(MSG)   call leave_new("COLL")
#define ABI_BUG(MSG)     call leave_new("COLL")
#else
#define ABI_ERROR(MSG)   call wrap_wrtout(std_out,MSG,"COLL") ; call leave_new("COLL")
#define ABI_BUG(MSG)     call wrap_wrtout(std_out,MSG,"COLL") ; call leave_new("COLL")
#endif

#define MEM_CHECK(ios)   if (ios/=0) call leave_new("MEM ERROR","COLL")
#define IO_CHECK(fname,unit,ios) if (ios/=0) call leave_new("Error during IO","COLL")

#else

/* Macro for basic messages (COLLECTIVE version) */
#define ABI_COMMENT(MSG) call msg_hndl(MSG, __FILE__, __LINE__, "COMMENT") 
#define ABI_WARNING(MSG) call msg_hndl(MSG, __FILE__, __LINE__, "WARNING")
#define ABI_ERROR(MSG)   call msg_hndl(MSG, __FILE__, __LINE__, "ERROR") 
#define ABI_BUG(MSG)     call msg_hndl(MSG, __FILE__, __LINE__, "BUG") 

/* Macro to check memory */
/*#define ABI_MEMCHECK(ios) if (ios/=0) call merr_hndl(ios,__FILE__, __LINE__) */

/* I/O check */
#define ABI_IOCHECK(unit,ios) if (ios/=0) call io_hndl(ios,unit, __FILE__, __LINE__)

#endif

