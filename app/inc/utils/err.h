/**
    @file       err.h
    @brief      error code table


    @version    0.1.0
    @date       2026-08-10
    @author     hjk

*/

#ifndef BBE_ERR_H__
#define BBE_ERR_H__

#include "NatureDSP_types.h"

typedef int32_t ErrorType;

#define ERR_OK				(0)		/* success 0 */
#define ERR_FAIL			(-1)	/* failure -1 */
#define ERR_BUSY			(1)	    /* busy 1 */
#define ERR_NOT_FOUND		(3)		/* param name not found in descriptor table 3 */
#define ERR_INVAL_PARAMS	(6)		/* invalid params 6 */
#define ERR_INVAL_STATE		(8)		/* invalid state 8 */
#define ERR_UNINIT			(12)	/* uninit 12 */
#define ERR_DATA_INTEG		(14)	/* data integrity: param file missing/format error 14 */


#endif /* BBE_ERR_H__ */
