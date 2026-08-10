#include "utils/debug_trace.h"


#if TRACE_ENABLED

#include <stddef.h>
#include <stdio.h>

ErrorType DumpArray(const void *const aArr,
                          const ArrayKindType aKind,
                          const uint16_t aCount,
                          const char *const aLabel)
{
	ErrorType ret = ERR_OK;
	uint16_t i;

	if ((NULL == aArr) || (aKind > ARRAY_KIND_CPLX16))
	{
		ret = ERR_INVAL_PARAMS;
	}
	else
	{
		printf("%s", (NULL != aLabel) ? aLabel : "");

		for (i = 0U; i < aCount; ++i)
		{
			if (ARRAY_KIND_U8 == aKind)
			{
				printf("%u, ", (unsigned) ((const uint8_t *) aArr)[i]);
			}
			else if (ARRAY_KIND_S16 == aKind)
			{
				printf("%d, ", (int) ((const int16_t *) aArr)[i]);
			}
			else if (ARRAY_KIND_U16 == aKind)
			{
				printf("%u, ", (unsigned) ((const uint16_t *) aArr)[i]);
			}
			else if (ARRAY_KIND_S32 == aKind)
			{
				printf("%d, ", (int) ((const int32_t *) aArr)[i]);
			}
			else
			{
				printf("%d+%di, ",
				       (int) ((const complex_fract16 *) aArr)[i].s.re,
				       (int) ((const complex_fract16 *) aArr)[i].s.im);
			}
		}

		printf("\n");
	}

	return ret;
}

#endif /* TRACE_ENABLED */