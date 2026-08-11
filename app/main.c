#include <stdio.h>

#include "system.h"


int main( int argc, char **argv )
{
	if (ERR_OK != system_init())
	{
		printf("system_init failed\n");
		return -1;
	}

	system_run();
	system_deinit();
	return 0;
}
