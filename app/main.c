#include <stdio.h>

#include "system.h"


int main( int argc, char **argv )
{
	system_init();
	system_run();
	system_deinit();
	return 0;
}
