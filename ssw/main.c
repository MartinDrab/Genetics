
#include <stdio.h>
#include "ssw.h"









int main(int argc, char *argv[])
{
	size_t opStringSize = 0;
	char *opString = NULL;
	ERR_VALUE ret;

	ret = utils_allocator_init(1);
	if (ret == ERR_SUCCESS) {
		ret = ssw_clever(argv[1], strlen(argv[1]), argv[2], strlen(argv[2]), 2, -1, -1, &opString, &opStringSize);
		if (ret == ERR_SUCCESS) {
			fprintf(stdout, "%.*s\n", (int)opStringSize, opString);
			utils_free(opString);
		}
	}

	return ret;
}
