
#include <stdio.h>
#include "err.h"
#include "utils.h"
#include "options.h"
#include "kmer.h"
#include "kmer-table.h"
#include "gassm.h"



static ERR_VALUE _set_default_values(void)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = option_add_UInt32(GASSM_OPTION_KMER_SIZE, 3);

	return ret;
}



int main(int argc, char *argv[])
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = options_module_init(37);
	if (ret == ERR_SUCCESS) {
		ret = _set_default_values();
		if (ret == ERR_SUCCESS) {
			ret = options_parse_command_line(argc - 1, argv + 1);
			if (ret == ERR_SUCCESS) {
				uint32_t kmerSize = 0;

				ret = option_get_UInt32(GASSM_OPTION_KMER_SIZE, &kmerSize);
				if (ret == ERR_SUCCESS) {
					PKMER_TABLE kmerTable = NULL;

					ret = kmer_table_create(kmerSize, 2, 37, &kmerTable);
					if (ret == ERR_SUCCESS) {

						kmer_table_destroy(kmerTable);
					}
				}
			}
		}

		options_module_finit();
	}

	return ret;
}
