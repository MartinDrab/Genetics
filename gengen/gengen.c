
#include <stdio.h>
#include "utils.h"
#include "options.h"
#include "gengen.h"


static ERR_VALUE _default_options(void)
{
	ERR_VALUE ret = ERR_SUCCESS;

	ret = option_add_UInt32(GENGEN_OPTION_REGION_MIN, 8);
	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GENGEN_OPTION_REGION_MAX, 8);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GENGEN_OPTION_READ_LENGTH, 0);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GENGEN_OPTION_READS_MIN, 100000);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GENGEN_OPTION_REGION_MAX, 200000);

	if (ret == ERR_SUCCESS)
		ret = option_add_String(GENGEN_OPTION_NUCLEOTIDES, "ACGT");

	if (ret == ERR_SUCCESS)
		ret = option_add_Float(GENGEN_OPTION_INDEL_PROB, 0.0);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GENGEN_OPTION_INDELS_MIN, 0);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GENGEN_OPTION_INDELS_MAX, 0);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt8(GENGEN_OPTION_DISABLE_INS, 0);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt8(GENGEN_OPTION_DISABLE_DELS, 0);

	if (ret == ERR_SUCCESS)
		ret = option_add_Float(GENGEN_OPTION_REPLACE_PROB, 0.0);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GENGEN_OPTION_REPLACE_MIN, 0);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GENGEN_OPTION_REPLACE_MAX, 0);

	return ret;
}


int main(int argc, char *argv[])
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = options_module_init(19);
	if (ret == ERR_SUCCESS) {
		ret = _default_options();
		if (ret == ERR_SUCCESS) {
			ret = options_parse_command_line(argc - 1, argv + 1);
			if (ret == ERR_SUCCESS) {

			}
		}

		options_module_finit();
	}

	return ret;
}