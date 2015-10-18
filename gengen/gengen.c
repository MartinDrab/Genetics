
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <strings.h>
#include "utils.h"
#include "options.h"
#include "generator.h"
#include "gengen.h"


static ERR_VALUE _default_options(void)
{
	ERR_VALUE ret = ERR_SUCCESS;

	ret = option_add_UInt32(GENGEN_OPTION_REGION_MIN, 10);
	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GENGEN_OPTION_REGION_MAX, 10);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GENGEN_OPTION_READ_LENGTH, 5);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GENGEN_OPTION_READS_MIN, 6);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GENGEN_OPTION_READS_MAX, 6);

	if (ret == ERR_SUCCESS)
		ret = option_add_String(GENGEN_OPTION_NUCLEOTIDES, "ACGT");

	if (ret == ERR_SUCCESS)
		ret = option_add_Float(GENGEN_OPTION_INDEL_PROB, 0.0);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GENGEN_OPTION_INDELS_MIN, 0);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GENGEN_OPTION_INDELS_MAX, 0);

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(GENGEN_OPTION_DISABLE_INS, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(GENGEN_OPTION_DISABLE_DELS, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_Float(GENGEN_OPTION_REPLACE_PROB, 0.0);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GENGEN_OPTION_REPLACE_MIN, 0);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GENGEN_OPTION_REPLACE_MAX, 0);

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(GENGEN_OPTION_HELP, FALSE);

	if (ret == ERR_SUCCESS) {
		ret = option_set_description_const(GENGEN_OPTION_REGION_MIN, GENGEN_OPTION_REGION_MIN_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GENGEN_OPTION_REGION_MAX, GENGEN_OPTION_REGION_MAX_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GENGEN_OPTION_READ_LENGTH, GENGEN_OPTION_READ_LENGTH_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GENGEN_OPTION_READS_MIN, GENGEN_OPTION_READS_MIN_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GENGEN_OPTION_READS_MAX, GENGEN_OPTION_READS_MAX_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GENGEN_OPTION_NUCLEOTIDES, GENGEN_OPTION_NUCLEOTIDES_DESC);
		assert(ret == ERR_SUCCESS);

		ret = option_set_description_const(GENGEN_OPTION_INDEL_PROB, GENGEN_OPTION_INDEL_PROB_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GENGEN_OPTION_INDELS_MIN, GENGEN_OPTION_INDELS_MIN_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GENGEN_OPTION_INDELS_MAX, GENGEN_OPTION_INDELS_MAX_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GENGEN_OPTION_DISABLE_INS, GENGEN_OPTION_DISABLE_INS_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GENGEN_OPTION_DISABLE_DELS, GENGEN_OPTION_DISABLE_DELS_DESC);
		assert(ret == ERR_SUCCESS);

		ret = option_set_description_const(GENGEN_OPTION_REPLACE_PROB, GENGEN_OPTION_REPLACE_PROB_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GENGEN_OPTION_REPLACE_MIN, GENGEN_OPTION_REPLACE_MIN_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GENGEN_OPTION_REPLACE_MAX, GENGEN_OPTION_REPLACE_MAX_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GENGEN_OPTION_HELP, GENGEN_OPTION_HELP_DESC);
		assert(ret == ERR_SUCCESS);

		ret = ERR_SUCCESS;
	}

	return ret;
}

static ERR_VALUE _fill_read_options(PGENERATOR_READ_OPTIONS Options)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = option_get_UInt32(GENGEN_OPTION_READ_LENGTH, &Options->ReadLength);
	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(GENGEN_OPTION_READS_MIN, &Options->MinReads);

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(GENGEN_OPTION_READS_MAX, &Options->MaxReads);

	return ret;
}

static ERR_VALUE _fill_indel_options(PGENERATOR_INDEL_OPTIONS Options)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = option_get_UInt32(GENGEN_OPTION_INDELS_MIN, &Options->MinIndels);
	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(GENGEN_OPTION_INDELS_MAX, &Options->MaxIndels);
	
	if (ret == ERR_SUCCESS)
		ret = option_get_Float(GENGEN_OPTION_INDEL_PROB, &Options->IndelProbability);
	
	if (ret == ERR_SUCCESS)
		ret = option_get_Boolean(GENGEN_OPTION_DISABLE_INS, &Options->DisableIns);
	
	if (ret == ERR_SUCCESS)
		ret = option_get_Boolean(GENGEN_OPTION_DISABLE_DELS, &Options->DisableDels);

	return ret;
}

ERR_VALUE _fill_replace_options(PGENERATOR_REPLACE_OPTIONS Options)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = option_get_UInt32(GENGEN_OPTION_REPLACE_MIN, &Options->MinReplace);
	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(GENGEN_OPTION_REPLACE_MAX, &Options->MaxReplace);

	if (ret == ERR_SUCCESS)
		ret = option_get_Float(GENGEN_OPTION_REPLACE_PROB, &Options->ReplaceProbability);

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
				if (argc == 1 || (argc == 2 && (strcasecmp(argv[1], "--help") == 0 ||
					strcasecmp(argv[1], "-h") == 0 || strcasecmp(argv[1], "/?") == 0 ||
					strcasecmp(argv[1], "-?") == 0))) {
					printf("Usage: gengen [options]");
					options_print_help();
				} else {
					uint32_t minRegion = 0;
					uint32_t maxRegion = 0;
					char *nucleotides = NULL;

					ret = option_get_UInt32(GENGEN_OPTION_REGION_MIN, &minRegion);
					if (ret == ERR_SUCCESS)
						ret = option_get_UInt32(GENGEN_OPTION_REGION_MAX, &maxRegion);

					if (ret == ERR_SUCCESS)
						ret = option_get_String(GENGEN_OPTION_NUCLEOTIDES, &nucleotides);

					if (ret == ERR_SUCCESS) {
						char *activeRegion = NULL;
						size_t activeRegionLength = 0;

						ret = generate_active_region(minRegion, maxRegion, nucleotides, &activeRegion, &activeRegionLength);
						if (ret == ERR_SUCCESS) {
							printf("%s\n", activeRegion);

							GENERATOR_READ_OPTIONS readOptions;
							GENERATOR_INDEL_OPTIONS indelOptions;
							GENERATOR_REPLACE_OPTIONS replaceOptions;

							ret = _fill_read_options(&readOptions);
							if (ret == ERR_SUCCESS) {
								ret = _fill_indel_options(&indelOptions);
								if (ret == ERR_SUCCESS) {
									ret = _fill_replace_options(&replaceOptions);
									if (ret == ERR_SUCCESS) {
										char **reads = NULL;
										size_t readCount = 0;

										ret = generate_reads(activeRegion, activeRegionLength, nucleotides, &readOptions, &indelOptions, &replaceOptions, &reads, &readCount);
										if (ret == ERR_SUCCESS) {
											for (size_t i = 0; i < readCount; ++i)
												printf("%s\n", reads[i]);

											free_reads(reads, readCount);
										}
									}
								}
							}

							free_active_region(activeRegion, activeRegionLength);
						}
					}
				}

			}
		}

		options_module_finit();
	}

	return ret;
}
