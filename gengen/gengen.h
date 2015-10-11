
#ifndef __GENGEN_H__
#define __GENGEN_H__


/************************************************************************/
/*                  PROGRAM OPTION NAMES                                */
/************************************************************************/

#define GENGEN_OPTION_REGION_MIN				"min-region"
#define GENGEN_OPTION_REGION_MAX				"max-region"
#define GENGEN_OPTION_READ_LENGTH				"read-length"
#define GENGEN_OPTION_READS_MIN					"min-reads"
#define GENGEN_OPTION_READS_MAX					"max-reads"
#define GENGEN_OPTION_NUCLEOTIDES				"nucleotides"

#define GENGEN_OPTION_INDEL_PROB				"indel-prob"
#define GENGEN_OPTION_INDELS_MIN				"min-indels"
#define GENGEN_OPTION_INDELS_MAX				"max-indels"
#define GENGEN_OPTION_DISABLE_INS				"disable-ins"
#define GENGEN_OPTION_DISABLE_DELS				"disable-dels"

#define GENGEN_OPTION_REPLACE_PROB				"replace-prob"
#define GENGEN_OPTION_REPLACE_MIN				"min-replace"
#define GENGEN_OPTION_REPLACE_MAX				"max-replace"

#define GENGEN_OPTION_HELP						"help"

/************************************************************************/
/*                  PROGRAM OPTION DESCRIPTIONS                         */
/************************************************************************/

#define GENGEN_OPTION_REGION_MIN_DESC				"Minimum active region length, in bases"
#define GENGEN_OPTION_REGION_MAX_DESC				"Maximum active region length, in bases"
#define GENGEN_OPTION_READ_LENGTH_DESC				"Length of a single read, in bases"
#define GENGEN_OPTION_READS_MIN_DESC				"Minimum number of reads to generate"
#define GENGEN_OPTION_READS_MAX_DESC				"Maximum number of reads to generate"
#define GENGEN_OPTION_NUCLEOTIDES_DESC				"Possible nucleotides within the generated data. Every character of the string is used as one nucleotide"

#define GENGEN_OPTION_INDEL_PROB_DESC				"Probability that a read contains any indels"
#define GENGEN_OPTION_INDELS_MIN_DESC				"Minimum indels in a read"
#define GENGEN_OPTION_INDELS_MAX_DESC				"Maximum indels in a read"
#define GENGEN_OPTION_DISABLE_INS_DESC				"Do not make insertions"
#define GENGEN_OPTION_DISABLE_DELS_DESC				"Do not make deletions"

#define GENGEN_OPTION_REPLACE_PROB_DESC				"Probability that a substitution occurs in a read"
#define GENGEN_OPTION_REPLACE_MIN_DESC				"Minimum number of replaced bases"
#define GENGEN_OPTION_REPLACE_MAX_DESC				"Maximum number of replaced bases"

#define GENGEN_OPTION_HELP_DESC						"See this help"


#endif 
