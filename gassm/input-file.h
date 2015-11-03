
#ifndef __GASSM_INPUT_FILE_H__
#define __GASSM_INPUT_FILE_H__





typedef struct _ACTIVE_REGION {
	uint64_t Offset;
	uint64_t Length;
	const char *Sequence;
} ACTIVE_REGION, *PACTIVE_REGION;



ERR_VALUE input_get_refseq(char **RefSeq, size_t *RefSeqLen);
void input_free_refseq(char *RefSeq, const size_t RefSeqLen);
ERR_VALUE input_get_reads(char ***Reads, size_t *ReadCount);
void input_free_reads(char **Reads, const size_t Count);
ERR_VALUE input_refseq_to_regions(const char *RefSeq, const size_t RefSeqLen, PACTIVE_REGION *Regions, size_t *Count);
void input_free_regions(PACTIVE_REGION Regions, const size_t Count);






#endif 
