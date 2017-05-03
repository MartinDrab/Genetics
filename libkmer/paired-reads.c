
#include <stdio.h>
#include <stdint.h>
#include "err.h"
#include "utils.h"
#include "khash.h"
#include "reads.h"
#include "paired-reads.h"


KHASH_MAP_INIT_STR(RP, PPOINTER_ARRAY_ONE_READ)


static khash_t(RP) *_table = NULL;




ERR_VALUE paired_reads_init(void)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	_table = kh_init(RP);
	ret = (_table != NULL) ? ERR_SUCCESS : ERR_OUT_OF_MEMORY;

	return ret;
}


void paired_reads_finit(void)
{
	PPOINTER_ARRAY_ONE_READ pairedReads = NULL;

	for (khiter_t it = kh_begin(_table); it != kh_end(_table); ++it) {
		if (kh_exist(_table, it)) {
			pairedReads = kh_val(_table, it);
			pointer_array_finit_ONE_READ(pairedReads);
			utils_free(pairedReads);
			kh_val(_table, it) = NULL;
		}
	}


	kh_destroy(RP, _table);
	_table = NULL;

	return;
}


ERR_VALUE paired_reads_insert(const ONE_READ *Read)
{
	khiter_t it = kh_end(_table);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PPOINTER_ARRAY_ONE_READ pairedReads = NULL;

	ret = ERR_SUCCESS;
	if (Read->Extension->TemplateName != NULL && *Read->Extension->TemplateName != '\0') {
		it = kh_get(RP, _table, Read->Extension->TemplateName);
		if (it == kh_end(_table)) {
			ret = utils_malloc(sizeof(POINTER_ARRAY_ONE_READ), &pairedReads);
			if (ret == ERR_SUCCESS) {
				int tmp = 0;

				pointer_array_init_ONE_READ(pairedReads, 140);
				it = kh_put(RP, _table, Read->Extension->TemplateName, &tmp);
				switch (tmp) {
					case -1: ret = ERR_OUT_OF_MEMORY;  break;
					case 0: ret = ERR_ALREADY_EXISTS; break;
					case 2: ret = ERR_NOT_FOUND; break;
				}

				if (ret == ERR_SUCCESS)
					kh_value(_table, it) = pairedReads;
			}
		} else pairedReads = kh_value(_table, it);

		if (ret == ERR_SUCCESS)
			ret = pointer_array_push_back_ONE_READ(pairedReads, Read);
	}

	return ret;
}


ERR_VALUE paired_reads_insert_array(const ONE_READ *Reads, const size_t Count)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const ONE_READ *r = Reads;

	ret = ERR_SUCCESS;
	for (size_t i = 0; i < Count; ++i) {
		ret = paired_reads_insert(r);
		if (ret != ERR_SUCCESS)
			break;

		++r;
	}

	return ret;
}

ERR_VALUE paired_reads_first(khiter_t *Iterator, PPOINTER_ARRAY_ONE_READ *Reads)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_NO_MORE_ENTRIES;
	for (khiter_t it = kh_begin(_table); it != kh_end(_table); ++it) {
		if (kh_exist(_table, it)) {
			*Iterator = it;
			*Reads = kh_value(_table, it);
			ret = ERR_SUCCESS;
			break;
		}
	}

	return ret;
}


ERR_VALUE paired_reads_next(khiter_t Iterator, khiter_t *NewIt, PPOINTER_ARRAY_ONE_READ *Reads)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_NO_MORE_ENTRIES;
	for (khiter_t it = Iterator + 1; it != kh_end(_table); ++it) {
		if (kh_exist(_table, it)) {
			*NewIt = it;
			*Reads = kh_value(_table, it);
			ret = ERR_SUCCESS;
			break;
		}
	}

	return ret;
}


void paired_reads_fix_overlaps(boolean Strip)
{
	khiter_t iter;
	ERR_VALUE err = ERR_INTERNAL_ERROR;
	PPOINTER_ARRAY_ONE_READ reads = NULL;
	size_t totalOverlaps = 0;
	size_t mismatches = 0;

	err = paired_reads_first(&iter, &reads);
	while (err == ERR_SUCCESS) {
		for (size_t i = 0; i < pointer_array_size(reads); ++i) {			
			for (size_t j = 0; j < pointer_array_size(reads); ++j) {
				if (i == j)
					continue;
				
				PONE_READ r1 = reads->Data[i];
				PONE_READ r2 = reads->Data[j];

				if (in_range(r1->Pos, r1->ReadSequenceLen, r2->Pos) &&
					!in_range(r1->Pos, r1->ReadSequenceLen, r2->Pos + r2->ReadSequenceLen)) {
					size_t mismatchCount = 0;
					boolean matches = FALSE;
					size_t overlapLength = r1->Pos + r1->ReadSequenceLen - r2->Pos;
					const char *or1 = r1->ReadSequence + r1->ReadSequenceLen - overlapLength;
					const uint8_t *oq1 = r1->Quality + r1->ReadSequenceLen - overlapLength;
					const char *or2 = r2->ReadSequence;
					const uint8_t *oq2 = r2->Quality;

					++totalOverlaps;
					if (overlapLength > r2->ReadSequenceLen)
						overlapLength = r2->ReadSequenceLen;

					matches = (strncmp(or1, or2, overlapLength) == 0);
					r1->NoEndStrip = matches;
					if (!matches) {
						++mismatches;
						for (size_t k = 0; k < overlapLength; ++k) {
							if (or1[k] != or2[k]) {
								++mismatchCount;
							}
						}

						if (Strip) {
							if (overlapLength < r2->ReadSequenceLen) {
								size_t r2Move = overlapLength;
								const char *tmp = r2->ReadSequence + overlapLength;
								const char *tmp2 = or1 + overlapLength;

								while (*tmp == *tmp2) {
									--r2Move;
									--tmp;
									--tmp2;
								}

								r1->ReadSequenceLen -= r2Move;
								r1->ReadSequence[r1->ReadSequenceLen] = '\0';
								r1->Quality[r1->ReadSequenceLen] = '\0';
							}

							r1->ReadSequenceLen -= (r1->Pos + r1->ReadSequenceLen - r2->Pos);
							while (*or1 == *or2) {
								++or1;
								++or2;
								++r1->ReadSequenceLen;
							}

							r1->ReadSequence[r1->ReadSequenceLen] = '\0';
						}
					}
				}
			}
		}

		err = paired_reads_next(iter, &iter, &reads);
	}

	fprintf(stderr, "Overlaps: %Iu, Mismatching: %Iu\n", totalOverlaps, mismatches);

	return;
}


void paired_reads_print(FILE *Stream)
{
	khiter_t it;
	PPOINTER_ARRAY_ONE_READ reads = NULL;
	ERR_VALUE err = ERR_INTERNAL_ERROR;
	size_t totalCount = 0;
	GEN_ARRAY_size_t counts;

	dym_array_init_size_t(&counts, 140);
	err = paired_reads_first(&it, &reads);
	while (err == ERR_SUCCESS) {
		if (pointer_array_size(reads) > gen_array_size(&counts)) {
			dym_array_reserve_size_t(&counts, pointer_array_size(reads));
			counts.ValidLength = pointer_array_size(reads);
			counts.Data[pointer_array_size(reads) - 1] = 0;
		}

		counts.Data[pointer_array_size(reads) - 1] += 1;
		totalCount += pointer_array_size(reads);
		fprintf(Stream, "%s\t --> %Iu reads\n", reads->Data[0]->Extension->TemplateName, pointer_array_size(reads));
		err = paired_reads_next(it, &it, &reads);
	}

	for (size_t i = 0; i < gen_array_size(&counts); ++i)
		fprintf(Stream, "%Iu, %Iu\n", i + 1, counts.Data[i]);

	dym_array_finit_size_t(&counts);
	fprintf(Stream, "%Iu reads, %u groups\n", totalCount, _table->n_occupied);

	return;
}
