
#include <inttypes.h>
#include "err.h"
#include "utils.h"
#include "refseq-storage.h"
#include "gen_dym_array.h"
#include "found-sequence.h"



ERR_VALUE variant_call_init(const char *Chrom, uint64_t Pos, const char *ID, const char *Ref, size_t RefLen, const char *Alt, size_t AltLen, const uint8_t Qual, const GEN_ARRAY_size_t *RefReads, const GEN_ARRAY_size_t *AltReads, PVARIANT_CALL VC)
{
	char *tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	
	while ((RefLen >= 2 && AltLen >= 2) && *Ref == *Alt) {
		++Pos;
		++Ref;
		++Alt;
		--RefLen;
		--AltLen;
	}

	VC->Valid = TRUE;
	VC->PhasedPos = 0;
	VC->PhaseType = vcptNone;
	ret = utils_copy_string(Chrom, &tmp);
	if (ret == ERR_SUCCESS) {
		VC->Chrom = tmp;
		VC->Pos = Pos;
		VC->Qual = Qual;
		ret = utils_copy_string(ID, &tmp);
		if (ret == ERR_SUCCESS) {
			VC->ID = tmp;
			ret = utils_calloc(RefLen + 1, sizeof(char), &tmp);
			if (ret == ERR_SUCCESS) {
				memcpy(tmp, Ref, RefLen*sizeof(char));
				tmp[RefLen] = '\0';
				VC->Ref = tmp;
				ret = utils_calloc(AltLen + 1, sizeof(char), &tmp);
				if (ret == ERR_SUCCESS) {
					memcpy(tmp, Alt, AltLen*sizeof(char));
					tmp[AltLen] = '\0';
					VC->Alt = tmp;
					dym_array_init_size_t(&VC->RefReads, 140);
					dym_array_init_size_t(&VC->AltReads, 140);
					ret = dym_array_reserve_size_t(&VC->RefReads, GEN_ARRAY_STATIC_ALLOC + 1);
					if (ret == ERR_SUCCESS)
						ret = dym_array_reserve_size_t(&VC->AltReads, GEN_ARRAY_STATIC_ALLOC + 1);

					if (ret == ERR_SUCCESS && RefReads != NULL)
						ret = dym_array_push_back_array_size_t(&VC->RefReads, RefReads);

					if (ret == ERR_SUCCESS && AltReads != NULL)
						dym_array_push_back_array_size_t(&VC->AltReads, AltReads);
						
					if (ret != ERR_SUCCESS) {
						dym_array_finit_size_t(&VC->AltReads);
						dym_array_finit_size_t(&VC->RefReads);
						utils_free(VC->Alt);
					}
				}
				
				if (ret != ERR_SUCCESS)
					utils_free(VC->Ref);
			}

			if (ret != ERR_SUCCESS)
				utils_free(VC->ID);
		}

		if (ret != ERR_SUCCESS)
			utils_free(VC->Chrom);
	}

	return ret;
}


void variant_call_finit(PVARIANT_CALL VC)
{
	dym_array_finit_size_t(&VC->AltReads);
	dym_array_finit_size_t(&VC->RefReads);
	if (VC->Alt != NULL)
		utils_free(VC->Alt);
	
	if (VC->Ref != NULL)
		utils_free(VC->Ref);
	
	if (VC->ID != NULL)
		utils_free(VC->ID);
	
	if (VC->Chrom != NULL)
		utils_free(VC->Chrom);

	return;
}


boolean variant_call_equal(const VARIANT_CALL *VC1, const VARIANT_CALL *VC2)
{
	return (
		VC1->Pos == VC2->Pos &&
		strcasecmp(VC1->Chrom, VC2->Chrom) == 0 &&
		strcasecmp(VC1->Ref, VC2->Ref) == 0 &&
		strcasecmp(VC1->Alt, VC2->Alt) == 0
		);
}


ERR_VALUE vc_array_add(PGEN_ARRAY_VARIANT_CALL Array, const VARIANT_CALL *VC)
{
	VARIANT_CALL *tmp = Array->Data;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_NOT_FOUND;
	for (size_t i = 0; i < gen_array_size(Array); ++i) {
		if (variant_call_equal(tmp, VC)) {
			if (VC->AltWeight > tmp->AltWeight) {
				variant_call_finit(tmp);
				memcpy(tmp, VC, sizeof(VARIANT_CALL));
				ret = ERR_SUCCESS;
			}  else ret = ERR_ALREADY_EXISTS;
			break;
		}

		++tmp;
	}

	if (ret == ERR_NOT_FOUND)
		ret = dym_array_push_back_VARIANT_CALL(Array, *VC);

	return ret;
}


void vc_array_finit(PGEN_ARRAY_VARIANT_CALL Array)
{
	vc_array_clear(Array);
	dym_array_finit_VARIANT_CALL(Array);

	return;
}


void vc_array_clear(PGEN_ARRAY_VARIANT_CALL Array)
{
	PVARIANT_CALL tmp = Array->Data;

	for (size_t i = 0; i < gen_array_size(Array); ++i) {
		variant_call_finit(tmp);
		++tmp;
	}

	dym_array_clear_VARIANT_CALL(Array);

	return;
}


void vc_array_print(FILE *Stream, const char *ReferenceFile, const GEN_ARRAY_VARIANT_CALL *Array)
{
	const size_t variantCount = gen_array_size(Array);

	fprintf(Stream, "##fileformat=VCFv4.1\n");
	fprintf(Stream, "##fileDate=20160525\n");
	fprintf(Stream, "##source=GASSMV2\n");
	fprintf(Stream, "##reference=%s\n", ReferenceFile);
	{
		uint64_t len = 0;
		const char *chr = NULL;
		const VARIANT_CALL *vc = Array->Data;

		for (size_t i = 0; i < variantCount; ++i) {

			if (i == variantCount - 1 || (chr != NULL && strcasecmp(chr, vc->Chrom) != 0)) {
				if (len > 0)
					fprintf(Stream, "##contig=<ID=%s,length=%" PRIu64 ">\n", chr, len);

				chr = vc->Chrom;
				if (i == variantCount - 1)
					break;
			}
			
			chr = vc->Chrom;
			len = vc->Pos + strlen(vc->Ref);
			++vc;
		}
	}

	fprintf(Stream, "##phasing=partial\n");
	fprintf(Stream, "##INFO=<ID=RW,Number=1,Type=Integer,Description=\"Reference weight\">\n");
	fprintf(Stream, "##INFO=<ID=AW,Number=A,Type=Integer,Description=\"Allele weight\">\n");
	fprintf(Stream, "##INFO=<ID=RC,Number=1,Type=Integer,Description=\"Reference read count\">\n");
	fprintf(Stream, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alelle read count\">\n");
	fprintf(Stream, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(Stream, "##FORMAT=<ID=PS,Number=1,Type=String,Description=\"Phase number\">\n");
	fprintf(Stream, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t13350_1\n");

	const VARIANT_CALL *tmp = Array->Data;

	for (size_t i = 0; i < variantCount; ++i) {
		if (tmp->Valid) {
			const char *genotype = NULL;
			
			switch (tmp->PhaseType) {
				case vcptNone: genotype = "0/1"; break;
				case vcptOneTwo: genotype = "0|1"; break;
				case vcptTwoOne: genotype = "1|0"; break;
				case vcptBothAlt: genotype = "1|1"; break;
				default: assert(FALSE); break;
			}

			size_t refLen = strlen(tmp->Ref);
			size_t altLen = strlen(tmp->Alt);

			if (refLen == altLen) {
				for (size_t j = 0; j < refLen; ++j)
					fprintf(Stream, "%s\t%" PRIu64  "\t%s\t%c\t%c\t60\tPASS\tRW=%u;RC=%u;AW=%u;AC=%u\tGT:PS\t%s:%" PRIu64 "\n", tmp->Chrom, tmp->Pos + j, tmp->ID, tmp->Ref[j], tmp->Alt[j], (uint32_t)tmp->RefWeight, (uint32_t)gen_array_size(&tmp->RefReads), (uint32_t)tmp->AltWeight, (uint32_t)gen_array_size(&tmp->AltReads), genotype, tmp->PhasedPos);
			} else fprintf(Stream, "%s\t%" PRIu64  "\t%s\t%s\t%s\t60\tPASS\tRW=%u;RC=%u;AW=%u;AC=%u\tGT:PS\t%s:%" PRIu64 "\n", tmp->Chrom, tmp->Pos, tmp->ID, tmp->Ref, tmp->Alt, (uint32_t)tmp->RefWeight, (uint32_t)gen_array_size(&tmp->RefReads), (uint32_t)tmp->AltWeight, (uint32_t)gen_array_size(&tmp->AltReads), genotype, tmp->PhasedPos);
		}

		++tmp;
	}

	return;
}


static int _vc_comparator(const VARIANT_CALL *VC1, const VARIANT_CALL *VC2)
{
	int ret = strcasecmp(VC1->Chrom, VC2->Chrom);

	if (ret == 0) {
		if (VC1->Pos < VC2->Pos)
			return -1;

		if (VC1->Pos > VC2->Pos)
			return 1;

		int ret = strcmp(VC1->Ref, VC2->Ref);
		if (ret == 0)
			ret = strcmp(VC1->Alt, VC2->Alt);
	}

	return ret;
}

static int _size_t_comparer(const size_t *A, const size_t *B)
{
	if (*A < *B)
		return -1;

	if (*A > *B)
		return 1;

	return 0;
}


void vc_array_sort(PGEN_ARRAY_VARIANT_CALL Array)
{
	PVARIANT_CALL vc = Array->Data;

	qsort(Array->Data, gen_array_size(Array), sizeof(VARIANT_CALL), _vc_comparator);
	for (size_t i = 0; i < gen_array_size(Array); ++i) {
		qsort(vc->RefReads.Data, gen_array_size(&vc->RefReads), sizeof(size_t), _size_t_comparer);
		qsort(vc->AltReads.Data, gen_array_size(&vc->AltReads), sizeof(size_t), _size_t_comparer);
		++vc;
	}

	return;
}


ERR_VALUE vc_array_merge(PGEN_ARRAY_VARIANT_CALL Dest, PGEN_ARRAY_VARIANT_CALL Sources, const size_t SourceCount)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	size_t *indices = NULL;
	size_t *counts = NULL;
	size_t remainingCount = SourceCount;

	for (size_t i = 0; i < SourceCount; ++i) {
		vc_array_sort(Sources + i);
		if (gen_array_size(Sources + i) == 0)
			--remainingCount;
	}

	ret = ERR_SUCCESS;
	if (remainingCount > 0) {
		ret = utils_calloc(SourceCount, sizeof(size_t), &indices);
		if (ret == ERR_SUCCESS) {
			memset(indices, 0, SourceCount*sizeof(size_t));
			ret = utils_calloc(SourceCount, sizeof(size_t), &counts);
			if (ret == ERR_SUCCESS) {
				for (size_t i = 0; i < SourceCount; ++i)
					counts[i] = gen_array_size(Sources + i);

				size_t minSourceIndex = 0;
				uint64_t minValue = (uint64_t)-1;
				PVARIANT_CALL minVc = NULL;
				
				while (ret == ERR_SUCCESS && remainingCount > 0) {
					minValue = (uint64_t)-1;
					for (size_t i = 0; i < SourceCount; ++i) {
						if (indices[i] < counts[i]) {
							PVARIANT_CALL vc = dym_array_item_VARIANT_CALL(Sources + i, indices[i]);

							if (vc->Pos < minValue) {
								minValue = vc->Pos;
								minSourceIndex = i;
								minVc = vc;
							}
						}
					}

					if (minValue != (uint64_t)-1) {
						ret = vc_array_add(Dest, minVc);
						if (ret == ERR_ALREADY_EXISTS) {
							variant_call_finit(minVc);
							ret = ERR_SUCCESS;
						}

						indices[minSourceIndex]++;
						if (indices[minSourceIndex] == counts[minSourceIndex])
							--remainingCount;
					}
				}

				for (size_t i = 0; i < SourceCount; ++i)
					dym_array_clear_VARIANT_CALL(Sources + i);

				utils_free(counts);
			}

			utils_free(indices);
		}
	}

	return ret;
}


void vc_array_map_to_edges(PGEN_ARRAY_VARIANT_CALL VCArray)
{
	for (size_t i = 0; i < gen_array_size(VCArray); ++i) {
		PKMER_EDGE e = NULL;
		PVARIANT_CALL var = VCArray->Data + i;

		e = (PKMER_EDGE)var->Context;
		pointer_array_push_back_VARIANT_CALL(&e->VCs, var);
	}


	return;
}


ERR_VALUE vc_array_intersection(GEN_ARRAY_VARIANT_CALL *A1, const GEN_ARRAY_VARIANT_CALL *A2, PGEN_ARRAY_VARIANT_CALL Intersection)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const VARIANT_CALL *v1 = A1->Data;
	
	for (size_t i = 0; i < gen_array_size(A1); ++i) {
		const VARIANT_CALL *v2 = A2->Data;

		for (size_t j = 00; j < gen_array_size(A2); ++j) {			
			if (variant_call_equal(v1, v2)) {
				VARIANT_CALL tmp;

				ret = variant_call_init(v1->Chrom, v1->Pos, v1->ID, v1->Ref, strlen(v1->Ref), v1->Alt, strlen(v1->Alt), v1->Qual, &v1->RefReads, &v1->AltReads, &tmp);
				if (ret == ERR_SUCCESS) {
					tmp.RefWeight = v1->RefWeight;
					tmp.AltWeight = v1->AltWeight;
					if (vc_array_add(Intersection, &tmp) != ERR_SUCCESS)
						variant_call_finit(&tmp);
				}
			}

			++v2;
		}

		++v1;
	}

	return ret;
}
