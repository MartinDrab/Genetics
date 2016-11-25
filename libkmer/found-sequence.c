
#include "err.h"
#include "utils.h"
#include "refseq-storage.h"
#include "gen_dym_array.h"
#include "found-sequence.h"


static void _variant_free(PFOUND_SEQUENCE_VARIANT Variant)
{
	if (Variant->Seq1 != NULL)
		utils_free(Variant->Seq1);
	
	if (Variant->Seq2 != NULL)
		utils_free(Variant->Seq2);

	return;
}



ERR_VALUE found_sequence_init(const char *Sequence, const size_t Length, const size_t VariantCount, PFOUND_SEQUENCE FS)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(Length + 1, sizeof(char), &FS->Sequence);
	if (ret == ERR_SUCCESS) {
		memcpy(FS->Sequence, Sequence, Length*sizeof(char));
		FS->Sequence[Length] = '\0';
		FS->Len = Length;
		dym_array_init_FOUND_SEQUENCE_VARIANT(&FS->Variants, 140);
		dym_array_init_FOUND_SEQUENCE_VARIANT(&FS->ReadVariants, 140);
		ret = dym_array_reserve_FOUND_SEQUENCE_VARIANT(&FS->Variants, VariantCount);
		if (ret != ERR_SUCCESS) {
			dym_array_finit_FOUND_SEQUENCE_VARIANT(&FS->ReadVariants);
			dym_array_finit_FOUND_SEQUENCE_VARIANT(&FS->Variants);
			utils_free(FS->Sequence);
		}
	}

	return ret;
}


void found_sequence_finit(PFOUND_SEQUENCE FS)
{
	const size_t rCount = gen_array_size(&FS->ReadVariants);
	const size_t vCount = gen_array_size(&FS->Variants);

	for (size_t i = 0; i < rCount; ++i)
		_variant_free(dym_array_item_FOUND_SEQUENCE_VARIANT(&FS->ReadVariants, i));

	dym_array_finit_FOUND_SEQUENCE_VARIANT(&FS->ReadVariants);
	for (size_t i = 0; i < vCount; ++i)
		_variant_free(dym_array_item_FOUND_SEQUENCE_VARIANT(&FS->Variants, i));

	dym_array_finit_FOUND_SEQUENCE_VARIANT(&FS->Variants);
	utils_free(FS->Sequence);

	return;
}


ERR_VALUE found_sequence_alloc(const char *Sequence, const size_t Length, const size_t VariantCount, PFOUND_SEQUENCE *FS)
{
	PFOUND_SEQUENCE tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_malloc(sizeof(FOUND_SEQUENCE), &tmp);
	if (ret == ERR_SUCCESS) {
		ret = found_sequence_init(Sequence, Length, VariantCount, tmp);
		if (ret == ERR_SUCCESS)
			*FS = tmp;

		if (ret != ERR_SUCCESS)
			utils_free(tmp);
	}

	return ret;
}


void found_sequence_free(PFOUND_SEQUENCE FS)
{
	found_sequence_finit(FS);
	utils_free(FS);

	return;
}


ERR_VALUE found_sequence_set_variant(PFOUND_SEQUENCE FS, const size_t Index, PFOUND_SEQUENCE_VARIANT Variant)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	FOUND_SEQUENCE_VARIANT tmp;

	tmp.Seq1Len = Variant->Seq1Len;
	ret = utils_copy_string(Variant->Seq1, &tmp.Seq1);
	if (ret == ERR_SUCCESS) {
		tmp.Seq2Len = Variant->Seq2Len;
		ret = utils_copy_string(Variant->Seq2, &tmp.Seq2);
		if (ret == ERR_SUCCESS) {
			tmp.RefSeqStart = Variant->RefSeqStart;
			tmp.RefSeqEnd = Variant->RefSeqEnd;
			tmp.Seq1Type = Variant->Seq1Type;
			tmp.Seq2Type = Variant->Seq2Type;
			tmp.Seq1Weight = Variant->Seq1Weight;
			tmp.Seq2Weight = Variant->Seq2Weight;
			dym_array_push_back_no_alloc_FOUND_SEQUENCE_VARIANT(&FS->Variants, tmp);
		}

		if (ret != ERR_SUCCESS)
			utils_free(tmp.Seq1);
	}

	return ret;
}


boolean found_sequence_match(const FOUND_SEQUENCE *FS, const char *Seq, const size_t Length)
{
	const char *fPos = FS->Sequence;
	const char *sPos = Seq;
	boolean ret = TRUE;
	size_t remaining = Length;
	size_t cvIndex = 0;
	PFOUND_SEQUENCE_VARIANT currentVariant = NULL;

	if (gen_array_size(&FS->Variants) > 0) {
		currentVariant = dym_array_item_FOUND_SEQUENCE_VARIANT(&FS->Variants, 0);
		currentVariant->Reserved = 0;
	}

	while (ret && remaining > 0) {
		if (*fPos == '?') {
			ret = (currentVariant->Reserved < 2);
			if (ret) {
				const size_t index = currentVariant->Reserved;
				const char *varSeq = (index == 1) ? currentVariant->Seq2 : currentVariant->Seq1;
				const size_t varSeqlen = (index == 1) ? currentVariant->Seq2Len : currentVariant->Seq1Len;
			
				++currentVariant->Reserved;
				if (remaining >= varSeqlen && memcmp(varSeq, sPos, varSeqlen) == 0) {
					currentVariant->LastSPos = sPos;
					currentVariant->LastFPos = fPos;
					sPos += varSeqlen;
					remaining -= varSeqlen;
					++fPos;
					++currentVariant;
					++cvIndex;
					if (cvIndex < gen_array_size(&FS->Variants))
						currentVariant->Reserved = 0;
				}
			} else {
				ret = (cvIndex > 0);
				if (ret) {
					--currentVariant;
					--cvIndex;
					fPos = currentVariant->LastFPos;
					sPos = currentVariant->LastSPos;
					remaining = sPos - Seq;
				}
			}
		} else {
			ret = (*fPos == *sPos);
			if (ret) {
				++fPos;
				++sPos;
				--remaining;
			} else {
				ret = (cvIndex > 0);
				if (ret) {
					--currentVariant;
					--cvIndex;
					fPos = currentVariant->LastFPos;
					sPos = currentVariant->LastSPos;
					remaining = sPos - Seq;
				}
			}
		}
	}

	return ret;
}


ERR_VALUE found_sequence_build_read_variants(PFOUND_SEQUENCE FS, const POINTER_ARRAY_KMER_EDGE *PathEdges)
{
	uint32_t refseqStart = 0;
	uint32_t refseqEnd = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const KMER_EDGE *e = NULL;
	size_t startIndex = (size_t)-1;
	REFSEQ_STORAGE seqStorage;
	size_t weight = 0;
	size_t rsWeight = 0;

	ret = ERR_SUCCESS;
	rs_storage_init(&seqStorage);
	for (size_t i = 0; i < pointer_array_size(PathEdges); ++i) {
		e = *pointer_array_const_item_KMER_EDGE(PathEdges, i);
		if (e->Type == kmetRead) {
			if ((e->Source->Type == kmvtRefSeqMiddle || e->Source->Type == kmvtRefSeqStart)) {
				const KMER_EDGE *rsEdge = _get_refseq_or_variant_edge(e->Source);
				
				refseqStart = e->Source->RefSeqPosition + 1;
				startIndex = i;
				rs_storage_reset(&seqStorage);
				weight = read_info_weight(&e->ReadInfo);
				switch (rsEdge->Type) {
					case kmetReference:
						rsWeight = read_info_weight(&rsEdge->ReadInfo);
						break;
					case kmetVariant:
						rsWeight = rsEdge->Seq1Weight;
						break;
					default:
						__debugbreak();
						break;
				}
			}

			if (startIndex != (size_t)-1) {
				ret = rs_storage_add_edge(&seqStorage, e, TRUE);
				if (ret != ERR_SUCCESS)
					break;
			}

			if (e->Dest->Type == kmvtRefSeqMiddle) {
				FOUND_SEQUENCE_VARIANT var;
				
				if (!e->Dest->Helper && !e->Dest->Type != kmvtRefSeqEnd)
					rs_storage_remove(&seqStorage, 1);
				
				var.RefSeqStart = refseqStart;
				var.RefSeqEnd = e->Dest->RefSeqPosition;
				if (var.RefSeqStart < var.RefSeqEnd) {
					var.Seq1Type = kmetRead;
					ret = rs_storage_create_string(&seqStorage, &var.Seq1);
					if (ret == ERR_SUCCESS) {
						var.Seq1Len = strlen(var.Seq1);
						var.Seq1Weight = weight;
						var.Seq2Type = kmetNone;
						var.Seq2 = NULL;
						var.Seq2Len = 0;
						var.Seq2Weight = rsWeight;
						ret = dym_array_push_back_FOUND_SEQUENCE_VARIANT(&FS->ReadVariants, var);
						if (ret != ERR_SUCCESS)
							utils_free(var.Seq1);
					}

					if (ret != ERR_SUCCESS)
						break;
				}

				startIndex = (size_t)-1;
			}
		}
	}

	rs_storage_finit(&seqStorage);

	return ret;
}


ERR_VALUE variant_call_init(const char *Chrom, const uint64_t Pos, const char *ID, const char *Ref, size_t RefLen, const char *Alt, size_t AltLen, const uint8_t Qual, PVARIANT_CALL VC)
{
	char *tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

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
					/*
					while (RefLen > 1 && altLen > 1 && VC->Alt[altLen - 1] == VC->Ref[RefLen - 1]) {
						--RefLen;
						VC->Ref[RefLen] = '\0';
						--altLen;
						VC->Alt[altLen] = '\0';
					}
					*/
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
	const VARIANT_CALL *tmp = Array->Data;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	for (size_t i = 0; i < gen_array_size(Array); ++i) {
		if (variant_call_equal(tmp, VC)) {
			ret = ERR_ALREADY_EXISTS;
			break;
		}

		++tmp;
	}

	if (ret == ERR_SUCCESS)
		ret = dym_array_push_back_VARIANT_CALL(Array, *VC);

	return ret;
}


void vc_array_finit(PGEN_ARRAY_VARIANT_CALL Array)
{
	PVARIANT_CALL tmp = Array->Data;

	for (size_t i = 0; i < gen_array_size(Array); ++i) {
		variant_call_finit(tmp);
		++tmp;
	}

	dym_array_finit_VARIANT_CALL(Array);

	return;
}


void vc_array_print(FILE *Stream, const GEN_ARRAY_VARIANT_CALL *Array)
{
	const VARIANT_CALL *tmp = Array->Data;

	fprintf(Stream, "##fileformat=VCFv4.0\n");
	fprintf(Stream, "##fileDate=20160525\n");
	fprintf(Stream, "##source=GASSMV2\n");
	fprintf(Stream, "##reference=1\n");
	fprintf(Stream, "##phasing=partial\n");
	fprintf(Stream, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
	fprintf(Stream, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
	fprintf(Stream, "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">\n");
	fprintf(Stream, "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n");
	fprintf(Stream, "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">\n");
	fprintf(Stream, "##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">\n");
	fprintf(Stream, "##FILTER=<ID=q10,Description=\"Quality below 10\">\n");
	fprintf(Stream, "##FILTER=<ID=s50,Description=\"Less than 50 %% of samples have data\">\n");
	fprintf(Stream, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(Stream, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
	fprintf(Stream, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
	fprintf(Stream, "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">\n");
	fprintf(Stream, "#CHROM\tPOS\tID\tREF ALT\tQUAL\tFILTER\tINFO\n");
	for (size_t i = 0; i < gen_array_size(Array); ++i) {
		fprintf(Stream, "%s\t%I64u\t%s\t%s\t%s\t60\tPASS\t\"%Iu %Iu\"\n", tmp->Chrom, tmp->Pos, tmp->ID, tmp->Ref, tmp->Alt, tmp->RefWeight, tmp->AltWeight);
		++tmp;
	}

	return;
}


static int _vc_comparator(const VARIANT_CALL *VC1, const VARIANT_CALL *VC2)
{
	if (VC1->Pos < VC2->Pos)
		return -1;
	
	if (VC1->Pos > VC2->Pos)
		return 1;

	return 0;
}


void vc_array_sort(PGEN_ARRAY_VARIANT_CALL Array)
{
	qsort(Array->Data, gen_array_size(Array), sizeof(VARIANT_CALL), _vc_comparator);

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
				const VARIANT_CALL *minVc = NULL;
				
				while (ret == ERR_SUCCESS && remainingCount > 0) {
					minValue = (uint64_t)-1;
					for (size_t i = 0; i < SourceCount; ++i) {
						if (indices[i] < counts[i]) {
							const VARIANT_CALL *vc = dym_array_const_item_VARIANT_CALL(Sources + i, indices[i]);

							if (vc->Pos < minValue) {
								minValue = vc->Pos;
								minSourceIndex = i;
								minVc = vc;
							}
						}
					}

					if (minValue != (uint64_t)-1) {
						ret = vc_array_add(Dest, minVc);
						if (ret == ERR_ALREADY_EXISTS)
							ret = ERR_SUCCESS;

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
