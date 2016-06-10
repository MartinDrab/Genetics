
#include "err.h"
#include "utils.h"
#include "gen_dym_array.h"
#include "found-sequence.h"


static void _variant_free(PFOUND_SEQUENCE_VARIANT Variant)
{
	utils_free(Variant->Seq1);
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
		ret = dym_array_reserve_FOUND_SEQUENCE_VARIANT(&FS->Variants, VariantCount);
		if (ret != ERR_SUCCESS)
			utils_free(FS->Sequence);
	}

	return ret;
}


void found_sequence_finit(PFOUND_SEQUENCE FS)
{
	const size_t vCount = gen_array_size(&FS->Variants);

	for (size_t i = 0; i < vCount; ++i)
		_variant_free(dym_array_item_FOUND_SEQUENCE_VARIANT(&FS->Variants, i));

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


ERR_VALUE variant_call_init(const char *Chrom, const uint64_t Pos, const char *ID, const char *Ref, const char *Alt, const uint8_t Qual, PVARIANT_CALL VC)
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
			ret = utils_copy_string(Ref, &tmp);
			if (ret == ERR_SUCCESS) {
				VC->Ref = tmp;
				ret = utils_copy_string(Alt, &tmp);
				if (ret == ERR_SUCCESS)
					VC->Alt = tmp;
				
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
	utils_free(VC->Alt);
	utils_free(VC->Ref);
	utils_free(VC->ID);
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
		dym_array_push_back_VARIANT_CALL(Array, *VC);

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
	PVARIANT_CALL tmp = Array->Data;

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
	fprintf(Stream, "##FILTER=<ID=s50,Description=\"Less than 50 % of samples have data\">\n");
	fprintf(Stream, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(Stream, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
	fprintf(Stream, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
	fprintf(Stream, "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">\n");
	fprintf(Stream, "#CHROM\tPOS\tID\tREF ALT\tQUAL\tFILTER\tINFO\n");
	for (size_t i = 0; i < gen_array_size(Array); ++i) {
		fprintf(Stream, "%s\t%I64u\t%s\t%s\t%s\t60\tPASS\t\"%li-vs-%li\"\n", tmp->Chrom, tmp->Pos, tmp->ID, tmp->Ref, tmp->Alt, tmp->RefWeight, tmp->AltWeight);
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
