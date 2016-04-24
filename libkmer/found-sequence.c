
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
		if (ret == ERR_SUCCESS)
			dym_array_push_back_no_alloc_FOUND_SEQUENCE_VARIANT(&FS->Variants, tmp);

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
	PFOUND_SEQUENCE_VARIANT currentVariant = dym_array_item_FOUND_SEQUENCE_VARIANT(&FS->Variants, 0);
	size_t cvIndex = 0;

	if (gen_array_size(&FS->Variants) > 0)
		currentVariant->Reserved = 0;

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
