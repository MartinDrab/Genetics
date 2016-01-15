
#include <assert.h>
#include <stdlib.h>
#include <emmintrin.h>
#include <smmintrin.h>
#include <stdio.h>
#include "err.h"
#include "utils.h"
#include "ssw.h"

/************************************************************************/
/*                             HELPER TYPES                             */
/************************************************************************/

typedef enum _EMatrixStep {
	msNone = 0,
	msDiagMatch,
	msDiagMisMatch,
	msUp,
	msLeft,
} EMatrixStep, *PEMatrixStep;

/************************************************************************/
/*                             HELPER FUNCTIONS                         */
/************************************************************************/

#define item_2d(aMatrix, aRowSize, aI, aJ)			(*((aMatrix) + (aRowSize)*(aI) + (aJ)))


static void _print_step_matrix(const EMatrixStep *Matrix, const size_t Rows, const size_t Cols)
{
	for (size_t i = 0; i < Rows; ++i) {
		for (size_t j = 0; j < Cols; ++j) {
			switch (item_2d(Matrix, Cols, i, j)) {
				case msNone: printf("N"); break;
				case msDiagMatch: 
				case msDiagMisMatch: printf("D"); break;
				case msUp: printf("U"); break;
				case msLeft: printf("L"); break;
				default: assert(FALSE); break;
			}
		}

		printf("\n");
	}

	return;
}


	
static void _print_score_matrix(const int32_t *Matrix, const size_t Rows, const size_t Cols)
{
	for (size_t i = 0; i < Rows; ++i) {
		for (size_t j = 0; j < Cols; ++j)
			printf("%i ", item_2d(Matrix, Cols, i, j));

		printf("\n");
	}

	return;
}


static ERR_VALUE _op_string_from_step_matrix(const EMatrixStep *StepMatrix, const size_t ColumnCount, size_t MaxValueRow, size_t MaxValueCol, char **OperationString, size_t *OperationStringLen)
{
	char *opString = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	size_t opStringMax = MaxValueCol + MaxValueRow;

	ret = utils_calloc(opStringMax + 1, sizeof(char), &opString);
	if (ret == ERR_SUCCESS) {
		size_t opStringIndex = opStringMax;

		opString[opStringMax] = '\0';
		while (item_2d(StepMatrix, ColumnCount, MaxValueRow, MaxValueCol) != msNone) {
			--opStringIndex;
			switch (item_2d(StepMatrix, ColumnCount, MaxValueRow, MaxValueCol)) {
				case msDiagMatch:
					opString[opStringIndex] = 'M';
					--MaxValueCol;
					--MaxValueRow;
					break;
				case msDiagMisMatch:
					opString[opStringIndex] = 'X';
					--MaxValueCol;
					--MaxValueRow;
					break;
				case msLeft:
					opString[opStringIndex] = 'I';
					--MaxValueCol;
					break;
				case msUp:
					opString[opStringIndex] = 'D';
					--MaxValueRow;
					break;
				default:
					assert(FALSE);
					break;
			}
		}

		while (MaxValueCol > 0) {
			--opStringIndex;
			opString[opStringIndex] = 'I';
			--MaxValueCol;
		}

		while (MaxValueRow > 0) {
			--opStringIndex;
			opString[opStringIndex] = 'D';
			--MaxValueRow;
		}

		memmove(opString, opString + opStringIndex, (opStringMax - opStringIndex + 1)*sizeof(char));
		*OperationString = opString;
		*OperationStringLen = opStringMax - opStringIndex;
		if (ret != ERR_SUCCESS)
			utils_free(opString);
	}

	return ret;
}

#define _compute_one_cell(aMatrix, aA, aB, aStepMatrix, aRowIndex, aColIndex, aColCount, aMatch, aMismatch, aIndel, aRowMaxes, aColMaxes, aResult)	\
{	\
	const boolean matches = ((aB)[(aRowIndex)] == (aA)[(aColIndex)]);	\
	const int32_t diag = item_2d((aMatrix), (aColCount), (aRowIndex) - 1, (aColIndex) - 1) + ((matches) ? (aMatch) : (aMismatch));			\
	const int32_t left = (aRowMaxes)[(aRowIndex)] + (aIndel);	\
	const int32_t up = (aColMaxes)[(aColIndex)] + (aIndel);		\
	\
	(aResult) = max(0, (max((max(up, left)), diag)));	\
	item_2d((aMatrix), (aColCount), (aRowIndex), (aColIndex)) = (aResult);	\
	if ((aResult) > up)	\
		(aColMaxes)[(aColIndex)] = (aResult);	\
																		\
	if ((aResult) > left)	\
		(aRowMaxes)[(aRowIndex)] = (aResult);	\
																		\
	if ((aResult) == diag)	\
		item_2d((aStepMatrix), (aColCount), (aRowIndex), (aColIndex)) = (matches) ? msDiagMatch : msDiagMisMatch;	\
	else if ((aResult) == left)	\
		item_2d((aStepMatrix), (aColCount), (aRowIndex), (aColIndex)) = msLeft;	\
	else if ((aResult) == up)	\
		item_2d((aStepMatrix), (aColCount), (aRowIndex), (aColIndex)) = msUp;	\
	else item_2d((aStepMatrix), (aColCount), (aRowIndex), (aColIndex)) = msNone;	\
	\
}	\

#define _update_maximum_cell(aMaxRow, aMaxCol, aMaxValue, aNewRow, aNewCol, aNewValue)	\
	{																					\
		if ((aNewValue) >= (aMaxValue)) {												\
			(aMaxValue) = (aNewValue);													\
			(aMaxRow) = (aNewRow);														\
			(aMaxCol) = (aNewCol);														\
		}																				\
	}																					\

static void _compute_one_stripe_sse32(int32_t *Matrix, const char *A, const char *B, EMatrixStep *StepMatrix, const size_t ColCount, const int32_t Match, const int32_t Mismatch, const int32_t Indel, int32_t *RowMaxes, int32_t *ColMaxes, size_t *MaxValueRow, size_t *MaxValueCol, int32_t *MaxValue)
{
	int32_t ret = 0;

	_compute_one_cell(Matrix, A, B, StepMatrix, 1, 1, ColCount, Match, Mismatch, Indel, RowMaxes, ColMaxes, ret);
	_update_maximum_cell(*MaxValueRow, *MaxValueCol, *MaxValue, 1, 1, ret);
	_compute_one_cell(Matrix, A, B, StepMatrix, 1, 2, ColCount, Match, Mismatch, Indel, RowMaxes, ColMaxes, ret);
	_update_maximum_cell(*MaxValueRow, *MaxValueCol, *MaxValue, 1, 2, ret);
	_compute_one_cell(Matrix, A, B, StepMatrix, 2, 1, ColCount, Match, Mismatch, Indel, RowMaxes, ColMaxes, ret);
	_update_maximum_cell(*MaxValueRow, *MaxValueCol, *MaxValue, 2, 1, ret);
	_compute_one_cell(Matrix, A, B, StepMatrix, 3, 1, ColCount, Match, Mismatch, Indel, RowMaxes, ColMaxes, ret);
	_update_maximum_cell(*MaxValueRow, *MaxValueCol, *MaxValue, 3, 1, ret);
	_compute_one_cell(Matrix, A, B, StepMatrix, 2, 2, ColCount, Match, Mismatch, Indel, RowMaxes, ColMaxes, ret);
	_update_maximum_cell(*MaxValueRow, *MaxValueCol, *MaxValue, 2, 2, ret);
	_compute_one_cell(Matrix, A, B, StepMatrix, 1, 3, ColCount, Match, Mismatch, Indel, RowMaxes, ColMaxes, ret);
	_update_maximum_cell(*MaxValueRow, *MaxValueCol, *MaxValue, 1, 3, ret);

	__declspec(align(16)) __m128i bv = _mm_setr_epi32(item_2d(Matrix, ColCount, 0, 3), item_2d(Matrix, ColCount, 1, 2), item_2d(Matrix, ColCount, 2, 1), item_2d(Matrix, ColCount, 3, 0));
	__declspec(align(16)) __m128i vA = _mm_setr_epi32(A[4], A[3], A[2], A[1]);
	__declspec(align(16)) __m128i vB = _mm_setr_epi32(B[1], B[2], B[3], B[4]);
	__declspec(align(16)) __m128i vMatch = _mm_setr_epi32(Match, Match, Match, Match);
	__declspec(align(16)) __m128i vMismatch = _mm_setr_epi32(Mismatch, Mismatch, Mismatch, Mismatch);
	__declspec(align(16)) __m128i vIndel = _mm_setr_epi32(Indel, Indel, Indel, Indel);
	__declspec(align(16)) __m128i vRowMaxes = _mm_setr_epi32(RowMaxes[1], RowMaxes[2], RowMaxes[3], RowMaxes[4]);
	__declspec(align(16)) __m128i vColMaxes = _mm_setr_epi32(ColMaxes[4], ColMaxes[3], ColMaxes[2], ColMaxes[1]);

	vRowMaxes = _mm_add_epi32(vRowMaxes, vIndel);
	vColMaxes = _mm_add_epi32(vColMaxes, vIndel);
	for (size_t i = 4; i < ColCount; ++i) {
		__declspec(align(16)) __m128i cmp;
		__declspec(align(16)) __m128i fv;
		__declspec(align(16)) int32_t tmp[4];
		__declspec(align(16)) int32_t up[4];
		__declspec(align(16)) int32_t diag[4];
		__declspec(align(16)) int32_t left[4];
		__declspec(align(16)) int32_t matches[4];

		// Compute new values for the score matrix
		cmp = _mm_cmpeq_epi32(vA, vB);
		_mm_store_si128((__m128i *)matches, cmp);
		fv = _mm_and_si128(cmp, vMatch);
		fv = _mm_add_epi32(fv, _mm_andnot_si128(cmp, vMismatch));
		fv = _mm_add_epi32(fv, bv);
		fv = _mm_max_epi32(fv, vRowMaxes);
		fv = _mm_max_epi32(fv, vColMaxes);

		// Store the values to the matrix
		_mm_store_si128((__m128i *)tmp, fv);
		item_2d(Matrix, ColCount, 1, i) = tmp[0];
		item_2d(Matrix, ColCount, 2, i - 1) = tmp[1];
		item_2d(Matrix, ColCount, 3, i - 2) = tmp[2];
		item_2d(Matrix, ColCount, 4, i - 3) = tmp[3];

		// Update the col and row of the maximum value in the matrix
		_update_maximum_cell(*MaxValueRow, *MaxValueCol, *MaxValue, 1, i, tmp[0]);
		_update_maximum_cell(*MaxValueRow, *MaxValueCol, *MaxValue, 2, i - 1, tmp[1]);
		_update_maximum_cell(*MaxValueRow, *MaxValueCol, *MaxValue, 3, i - 2, tmp[2]);
		_update_maximum_cell(*MaxValueRow, *MaxValueCol, *MaxValue, 4, i - 3, tmp[3]);
		
		// Compute the step matrix
		_mm_store_si128((__m128i *)diag, bv);
		_mm_store_si128((__m128i *)left, vRowMaxes);
		_mm_store_si128((__m128i *)up, vColMaxes);
		for (size_t k = 0; k < 4; ++k) {
			if (tmp[k] == diag[k])
				item_2d(StepMatrix, ColCount, k + 1, i - k) = (matches[k]) ? msDiagMatch : msDiagMisMatch;
			else if (tmp[k] == left[k])
				item_2d(StepMatrix, ColCount, k + 1, i - k) = msLeft;
			else if (tmp[k] == up[k])
				item_2d(StepMatrix, ColCount, k + 1, i - k) = msUp;
			else item_2d(StepMatrix, ColCount, k + 1, i - k) = msNone;
		}

		vColMaxes = _mm_shuffle_epi32(vColMaxes, _MM_SHUFFLE(0, 0, 1, 2));
		vColMaxes = _mm_insert_epi32(vColMaxes, ColMaxes[i + 1] + Indel, 0);
		// Move vA, vB, vb
		vA = _mm_shuffle_epi32(vA, _MM_SHUFFLE(0, 0, 1, 2));
		vA = _mm_insert_epi32(vA, A[i + 1], 0);
		vB = _mm_shuffle_epi32(vB, _MM_SHUFFLE(1, 2, 3, 3));
		vB = _mm_insert_epi32(vB, B[i + 1], 3);
		bv = fv;
		bv = _mm_shuffle_epi32(bv, _MM_SHUFFLE(0, 0, 1, 2));
		bv = _mm_insert_epi32(bv, item_2d(Matrix, ColCount, 0, i), 0);
	}
	
	_compute_one_cell(Matrix, A, B, StepMatrix, 2, ColCount - 1, ColCount, Match, Mismatch, Indel, RowMaxes, ColMaxes, ret);
	_update_maximum_cell(*MaxValueRow, *MaxValueCol, *MaxValue, 2, ColCount - 1, ret);
	_compute_one_cell(Matrix, A, B, StepMatrix, 3, ColCount - 2, ColCount, Match, Mismatch, Indel, RowMaxes, ColMaxes, ret);
	_update_maximum_cell(*MaxValueRow, *MaxValueCol, *MaxValue, 3, ColCount - 2, ret);
	_compute_one_cell(Matrix, A, B, StepMatrix, 4, ColCount - 3, ColCount, Match, Mismatch, Indel, RowMaxes, ColMaxes, ret);
	_update_maximum_cell(*MaxValueRow, *MaxValueCol, *MaxValue, 4, ColCount - 3, ret);
	_compute_one_cell(Matrix, A, B, StepMatrix, 3, ColCount - 1, ColCount, Match, Mismatch, Indel, RowMaxes, ColMaxes, ret);
	_update_maximum_cell(*MaxValueRow, *MaxValueCol, *MaxValue, 3, ColCount - 1, ret);
	_compute_one_cell(Matrix, A, B, StepMatrix, 4, ColCount - 2, ColCount, Match, Mismatch, Indel, RowMaxes, ColMaxes, ret);
	_update_maximum_cell(*MaxValueRow, *MaxValueCol, *MaxValue, 4, ColCount - 2, ret);
	_compute_one_cell(Matrix, A, B, StepMatrix, 4, ColCount - 1, ColCount, Match, Mismatch, Indel, RowMaxes, ColMaxes, ret);
	_update_maximum_cell(*MaxValueRow, *MaxValueCol, *MaxValue, 4, ColCount - 1, ret);

	return;
}



/************************************************************************/
/*                       PUBLIC FUNCTIONS                               */
/************************************************************************/


ERR_VALUE ssw_simple(const char *A, const size_t ALen, const char *B, const size_t BLen, const int Match, const int Mismatch, const int Indel, char **OperationString, size_t *OperationStringLen)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	int32_t *matrix = NULL;
	EMatrixStep *steps = NULL;
	size_t rows = BLen + 1;
	size_t cols = ALen + 1;
	int32_t maxValue = 0;
	size_t maxValueRow = 0;
	size_t maxValueCol = 0;

	ret = utils_calloc(rows*cols, sizeof(int32_t), (void **)&matrix);
	if (ret == ERR_SUCCESS) {
		ret = utils_calloc(rows*cols, sizeof(EMatrixStep), (void **)&steps);
		if (ret == ERR_SUCCESS) {
			for (size_t j = 0; j < cols; ++j) {
				item_2d(matrix, cols, 0, j) = 0;
				item_2d(steps, cols, 0, j) = msNone;
			}

			for (size_t i = 0; i < rows; ++i) {
				item_2d(matrix, cols, i, 0) = 0;
				item_2d(steps, cols, i, 0) = msNone;
			}

			--A;
			--B;
			for (size_t i = 1; i < rows; ++i) {
				for (size_t j = 1; j < cols; ++j) {
					int32_t left = 0;
					int32_t up = 0;
					const boolean matches = (B[i] == A[j]);
					const int32_t diag = item_2d(matrix, cols, i - 1, j - 1) + ((matches) ? Match : Mismatch);
					int32_t newValue = 0;

					for (size_t k = 0; k < i; ++k) {
						if (up <= item_2d(matrix, cols, k, j) + Indel)
							up = item_2d(matrix, cols, k, j) + Indel;
					}

					for (size_t k = 0; k < j; ++k) {
						if (left <= item_2d(matrix, cols, i, k) + Indel)
							left = item_2d(matrix, cols, i, k) + Indel;
					}

					newValue = max(0, (max((max(up, left)), diag)));
					item_2d(matrix, cols, i, j) = newValue;
					if (newValue == diag)
						item_2d(steps, cols, i, j) = (matches) ? msDiagMatch : msDiagMisMatch;
					else if (newValue == left)
						item_2d(steps, cols, i, j) = msLeft;
					else if (newValue == up)
						item_2d(steps, cols, i, j) = msUp;
					else item_2d(steps, cols, i, j) = msNone;

					if (newValue >= maxValue) {
						maxValue = newValue;
						maxValueRow = i;
						maxValueCol = j;
					}
				}
			}

			ret = _op_string_from_step_matrix(steps, cols, maxValueRow, maxValueCol, OperationString, OperationStringLen);
			utils_free(steps);
		}

		utils_free(matrix);
	}

	return ret;
}


ERR_VALUE ssw_clever(const char *A, const size_t ALen, const char *B, const size_t BLen, const int Match, const int Mismatch, const int Indel, char **OperationString, size_t *OperationStringLen)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	int32_t *matrix = NULL;
	EMatrixStep *steps = NULL;
	size_t rows = BLen + 1;
	size_t cols = ALen + 1;
	int32_t maxValue = 0;
	size_t maxValueRow = 0;
	size_t maxValueCol = 0;
	int32_t *rowMaxes = NULL;
	int32_t *colMaxes = NULL;

	ret = utils_calloc(rows*cols, sizeof(int32_t), (void **)&matrix);
	if (ret == ERR_SUCCESS) {
		ret = utils_calloc(rows*cols, sizeof(EMatrixStep), (void **)&steps);
		if (ret == ERR_SUCCESS) {
			ret = utils_calloc(cols + rows, sizeof(int32_t), &rowMaxes);
			if (ret == ERR_SUCCESS) {
				colMaxes = rowMaxes + rows;
				memset(rowMaxes, 0, (cols + rows)*sizeof(int32_t));
				for (size_t j = 0; j < cols; ++j) {
					item_2d(matrix, cols, 0, j) = 0;
					item_2d(steps, cols, 0, j) = msNone;
				}

				for (size_t i = 0; i < rows; ++i) {
					item_2d(matrix, cols, i, 0) = 0;
					item_2d(steps, cols, i, 0) = msNone;
				}

				--A;
				--B;
				for (size_t i = 1; i < rows; ++i) {
					for (size_t j = 1; j < cols; ++j) {
						int32_t newValue = 0;
						
						_compute_one_cell(matrix, A, B, steps, i, j, cols, Match, Mismatch, Indel, rowMaxes, colMaxes, newValue);
						_update_maximum_cell(maxValueRow, maxValueCol, maxValue, i, j, newValue);
					}
				}

				ret = _op_string_from_step_matrix(steps, cols, maxValueRow, maxValueCol, OperationString, OperationStringLen);
				utils_free(rowMaxes);
			}

			utils_free(steps);
		}

		utils_free(matrix);
	}

	return ret;
}


ERR_VALUE ssw_sse(const char *A, const size_t ALen, const char *B, const size_t BLen, const int Match, const int Mismatch, const int Indel, char **OperationString, size_t *OperationStringLen)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	int32_t *matrix = NULL;
	EMatrixStep *steps = NULL;
	size_t rows = BLen + 1;
	size_t cols = ALen + 1;
	int32_t maxValue = 0;
	size_t maxValueRow = 0;
	size_t maxValueCol = 0;
	int32_t *rowMaxes = NULL;
	int32_t *colMaxes = NULL;

	ret = utils_calloc(rows*cols, sizeof(int32_t), (void **)&matrix);
	if (ret == ERR_SUCCESS) {
		ret = utils_calloc(rows*cols, sizeof(EMatrixStep), (void **)&steps);
		if (ret == ERR_SUCCESS) {
			ret = utils_calloc(cols + rows, sizeof(int32_t), &rowMaxes);
			if (ret == ERR_SUCCESS) {
				colMaxes = rowMaxes + rows;
				memset(rowMaxes, 0, (cols + rows)*sizeof(int32_t));
				for (size_t j = 0; j < cols; ++j) {
					item_2d(matrix, cols, 0, j) = 0;
					item_2d(steps, cols, 0, j) = msNone;
				}

				for (size_t i = 0; i < rows; ++i) {
					item_2d(matrix, cols, i, 0) = 0;
					item_2d(steps, cols, i, 0) = msNone;
				}

				--A;
				--B;

				ret = _op_string_from_step_matrix(steps, cols, maxValueRow, maxValueCol, OperationString, OperationStringLen);
				utils_free(rowMaxes);
			}

			utils_free(steps);
		}

		utils_free(matrix);
	}

	return ret;
}


ERR_VALUE write_differences(const char *RefSeq, const size_t RegionStart, const size_t RegionLength, const char *OpString, const char *AltSeq)
{
	size_t pos = RegionStart;
	const char *rsPos = RefSeq + RegionStart;
	size_t remainingLength = RegionLength;
	boolean isDiff = FALSE;
	size_t diffStart = 0;
	char *rsContent = NULL;
	size_t rsIndex = 0;
	char *altContent = NULL;
	size_t altIndex = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(RegionLength + 1, sizeof(char), &rsContent);
	if (ret == ERR_SUCCESS) {
		ret = utils_calloc(strlen(AltSeq) + 1, sizeof(char), &altContent);
		if (ret == ERR_SUCCESS) {
			while (remainingLength > 0) {
				switch (*OpString) {
				case 'M':
					if (isDiff) {
						rsContent[rsIndex] = '\0';
						altContent[altIndex] = '\0';
						if (diffStart > RegionStart) {
							--diffStart;
							printf("1\t%u\t%c%s\t%c%s\n", (uint32_t)diffStart, RefSeq[diffStart], rsContent, RefSeq[diffStart], altContent);
						} else printf("1\t%u\t%s\t%s\n", (uint32_t)diffStart, rsContent, altContent);
						
						rsIndex = 0;
						altIndex = 0;
						diffStart = 0;
						isDiff = FALSE;
					}

					++rsPos;
					--remainingLength;
					++OpString;
					++AltSeq;
					++pos;
					break;
				case 'X':
					diffStart = pos;
					isDiff = TRUE;
					rsContent[rsIndex] = *rsPos;
					++rsIndex;
					altContent[altIndex] = *AltSeq;
					++altIndex;

					++rsPos;
					--remainingLength;
					++OpString;
					++AltSeq;
					++pos;
					break;
				case 'I':
					diffStart = pos;
					isDiff = TRUE;
					altContent[altIndex] = *AltSeq;
					++altIndex;

					++OpString;
					++AltSeq;
					break;
				case 'D':
					diffStart = pos;
					isDiff = TRUE;
					rsContent[rsIndex] = *rsPos;
					++rsIndex;

					++rsPos;
					--remainingLength;
					++OpString;
					++pos;
					break;
				case '\0':
					remainingLength = 0;
					break;
				default:
					assert(FALSE);
					break;
				}
			}

			if (isDiff) {
				rsContent[rsIndex] = '\0';
				altContent[altIndex] = '\0';
				if (diffStart > RegionStart) {
					--diffStart;
					printf("1\t%u\t%c%s\t%c%s\n", (uint32_t)diffStart, RefSeq[diffStart], rsContent, RefSeq[diffStart], altContent);
				} else printf("1\t%u\t%s\t%s\n", (uint32_t)diffStart, rsContent, altContent);

				rsIndex = 0;
				altIndex = 0;
				diffStart = 0;
				isDiff = FALSE;
			}

			utils_free(altContent);
		}

		utils_free(rsContent);
	}


	return ret;
}


void opstring_statistics(const char *OpString, const size_t OpStringLen, SSW_STATISTICS *Statistics)
{
	EGapType gt = gtMatch;

	memset(Statistics, 0, sizeof(SSW_STATISTICS));
	for (size_t i = 0; i < OpStringLen; ++i) {
		switch (*OpString) {
			case 'X':
				switch (gt) {
					case gtDeletion: ++Statistics->DeleteGapCount; break;
					case gtInsertion: ++Statistics->InsertGapCount; break;
				}

				++Statistics->TotalMismatches;
				gt = gtMismatch;
				break;
			case 'I':
				switch (gt) {
					case gtMismatch: ++Statistics->MismatchGapCount; break;
					case gtDeletion: ++Statistics->DeleteGapCount; break;
				}

				++Statistics->TotalInsertions;
				gt = gtInsertion;
				break;
			case 'D':
				switch (gt) {
					case gtMismatch: ++Statistics->MismatchGapCount; break;
					case gtInsertion: ++Statistics->InsertGapCount; break;
				}

				++Statistics->TotalDeletions;
				gt = gtDeletion;
				break;
			case 'M':
				switch (gt) {
					case gtMismatch: ++Statistics->MismatchGapCount; break;
					case gtInsertion: ++Statistics->InsertGapCount; break;
					case gtDeletion: ++Statistics->DeleteGapCount; break;
				}

				gt = gtMatch;
				break;
		}

		++OpString;
	}

	return;
}
