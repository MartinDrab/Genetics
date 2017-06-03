
#include <assert.h>
#include <stdlib.h>
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
					opString[opStringIndex] = 'D';
					--MaxValueCol;
					break;
				case msUp:
					opString[opStringIndex] = 'I';
					--MaxValueRow;
					break;
				default:
					assert(FALSE);
					break;
			}
		}

		while (MaxValueCol > 0) {
			--opStringIndex;
			opString[opStringIndex] = 'D';
			--MaxValueCol;
		}

		while (MaxValueRow > 0) {
			--opStringIndex;
			opString[opStringIndex] = 'I';
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
	const int32_t diag = max(0, item_2d((aMatrix), (aColCount), (aRowIndex) - 1, (aColIndex) - 1) + ((matches) ? (aMatch) : (aMismatch)));			\
	const int32_t left = max(0, (aRowMaxes)[(aRowIndex)] + (aIndel));	\
	const int32_t up = max(0, (aColMaxes)[(aColIndex)] + (aIndel));		\
	\
	(aResult) = max(0, max(max(up, left), diag));	\
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

			ret = _op_string_from_step_matrix(steps, cols, rows - 1, cols - 1, OperationString, OperationStringLen);
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

	if (ALen == 0 || BLen == 0) {
		char *tmpOpString = NULL;
		size_t tmpOpStringLen = max(ALen, BLen);

		ret = utils_calloc(tmpOpStringLen + 1, sizeof(char), &tmpOpString);
		if (ret == ERR_SUCCESS) {
			char zn = (ALen == 0) ? 'I' : 'D';
			
			tmpOpString[tmpOpStringLen] = '\0';
			for (size_t i = 0; i < tmpOpStringLen; ++i)
				tmpOpString[i] = zn;

			*OperationString = tmpOpString;
			*OperationStringLen = tmpOpStringLen;
		}

		return ret;
	}

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
