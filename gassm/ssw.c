
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
	msDiag,
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
				case msDiag: printf("D"); break;
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
	int maxValue = 0;
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
					int32_t diag = item_2d(matrix, cols, i - 1, j - 1) + ((B[i] == A[j]) ? Match : Mismatch);
					int32_t newValue = 0;

					for (size_t k = 0; k < i; ++k) {
						if (up < item_2d(matrix, cols, k, j) + Indel)
							up = item_2d(matrix, cols, k, j) + Indel;
					}

					for (size_t k = 0; k < j; ++k) {
						if (left < item_2d(matrix, cols, i, k) + Indel)
							left = item_2d(matrix, cols, i, k) + Indel;
					}

					newValue = max(0, (max((max(up, left)), diag)));
					item_2d(matrix, cols, i, j) = newValue;
					if (newValue == diag)
						item_2d(steps, cols, i, j) = msDiag;
					else if (newValue == left)
						item_2d(steps, cols, i, j) = msLeft;
					else if (newValue == up)
						item_2d(steps, cols, i, j) = msUp;
					else item_2d(steps, cols, i, j) = msNone;

					if (newValue > maxValue) {
						maxValue = newValue;
						maxValueRow = i;
						maxValueCol = j;
					}
				}
			}

//			printf("MAX: [%u:%u] = %i\n", maxValueRow, maxValueCol, maxValue);
//			_print_score_matrix(matrix, rows, cols);
//			printf("\n");
//			_print_step_matrix(steps, rows, cols);
			{
				char *opString = NULL;
				size_t opStringMax = maxValueCol + maxValueRow;

				ret = utils_calloc(opStringMax + 1, sizeof(char), &opString);
				if (ret == ERR_SUCCESS) {
					size_t opStringIndex = opStringMax;
					
					opString[opStringMax] = '\0';
					while (item_2d(steps, cols, maxValueRow, maxValueCol) != msNone) {
						--opStringIndex;
						switch (item_2d(steps, cols, maxValueRow, maxValueCol)) {
							case msDiag:
								opString[opStringIndex] = (A[maxValueCol] == B[maxValueRow]) ? 'M' : 'X';
								--maxValueCol;
								--maxValueRow;
								break;
							case msLeft:
								opString[opStringIndex] = 'I';
								--maxValueCol;
								break;
							case msUp:
								opString[opStringIndex] = 'D';
								--maxValueRow;
								break;
							default:
								assert(FALSE);
								break;
						}
					}

					while (maxValueCol > 0) {
						--opStringIndex;
						opString[opStringIndex] = 'I';
						--maxValueCol;
					}

					while (maxValueRow > 0) {
						--opStringIndex;
						opString[opStringIndex] = 'D';
						--maxValueRow;
					}

					memmove(opString, opString + opStringIndex, (opStringMax - opStringIndex + 1)*sizeof(char));
					*OperationString = opString;
					*OperationStringLen = opStringMax - opStringIndex;
					/*
					{
						uint32_t count = 1;
						char ch = opString[opStringIndex - 1];
						
						for (size_t i = opStringIndex - 1; i > 0; --i) {
							char tmp = opString[i - 1];

							if (tmp != ch) {
								printf("%u%c", count, ch);
								count = 1;
								ch = tmp;
							} else ++count;
						}

						if (count > 0)
							printf("%u%c\n", count, ch);
					}
					*/
				}
			}

			utils_free(steps);
		}

		utils_free(matrix);
	}

	return ret;
}
