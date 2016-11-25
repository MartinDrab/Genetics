
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <algorithm>
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
				default: assert(false); break;
			}
		}

		printf("\n");
	}

	return;
}


static void _print_score_matrix(const int *Matrix, const size_t Rows, const size_t Cols)
{
	for (size_t i = 0; i < Rows; ++i) {
		for (size_t j = 0; j < Cols; ++j)
			printf("%i ", item_2d(Matrix, Cols, i, j));

		printf("\n");
	}

	return;
}


static void _op_string_from_step_matrix(const EMatrixStep *StepMatrix, const size_t RowCount, const size_t ColumnCount, size_t MaxValueRow, size_t MaxValueCol, char **OperationString, size_t *OperationStringLen)
{
	char *opString = NULL;
	size_t opStringMax = MaxValueCol + MaxValueRow;

	opString = new char[opStringMax + 1];
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
				assert(false);
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

	return;
}


#define _compute_one_cell(aMatrix, aA, aB, aStepMatrix, aRowIndex, aColIndex, aColCount, aMatch, aMismatch, aIndel, aRowMaxes, aColMaxes, aResult)	\
{	\
	const bool matches = ((aB)[(aRowIndex)] == (aA)[(aColIndex)]);	\
	const int diag = std::max<int>(0, item_2d((aMatrix), (aColCount), (aRowIndex) - 1, (aColIndex) - 1) + ((matches) ? (aMatch) : (aMismatch)));			\
	const int left = std::max<int>(0, (aRowMaxes)[(aRowIndex)] + (aIndel));	\
	const int up = std::max<int>(0, (aColMaxes)[(aColIndex)] + (aIndel));		\
	\
	(aResult) = std::max<int>(0, std::max<int>(std::max<int>(up, left), diag));	\
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


void ssw_clever(const char *A, const size_t ALen, const char *B, const size_t BLen, const int Match, const int Mismatch, const int Indel, char **OperationString, size_t *OperationStringLen)
{
	int *matrix = NULL;
	EMatrixStep *steps = NULL;
	size_t rows = BLen + 1;
	size_t cols = ALen + 1;
	int maxValue = 0;
	size_t maxValueRow = 0;
	size_t maxValueCol = 0;
	int *rowMaxes = NULL;
	int *colMaxes = NULL;

	if (ALen == 0 || BLen == 0) {
		char *tmpOpString = NULL;
		size_t tmpOpStringLen = std::max<int>(ALen, BLen);

		tmpOpString = new char[tmpOpStringLen + 1];
		char zn = (ALen == 0) ? 'I' : 'D';
			
		tmpOpString[tmpOpStringLen] = '\0';
		for (size_t i = 0; i < tmpOpStringLen; ++i)
			tmpOpString[i] = zn;

		*OperationString = tmpOpString;
		*OperationStringLen = tmpOpStringLen;

		return;
	}

	matrix = new int[rows*cols];
	steps = new EMatrixStep[rows*cols];
	rowMaxes = new int[cols + rows];
	colMaxes = rowMaxes + rows;
	memset(rowMaxes, 0, (cols + rows)*sizeof(int));
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
			int newValue = 0;
						
			_compute_one_cell(matrix, A, B, steps, i, j, cols, Match, Mismatch, Indel, rowMaxes, colMaxes, newValue);
			_update_maximum_cell(maxValueRow, maxValueCol, maxValue, i, j, newValue);
		}
	}

//	_op_string_from_step_matrix(steps, rows, cols, maxValueRow, maxValueCol, OperationString, OperationStringLen);
	_op_string_from_step_matrix(steps, rows, cols, rows - 1, cols - 1, OperationString, OperationStringLen);
	delete[] rowMaxes;
	delete []steps;
	delete [] matrix;

	return;
}
