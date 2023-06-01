/* C code document by David D. de Best for the master thesis 
   'Approximating the maximum weight matching problem with 
   positive-gain k-augmentations using the GraphBLAS standard', 
   supervised by prof. Rob H. Bisseling, Utrecht University, 2023 */

#include "generalheader.h"


// ----- TIME -----

// Calculate the difference between two time values in seconds and milliseconds.
// O(1) runtime.
double TimeDifference (struct timespec *startTime, struct timespec *endTime)
{
    return (endTime->tv_sec - startTime->tv_sec) + 1e-9 * (endTime->tv_nsec - startTime->tv_nsec);
}


// ----- CALCULATION -----

// Sum all values in a floating point matrix.
// O(n+m) runtime.
GrB_Info SumAllMatrixValues_FP64 (GrB_Matrix *Input, double *Result)
{
    GrB_Index matrixSize;
    GrB_Matrix_nrows(&matrixSize, *Input);

    GrB_Vector a;
    GrB_Vector_new(&a, GrB_FP64, matrixSize);

    GrB_reduce(a, NULL, NULL, GrB_PLUS_MONOID_FP64, *Input, NULL);
    GrB_reduce(Result, NULL, GrB_PLUS_MONOID_FP64, a, NULL);

    GrB_Vector_free(&a);

    return GrB_SUCCESS;
}

// Matrix-matrix multiplication with the dot-product method forced and a mask for
// guaranteed O(n+m) runtime. Second input is transposed.
// Permanently changes result matrix.
// O(n+m) runtime.
GrB_Info MxMForcedDotMethod(GrB_Matrix A, GrB_Matrix B, GrB_Matrix Mask, GrB_Semiring Semiring, GrB_Matrix *Result)
{
    GrB_Descriptor desc;
    GrB_Descriptor_new(&desc);
    GrB_Descriptor_set(desc, GxB_AxB_METHOD, GxB_AxB_DOT);
    GrB_Descriptor_set(desc, GrB_MASK, GrB_STRUCTURE);
    GrB_Descriptor_set(desc, GrB_INP1, GrB_TRAN);

    GrB_mxm(*Result, Mask, NULL, Semiring, A, B, desc);

    GrB_Descriptor_free(&desc);

    return GrB_SUCCESS;
}

// Matrix-matrix multiplication with the dot-product method forced and a mask for
// guaranteed O(n+m) runtime. Second input is transposed. Uses accumulator operator.
// Permanently changes result matrix.
// O(n+m) runtime.
GrB_Info MxMForcedDotMethodAccum(GrB_Matrix A, GrB_Matrix B, GrB_Matrix Mask, GrB_Semiring Semiring, GrB_BinaryOp Accumulator, GrB_Matrix *Result)
{
    GrB_Descriptor desc;
    GrB_Descriptor_new(&desc);
    GrB_Descriptor_set(desc, GxB_AxB_METHOD, GxB_AxB_DOT);
    GrB_Descriptor_set(desc, GrB_MASK, GrB_STRUCTURE);
    GrB_Descriptor_set(desc, GrB_INP1, GrB_TRAN);

    GrB_mxm(*Result, Mask, Accumulator, Semiring, A, B, desc);

    GrB_Descriptor_free(&desc);

    return GrB_SUCCESS;
}

// Matrix-matrix multiplication with the dot-product method forced and a mask for
// guaranteed O(n+m) runtime. Second input is not transposed. Result is completely replaced.
// Permanently changes result matrix.
// O(n+m) runtime.
GrB_Info MxMForcedDotMethod2(GrB_Matrix A, GrB_Matrix B, GrB_Matrix Mask, GrB_Semiring Semiring, GrB_Matrix *Result)
{
    GrB_Descriptor desc;
    GrB_Descriptor_new(&desc);
    GrB_Descriptor_set(desc, GxB_AxB_METHOD, GxB_AxB_DOT);
    GrB_Descriptor_set(desc, GrB_MASK, GrB_STRUCTURE);
    GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);

    GrB_mxm(*Result, Mask, NULL, Semiring, A, B, desc);

    GrB_Descriptor_free(&desc);

    return GrB_SUCCESS;
}

// Calculate a basic upper bound for the MWM of an input graph.
// O(n+m) runtime.
GrB_Info MWMUpperBound(GrB_Matrix Input, double *Result)
{
    GrB_Index matrixSize;
    GrB_Matrix_nrows(&matrixSize, Input);

    GrB_Vector a;
    GrB_Vector_new(&a, GrB_FP64, matrixSize);

    GrB_reduce(a, NULL, NULL, GrB_MAX_MONOID_FP64, Input, NULL);
    GrB_reduce(Result, NULL, GrB_PLUS_MONOID_FP64, a, NULL);

    *Result /= 2;

    GrB_Vector_free(&a);

    return GrB_SUCCESS;
}


// ----- PERMUTATION -----

// Permute the rows of a matrix using a permutation matrix. No check is performed whether or not the used matrix is suitable 
// as a permutation matrix. Each explicit value in the permutation matrix is treated as a 1 implicitly, by use of the max.second semiring.
// If the result matrix is NULL, the original matrix is replaced.
// Permanently changes input/result matrix.
// O(n+m) runtime when the permutation matrix is valid.
GrB_Info PermuteMatrixRows (GrB_Matrix *Input, GrB_Matrix *PermutationMatrix, GrB_Matrix *Mask, GrB_Matrix *Result)
{
    return GrB_mxm(((Result == NULL) ? *Input : *Result), 
                   ((Mask == NULL) ? NULL : *Mask), 
                   NULL, GrB_MAX_SECOND_SEMIRING_FP64, *PermutationMatrix, *Input, GrB_DESC_RS);
}

// Permute the columns of a matrix using a permutation matrix. No check is performed whether or not the used matrix is suitable 
// as a permutation matrix. Each explicit value in the permutation matrix is treated as a 1 implicitly, by use of the max.first semiring.
// If the result matrix is NULL, the original matrix is replaced.
// Permanently changes input/result matrix.
// O(n+m) runtime when the permutation matrix is valid.
GrB_Info PermuteMatrixColumns (GrB_Matrix *Input, GrB_Matrix *PermutationMatrix, GrB_Matrix *Mask, GrB_Matrix *Result)
{
    return GrB_mxm(((Result == NULL) ? *Input : *Result), 
                   ((Mask == NULL) ? NULL : *Mask), 
                   NULL, GrB_MAX_FIRST_SEMIRING_FP64, *Input, *PermutationMatrix, GrB_DESC_RS);
}


// Permute a vector using a permutation matrix. The permutation is done via a premultiplication with the matrix. 
// No check is performed whether or not the used matrix is suitable as a permutation matrix. 
// Each explicit value in the permutation matrix is treated as a 1 implicitly, by use of the max.second semiring.
// If the result vector is NULL, the original vector is replaced.
// Permanently changes input/result vector.
// O(n) runtime when the permutation matrix is valid.
GrB_Info PermuteVector (GrB_Vector *Input, GrB_Matrix *PermutationMatrix, GrB_Vector *Mask, GrB_Vector *Result)
{
    return GrB_mxv(((Result == NULL) ? *Input : *Result), 
                   ((Mask == NULL) ? NULL : *Mask), 
                   NULL, GrB_MAX_SECOND_SEMIRING_FP64, *PermutationMatrix, *Input, GrB_DESC_RS);
}


// ----- SELECTION -----

// Select all values in a matrix indicated with explicit values in a given mask, or its complement.
// If the result matrix is NULL, the original matrix is replaced.
// Permanently changes the input/result matrix.
// O(n+m) runtime.
GrB_Info ApplyMaskToMatrix (GrB_Matrix *Input, GrB_Matrix *Mask, GrB_Matrix *Result, bool UseComplement)
{
    // We leave the input values unchanged by applying the identity operator, but using the structure of the mask.
    GrB_Matrix_apply(((Result == NULL) ? *Input : *Result),
              *Mask, NULL, GrB_IDENTITY_FP64, *Input, (UseComplement ? GrB_DESC_RSC : GrB_DESC_RS));        

    return GrB_SUCCESS;
}

// Select all values in a vector indicated with explicit values in a given mask, or its complement.
// If the result vector is NULL, the original vector is replaced.
// Permanently changes the input/result vector.
// O(n+m) runtime.
GrB_Info ApplyMaskToVector (GrB_Vector *Input, GrB_Vector *Mask, GrB_Vector *Result, bool UseComplement)
{
    // We leave the input values unchanged by applying the identity operator, but using the structure of the mask.
    GrB_Vector_apply(((Result == NULL) ? *Input : *Result),
              *Mask, NULL, GrB_IDENTITY_FP64, *Input, (UseComplement ? GrB_DESC_RSC : GrB_DESC_RS));               

    return GrB_SUCCESS;
}


// Selects the highest value in each row of the input matrix, using the highest column index in case of ties.
// If the result matrix is NULL, the original matrix is replaced.
// Permanently changes the input matrix.
// O(n+m) runtime.
GrB_Info SelectHighestValuePerRow (GrB_Matrix *Input, GrB_Matrix *Result, GrB_Index n)
{
    GrB_Vector rowMax, FullVector, maxIndex;
    GrB_Vector_new(&rowMax, GrB_FP64, n);
    GrB_Vector_new(&FullVector, GrB_FP64, n);
    GrB_Vector_new(&maxIndex, GrB_INT64, n);

    GrB_Matrix indices, maxIndexMatrix, rowMaxDiag, rowMaxMask, columnMaxMask;
    GrB_Matrix_new(&indices, GrB_INT64, n, n);
    GrB_Matrix_new(&maxIndexMatrix, GrB_INT64, n, n);
    GrB_Matrix_new(&rowMaxDiag, GrB_FP64, n, n);
    GrB_Matrix_new(&rowMaxMask, GrB_FP64, n, n);
    GrB_Matrix_new(&columnMaxMask, GrB_FP64, n, n);

    GrB_assign(FullVector, NULL, NULL, 1, GrB_ALL, n, NULL);

    GrB_Matrix *A = (Result == NULL) ? Input : Result;

    // First we save the maximum value in each row in a vector, and then diagonalize that vector.
    GrB_reduce(rowMax, NULL, NULL, GrB_MAX_MONOID_FP64, *Input, NULL);
    GrB_Matrix_diag(&rowMaxDiag, rowMax, 0);

    // We build a matrix with a value of 1 for each explicit value in the input matrix that is equal to the maximum in that row.
    // By applying this matrix as a mask to the input matrix, each explicit value in the input matrix 
    // that is equal to its row maximum is maintained, and all other values become explicit zeros.
    // O(n+m) runtime because of the mask and single value per row and column of the diagonal matrix.
    GrB_mxm(rowMaxMask, *Input, NULL, GxB_ANY_ISEQ_FP64, rowMaxDiag, *Input, GrB_DESC_S);
    GrB_select(rowMaxMask, NULL, NULL, GrB_VALUENE_FP64, rowMaxMask, 0, NULL);
    ApplyMaskToMatrix(Input, &rowMaxMask, A, false);

    // We save the highest column index of the remaining values in the input matrix for each row.
    GrB_mxv(maxIndex, NULL, NULL, GxB_MAX_SECONDI_INT64, *A, FullVector, NULL);
    // We build a matrix that sets for each explicit value in the input matrix the saved highest column index in that row.
    MxMForcedDotMethod((GrB_Matrix)maxIndex, (GrB_Matrix)FullVector, *A, GxB_MAX_FIRST_INT64, &maxIndexMatrix);

    // We build another matrix that has the column index for each explicit value in the input matrix.
    GrB_apply(indices, NULL, NULL, GrB_COLINDEX_INT64, *A, 0, NULL);
    // We only select the explicit values in the input matrix for which the actual column index coincides 
    // with the maximum column index value in the row.
    GrB_eWiseMult(columnMaxMask, *A, NULL, GrB_EQ_FP64, indices, maxIndexMatrix, GrB_DESC_S);   
    GrB_select(columnMaxMask, NULL, NULL, GrB_VALUENE_FP64, columnMaxMask, 0, NULL);
    ApplyMaskToMatrix(A, &columnMaxMask, NULL, false);

    GrB_Vector_free(&rowMax);
    GrB_Vector_free(&FullVector);
    GrB_Vector_free(&maxIndex);
    
    GrB_Matrix_free(&indices);
    GrB_Matrix_free(&maxIndexMatrix);
    GrB_Matrix_free(&rowMaxDiag);
    GrB_Matrix_free(&rowMaxMask);
    GrB_Matrix_free(&columnMaxMask);

    return GrB_SUCCESS;
}

// Selects the highest value in each column of the input matrix, using the highest row index in case of ties,
// by applying the row selection method to the transposed input matrix.
// If the result matrix is NULL, the original matrix is replaced.
// Permanently changes the input/result matrix.
// O(n+m) runtime.
GrB_Info SelectHighestValuePerColumn (GrB_Matrix *Input, GrB_Matrix *Result, GrB_Index n)
{
    GrB_Matrix *A = (Result == NULL) ? Input : Result;

    GrB_transpose(*A, NULL, NULL, *Input, NULL);
    GrB_Info result = SelectHighestValuePerRow(A, NULL, n);
    GrB_transpose(*A, NULL, NULL, *A, NULL);
    
    return result;
}

// Selection of rows in an input matrix by construction of a matrix that functions as a mask.
// If the result matrix is NULL, the original matrix is replaced.
// Permanently changes input/result matrix.
// O(n+m) runtime.
GrB_Info SelectRowsInMatrixFromVector (GrB_Matrix *Input, GrB_Vector *RowIndices, GrB_Matrix *Result, GrB_Index n)
{
    GrB_Matrix FullVector, Mask;
    GrB_Matrix_new(&FullVector, GrB_FP64, n, 1);
    GrB_Matrix_new(&Mask, GrB_FP64, n, n);
    GrB_assign(FullVector, NULL, NULL, 1, GrB_ALL, n, GrB_ALL, 1, NULL);

    // We produce a matrix with values 1 for each explicit value of the input matrix, 
    // but only in the rows indicated in the row indices vector.
    MxMForcedDotMethod((GrB_Matrix)*RowIndices, FullVector, *Input, GrB_MAX_SECOND_SEMIRING_FP64, &Mask);

    // Now we select only those values indicated by the mask in the input matrix.
    GrB_eWiseMult(((Result == NULL) ? *Input : *Result), NULL, NULL, GrB_FIRST_FP64, *Input, Mask, NULL);
    
    GrB_Matrix_free(&FullVector);
    GrB_Matrix_free(&Mask);

    return GrB_SUCCESS;
}

// Selection of columns in an input matrix by construction of a matrix that functions as a mask.
// If the result matrix is NULL, the original matrix is replaced.
// Permanently changes input/result matrix.
// O(n+m) runtime.
GrB_Info SelectColumnsInMatrixFromVector (GrB_Matrix *Input, GrB_Vector *ColumnIndices, GrB_Matrix *Result, GrB_Index n)
{
    GrB_Matrix FullVector, Mask;
    GrB_Matrix_new(&FullVector, GrB_FP64, n, 1);
    GrB_Matrix_new(&Mask, GrB_FP64, n, n);
    GrB_assign(FullVector, NULL, NULL, 1, GrB_ALL, n, GrB_ALL, 1, NULL);

    // We produce a matrix with values 1 for each explicit value of the input matrix, 
    // but only in the columns indicated in the row indices vector.
    MxMForcedDotMethod(FullVector, (GrB_Matrix)*ColumnIndices, *Input, GrB_MAX_FIRST_SEMIRING_FP64, &Mask);

    // Now we select only those values indicated by the mask in the input matrix.
    GrB_eWiseMult(((Result == NULL) ? *Input : *Result), NULL, NULL, GrB_FIRST_FP64, *Input, Mask, NULL);
    
    GrB_Matrix_free(&FullVector);
    GrB_Matrix_free(&Mask);

    return GrB_SUCCESS;
}


// ----- TRANSFORMATION -----

// Round off all values in an input matrix to a given decimal. Usable for throwing away round-off errors.
// Consider removing explicit zeroes afterwards, this is not done in the method.
// If the result matrix is NULL, the original matrix is replaced.
// Permanently changes input/result matrix.
// O(n+m) runtime.
GrB_Info RoundOffAllValuesInMatrix (GrB_Matrix *Input, int ToDecimal, GrB_Matrix *Result)
{
    GrB_Index power = pow(10, ToDecimal);
    GrB_Matrix *A = (Result == NULL) ? Input : Result;

    // We multiply all values by a given power of 10, round off to the nearest integer and divide by the power of 10.
    GrB_Matrix_apply_BinaryOp1st_FP64(*A, NULL, NULL, GrB_TIMES_FP64, power, *Input, NULL);
    GrB_apply(*A, NULL, NULL, GxB_ROUND_FP64, *A, NULL);
    GrB_Matrix_apply_BinaryOp1st_FP64(*A, NULL, NULL, GxB_RDIV_FP64, power, *A, NULL);

    return GrB_SUCCESS;
}

// Round off all values in an input vector to a given decimal. Usable for throwing away round-off errors.
// Consider removing explicit zeroes afterwards, this is not done in the method.
// If the result vector is NULL, the original matrix is replaced.
// Permanently changes input/result vector.
// O(n+m) runtime.
GrB_Info RoundOffAllValuesInVector (GrB_Vector *Input, int ToDecimal, GrB_Vector *Result)
{
    GrB_Index power = pow(10, ToDecimal);
    GrB_Vector *A = (Result == NULL) ? Input : Result;

    // We multiply all values by a given power of 10, round off to the nearest integer and divide by the power of 10.
    GrB_Vector_apply_BinaryOp1st_FP64(*A, NULL, NULL, GrB_TIMES_FP64, power, *Input, NULL);
    GrB_apply(*A, NULL, NULL, GxB_ROUND_FP64, *A, NULL);
    GrB_Vector_apply_BinaryOp1st_FP64(*A, NULL, NULL, GxB_RDIV_FP64, power, *A, NULL);

    return GrB_SUCCESS;
}

// Create a matrix with a 1 value per row at a given column index from the input vector,
// so if Input[i] = x then Result[i,x] = 1. 
// Permanently changes result matrix.
// O(n) runtime.
GrB_Info CreateMatrixFromColumnIndices (GrB_Vector *Input, GrB_Matrix *Result)
{
    GrB_Index *columnIndices, *rowIndices, nvals;
    GrB_Vector_nvals(&nvals, *Input);
    columnIndices = malloc(nvals * sizeof(GrB_Index));
    rowIndices = malloc(nvals * sizeof(GrB_Index));

    int FullVector[nvals];
    for (size_t i = 0; i < nvals; i++)
        FullVector[i] = 1;

    GrB_Vector_extractTuples(rowIndices, columnIndices, &nvals, *Input);
    GrB_Matrix_build(*Result, rowIndices, columnIndices, FullVector, nvals, GrB_PLUS_INT64);

    free(columnIndices);
    free(rowIndices);
    
    return GrB_SUCCESS;
}

// Sets each value in the input matrix on the diagonal as an implicit value in the result
// matrix and copies over all other explicit values.
// Permanently changes result matrix.
// O(n+m) runtime.
GrB_Info RemoveDiagonalFromMatrix(GrB_Matrix *Input, GrB_Matrix *Result)
{
    GrB_Index matrixSize;
    GrB_Matrix_nrows(&matrixSize, *Input);

    GrB_Matrix Mask;
    GrB_Matrix_new(&Mask, GrB_INT64, matrixSize, matrixSize);

    GrB_select(Mask, NULL, NULL, GrB_DIAG, *Input, 0, NULL);
    ApplyMaskToMatrix(Input, &Mask, Result, true);

    GrB_Matrix_free(&Mask);

    return GrB_SUCCESS;
}