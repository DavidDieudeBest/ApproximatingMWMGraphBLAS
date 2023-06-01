/* C code document by David D. de Best for the master thesis 
   'Approximating the maximum weight matching problem with 
   positive-gain k-augmentations using the GraphBLAS standard', 
   supervised by prof. Rob H. Bisseling, Utrecht University, 2023 */


#include "generalheader.h"


// ----- CUSTOM TYPE CREATION/OPERATORS -----

// Create a new struct of two integers.
void MakeX1X2Struct (X1X2Struct *result, int64_t *x1, int64_t *x2)
{
    (*result).x1 = *x1;
    (*result).x2 = *x2;
};

// Create a new struct of two integers and a double.
void MakeVX1X2Struct (VX1X2Struct *result, double *value, X1X2Struct *x1x2)
{
    (*result).value = *value;
    (*result).x1 = (*x1x2).x1;
    (*result).x2 = (*x1x2).x2;
};

// Identity operator for custom type.
void IdentityMethod (VX1X2Struct *output, VX1X2Struct *input)
{
    *output = *input;
};

// Get value from custom type.
void GetValueFromTupleMethod (double *value, VX1X2Struct *tuple)
{
    *value = (*tuple).value;
};

// Get x1 from custom type.
void GetX1FromTupleMethod (int64_t *x1, VX1X2Struct *tuple)
{
    *x1 = (*tuple).x1;
};

// Get x2 from custom type.
void GetX2FromTupleMethod (int64_t *x2, VX1X2Struct *tuple)
{
    *x2 = (*tuple).x2;
};

// Add custom types in mixed fashion.
void PlusMixVX1X2Method (VX1X2Struct *result, VX1X2Struct *in1, VX1X2Struct *in2)
{
    (*result).value = (*in1).value + (*in2).value;
    (*result).x1 = (*in1).x1;
    (*result).x2 = (*in2).x2;
};

// Add all entries in custom type.
void PlusAllVX1X2Method (VX1X2Struct *result, VX1X2Struct *in1, VX1X2Struct *in2)
{
    (*result).value = (*in1).value + (*in2).value;
    (*result).x1 = (*in1).x1 + (*in2).x1;
    (*result).x2 = (*in2).x2 + (*in2).x2;
};

// Maximum operator for custom type based on value entry.
void MaxVX1X2Method (VX1X2Struct *result, VX1X2Struct *in1, VX1X2Struct *in2)
{
    if ((*in1).value > (*in2).value)
        *result = *in1;
    else
        *result = *in2;
};

// Second operator for custom type.
void SecondVX1X2Method (VX1X2Struct *result, VX1X2Struct *in1, VX1X2Struct *in2)
{
    *result = *in2;
};


// ----- CUSTOM TYPE MATRIX CREATION/TRANSFORMATION -----

// Split a matrix of custom type into three individual matrices.
// Permanently changes the value, x1 and x2 matrix.
// O(n+m) runtime.
GrB_Info SplitVX1X2Matrix (GrB_Matrix *VX1X2, GrB_Type ValueX1X2, GrB_Matrix *Value, GrB_Matrix *X1, GrB_Matrix *X2)
{
    GrB_UnaryOp GetValueFromTuple, GetX1FromTuple, GetX2FromTuple;
    GrB_UnaryOp_new(&GetValueFromTuple, (GxB_unary_function) GetValueFromTupleMethod, GrB_FP64, ValueX1X2);
    GrB_UnaryOp_new(&GetX1FromTuple, (GxB_unary_function) GetX1FromTupleMethod, GrB_INT64, ValueX1X2);
    GrB_UnaryOp_new(&GetX2FromTuple, (GxB_unary_function) GetX2FromTupleMethod, GrB_INT64, ValueX1X2);

    GrB_apply(*Value, NULL, NULL, GetValueFromTuple, *VX1X2, NULL);
    GrB_apply(*X1, NULL, NULL, GetX1FromTuple, *VX1X2, NULL);
    GrB_apply(*X2, NULL, NULL, GetX2FromTuple, *VX1X2, NULL);

    GrB_UnaryOp_free(&GetValueFromTuple);
    GrB_UnaryOp_free(&GetX1FromTuple);
    GrB_UnaryOp_free(&GetX2FromTuple);

    return GrB_SUCCESS;
}

// Multiply two vectors of length n of custom type to an nxn-matrix using a mask.
// Permanently changes the result matrix.
// O(n+m) runtime, by use of forced dot-product method.
GrB_Info MultiplyVX1X2VecsToMatr (GrB_Vector *V1, GrB_Vector *V2, GrB_Matrix *Mask, GrB_Type ValueX1X2, GrB_Matrix *Result)
{
    GrB_BinaryOp PlusMixVX1X2, PlusAllVX1X2;
    GrB_BinaryOp_new(&PlusMixVX1X2, (GxB_binary_function) PlusMixVX1X2Method, ValueX1X2, ValueX1X2, ValueX1X2);
    GrB_BinaryOp_new(&PlusAllVX1X2, (GxB_binary_function) PlusMixVX1X2Method, ValueX1X2, ValueX1X2, ValueX1X2);

    GrB_Scalar Identity;
    GrB_Scalar_new(&Identity, ValueX1X2);
    VX1X2Struct identity = {0, 0, 0};
    GrB_Scalar_setElement_UDT(Identity, &identity);

    GrB_Monoid PlusAllVX1X2_Monoid;
    GrB_Monoid_new_UDT(&PlusAllVX1X2_Monoid, PlusAllVX1X2, &Identity);

    GrB_Semiring MultiplyV1V2_Semiring;
    GrB_Semiring_new(&MultiplyV1V2_Semiring, PlusAllVX1X2_Monoid, PlusMixVX1X2);

    MxMForcedDotMethodAccum((GrB_Matrix)*V1, (GrB_Matrix)*V2, ((Mask == NULL) ? NULL : *Mask), MultiplyV1V2_Semiring, NULL, Result);
    
    GrB_Semiring_free(&MultiplyV1V2_Semiring);
    GrB_Monoid_free(&PlusAllVX1X2_Monoid);
    GrB_Scalar_free(&Identity);
    GrB_BinaryOp_free(&PlusMixVX1X2);
    GrB_BinaryOp_free(&PlusAllVX1X2);

    return GrB_SUCCESS;
}

// Construct a vector of custom type from three regular vectors.
// Permanently changes the result vector.
// O(n) runtime.
GrB_Info ConstructVX1X2Vector (GrB_Vector *Value, GrB_Vector *X1, GrB_Vector *X2, GrB_Vector *Result, GrB_Type ValueX1X2, GrB_Index n)
{
    GrB_Type X1X2;
    GrB_Type_new(&X1X2, sizeof(struct X1X2Struct));

    GrB_BinaryOp MakeX1X2, MakeVX1X2;
    GrB_BinaryOp_new(&MakeX1X2, (GxB_binary_function) MakeX1X2Struct, X1X2, GrB_INT64, GrB_INT64);
    GrB_BinaryOp_new(&MakeVX1X2, (GxB_binary_function) MakeVX1X2Struct, ValueX1X2, GrB_FP64, X1X2);

    GrB_Vector X1X2s;
    GrB_Vector_new(&X1X2s, X1X2, n);

    // If either element is not present, the result is implicit because of eWiseMult.
    GrB_eWiseMult(X1X2s, NULL, NULL, MakeX1X2, *X1, *X2, NULL);
    GrB_eWiseMult(*Result, NULL, NULL, MakeVX1X2, *Value, X1X2s, NULL);

    GrB_Vector_free(&X1X2s);

    GrB_BinaryOp_free(&MakeX1X2);
    GrB_BinaryOp_free(&MakeVX1X2);
    
    GrB_Type_free(&X1X2);

    return GrB_SUCCESS;
}

// Select all values in a matrix of custom type indicated with explicit values in a given mask, or its complement.
// If the result matrix is NULL, the original matrix is replaced.
// Permanently changes the input/result matrix.
// O(n+m) runtime.
GrB_Info ApplyMaskToMatrix_CustomType (GrB_Matrix *Input, GrB_Matrix *Mask, GrB_Matrix *Result, bool UseComplement, GrB_Type ValueX1X2)
{
    GrB_UnaryOp Identity;
    GrB_UnaryOp_new(&Identity, (GxB_unary_function) IdentityMethod, ValueX1X2, ValueX1X2);

    // We leave the input values unchanged by applying the identity operator, but using the structure of the mask.
    GrB_Matrix_apply(((Result == NULL) ? *Input : *Result),
              *Mask, NULL, Identity, *Input, (UseComplement ? GrB_DESC_RSC : GrB_DESC_RS));    

    GrB_UnaryOp_free(&Identity);

    return GrB_SUCCESS;
}