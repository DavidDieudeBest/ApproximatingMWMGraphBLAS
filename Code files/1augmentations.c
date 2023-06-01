/* C code document by David D. de Best for the master thesis 
   'Approximating the maximum weight matching problem with 
   positive-gain k-augmentations using the GraphBLAS standard', 
   supervised by prof. Rob H. Bisseling, Utrecht University, 2023 */


#include "generalheader.h"

// The procedure for calculating a set of positive-gain 1-augmentations.
// Permanently changes augmentation matrix.
// O(n+m) runtime.
GrB_Info Search1AugsProcedure (GrB_Matrix *Adjacency, GrB_Matrix *Matching, GrB_Matrix *Augmentation, 
                               LogEntry *logEntry, GrB_Index n)
{
    // --- Allocate space for data.
    GrB_Vector m_w;
    GrB_Vector_new(&m_w, GrB_FP64, n);

    GrB_Matrix C, AwoM, G1, D1;
    GrB_Matrix_new(&C, GrB_FP64, n, n);
    GrB_Matrix_new(&AwoM, GrB_FP64, n, n);
    GrB_Matrix_new(&G1, GrB_FP64, n, n);
    GrB_Matrix_new(&D1, GrB_FP64, n, n);

    GrB_Index g1_nvals;


    // --- Set up basic data.

    // We create vector m_w to contain the weights of the matched edges.
    // This is equal to the sum of of the rows of our Matching matrix, 
    // which contains at most one element per row.
    // Vector m_w needs to be full, so we set it up with all zeros.
    GrB_assign(m_w, NULL, NULL, 0, GrB_ALL, n, NULL);
    GrB_reduce(m_w, NULL, GrB_PLUS_FP64, GrB_PLUS_MONOID_FP64, *Matching, NULL);


    // --- Algorithm start.

    // We produce matrix [A \ M] containing all unmatched edges.
    ApplyMaskToMatrix(Adjacency, Matching, &AwoM, true);

    // For each potential edge that can be added to the matching, we subtract from 
    // its weight the weight of the already matched edge that needs to be removed,
    // and take the maximum for each vertex. We have thus found the maximum gain of 
    // part of a 1-augmentation for each vertex.
    MxMForcedDotMethod((GrB_Matrix)m_w, (GrB_Matrix)m_w, AwoM, GxB_PLUS_PLUS_FP64, &C);

    // We subtract for each vertex the weight of its adjacent edge that is matched,
    // resulting in a vector gain containing the maximum gain of 1-augmentations 
    // starting in each vertex.
    GrB_eWiseAdd(G1, NULL, NULL, GrB_MINUS_FP64, AwoM, C, NULL);
    GrB_select(G1, NULL, NULL, GrB_VALUEGT_FP64, G1, 0, NULL);
    GrB_Matrix_nvals(&g1_nvals, G1);

    // --- Early stop when no positive-gain augmentations exist.
    if (g1_nvals > 0)
    {   
        // We select the best 1-augmentation per vertex from G1 and save these in D1,
        // using the highest column index for ties in gain values.
        SelectHighestValuePerRow(&G1, &D1, n);

        // --- Solve conflicts.
        
        // Two types of conflicts are solved by multiplication of D1 with its own transpose.
        GrB_eWiseMult(*Augmentation, NULL, NULL, GxB_PAIR_FP64, D1, D1, GrB_DESC_T0);
        // We only select the explicit values indicating edges that need to be added to the matching.
        GrB_select(*Augmentation, NULL, NULL, GrB_VALUENE_FP64, *Augmentation, 0, NULL);
    }
    else
        return GrB_NO_VALUE;
    

    // --- Diagnostics.
    GrB_Index nPos1Augs, n1AugsNoConflicts;
    GrB_Matrix_nvals(&nPos1Augs, D1);
    GrB_Matrix_nvals(&n1AugsNoConflicts, *Augmentation); n1AugsNoConflicts /= 2; // Divide by two for augmentations occur twice in symmetric augmentation matrix.
    (*logEntry).nAugsFound = nPos1Augs;
    (*logEntry).nAugsNoConflicts = n1AugsNoConflicts;

    // --- Free data.
    GrB_Vector_free(&m_w);

    GrB_Matrix_free(&C);
    GrB_Matrix_free(&AwoM);
    GrB_Matrix_free(&G1);
    GrB_Matrix_free(&D1);

    return GrB_SUCCESS;
}
