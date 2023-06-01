/* C code document by David D. de Best for the master thesis 
   'Approximating the maximum weight matching problem with 
   positive-gain k-augmentations using the GraphBLAS standard', 
   supervised by prof. Rob H. Bisseling, Utrecht University, 2023 */


#include "generalheader.h"


// Calculates a set of positive-gain k-augmentations for a given k = 1,2,3, without conflicts and 
// stores the edges that need to be added to the matching in the augmentation matrix. 
// Also calculate diagnostics on runtime, number of found augmentations and conflicts.
// Permanently changes augmentation matrix.
// O(n+m) runtime for all k.
GrB_Info SearchKAugmentations (GrB_Matrix *Adjacency, GrB_Matrix *Matching, GrB_Matrix *Augmentation, 
                               GrB_Index k, LogEntry *logEntry, GrB_Index n)
{
    if (GIVEOUTPUT == 2)
        printf("\n------------------------------\nStarted %llu-augmentation search\n------------------------------\n", k);

    GrB_Info result;

    struct timespec runtimeStart, runtimeEnd;
    clock_gettime(CLOCK_REALTIME, &runtimeStart);

    // Perform the correct algorithm according to the given k.
    switch(k)
    {
        case 1: result = Search1AugsProcedure(Adjacency, Matching, Augmentation, logEntry, n); break;
        case 2: result = Search2AugsProcedure(Adjacency, Matching, Augmentation, logEntry, n); break;
        case 3: result = Search3AugsProcedure(Adjacency, Matching, Augmentation, logEntry, n); break;
    } 

    clock_gettime(CLOCK_REALTIME, &runtimeEnd);
    double runtime = TimeDifference(&runtimeStart, &runtimeEnd);

    (*logEntry).iterationType = 'S';
    (*logEntry).runtime = runtime;
    (*logEntry).kAugSearch = k;

    // Print search results.
    if (GIVEOUTPUT == 2)
    {
        printf("Searching %d-augmentations results:\n", (*logEntry).kAugSearch);

        if (result == GrB_NO_VALUE)
            printf(" - No positive-gain %d-augmentations found\n", (*logEntry).kAugSearch);
        else
        {
            printf(" - %d positive-gain %d-augmentations found\n", (*logEntry).nAugsFound, (*logEntry).kAugSearch);
            printf(" - %d %d-augmentations after removing conflicts\n", (*logEntry).nAugsNoConflicts, (*logEntry).kAugSearch);
        }

        printf("Runtime of the %d-augmentation search: %f seconds\n", (*logEntry).kAugSearch, (*logEntry).runtime);
        printf("\n----------------------------\nEnded %llu-augmentation search\n----------------------------\n", k);
    }

    return result;
}


// This method adds all edges that have an explicit value in the augmentation matrix to the matching,
// and removes the adjacent matched edges that cannot be in the matching anymore.
// Permanently changes matching and augmentation matrix.
// O(n) runtime.
GrB_Info FlipAugmentations (GrB_Matrix *Adjacency, GrB_Matrix *Matching, GrB_Matrix *Augmentation, 
                            double *CurrentMatchingWeight, LogEntry *logEntry, GrB_Index n)
{
    if (GIVEOUTPUT == 2)
        printf("\n-----------------------\nStarted flipping phase\n-----------------------\n");

    GrB_Vector m, a, r;
    GrB_Vector_new(&m, GrB_FP64, n);
    GrB_Vector_new(&a, GrB_FP64, n);
    GrB_Vector_new(&r, GrB_FP64, n);

    GrB_Matrix R;
    GrB_Matrix_new(&R, GrB_FP64, n, n);

    struct timespec runtimeStart, runtimeEnd;
    clock_gettime(CLOCK_REALTIME, &runtimeStart);


    // --- CHECKING DATA --- 

    if (DEBUGMODE)
    {
        GrB_Matrix T;
        GrB_Matrix_new(&T, GrB_FP64, n, n);

        GrB_Index t;
        // A valid augmentation cannot indicate to add any edges already in the matching, so
        // if the pair operator gives a positive result there must exist an error in a search algorithm.
        GrB_eWiseMult(T, NULL, NULL, GxB_PAIR_FP64, *Matching, *Augmentation, NULL);
        GrB_Matrix_nvals(&t, T);
        if (t > 0 && GIVEOUTPUT)
        {
            printf("Augmentation invalid: %llu edges in augmentation already in matching.\n Incorrect values indicated by matrix T:\n", t);
            GxB_print(T, GxB_SHORT);
        }

        GrB_Matrix_free(&T);
    }


    // --- START FLIPPING PHASE ALGORITHM ---  

    // We reduce the matching and augmentation matrices to vectors indicating which rows have explicit values.
    GrB_reduce(m, NULL, NULL, GxB_ANY_FP64, *Matching, NULL);
    GrB_reduce(a, NULL, NULL, GxB_ANY_FP64, *Augmentation, NULL);

    // We produce a vector indicating the intersection of m and a, containing ones at indices where both vectors have explicit values.
    GrB_eWiseMult(r, NULL, NULL, GxB_PAIR_FP64, m, a, NULL);

    // We create a matrix R with explicit values of values in the matching matrix that need to be removed.
    // We do this by selecting the correct rows from the matching indicated by vector r and adding its own transpose.
    // This is equivalent to selecting the correct columns as well and adding it to R, 
    // because of the structure and symmetry of the matching matrix.
    // Since the matching has O(n) values, obtaining the transpose is fast.
    SelectRowsInMatrixFromVector(Matching, &r, &R, n);
    GrB_eWiseAdd(R, NULL, NULL, GxB_ANY_FP64, R, R, GrB_DESC_T1);

    // Using an element-wise multiplication, we create the matrix containing the values 
    // that need to be added to the matching as indicated in the augmentation matrix.
    GrB_eWiseMult(*Augmentation, NULL, NULL, GrB_SECOND_FP64, *Augmentation, *Adjacency, NULL);

    // We apply the changes to the matching.
    GrB_eWiseAdd(*Matching, NULL, NULL, GrB_PLUS_FP64, *Matching, *Augmentation, NULL); 
    GrB_eWiseAdd(*Matching, NULL, NULL, GrB_MINUS_FP64, *Matching, R, NULL);
    GrB_select(*Matching, NULL, NULL, GrB_VALUENE_FP64, *Matching, 0, NULL);


    // --- END FLIPPING PHASE ALGORITHM ---

    // --- CHECKING DATA --- 

    if (DEBUGMODE)
    {
        GrB_Vector FullVector, T1, T2;
        GrB_Vector_new(&FullVector, GrB_INT64, n);
        GrB_Vector_new(&T1, GrB_INT64, n);
        GrB_Vector_new(&T2, GrB_INT64, n);

        // Full vector constant size in bytes: 224 | column vector
        GrB_assign(FullVector, NULL, NULL, 1, GrB_ALL, n, NULL);
    

        GrB_Index t1, t2;
        // A valid matching must have all pair-wise distinct edges, so each row and column can have at most one value.
        GrB_mxv(T1, NULL, NULL, GxB_PLUS_SECOND_INT64, *Matching, FullVector, NULL);
        GrB_vxm(T2, NULL, NULL, GxB_PLUS_FIRST_INT64, FullVector, *Matching, NULL);
        GrB_select(T1, NULL, NULL, GrB_VALUEGT_INT64, T1, 1, NULL);
        GrB_select(T2, NULL, NULL, GrB_VALUEGT_INT64, T2, 1, NULL);
        GrB_Vector_nvals(&t1, T1);
        GrB_Vector_nvals(&t2, T2);

        if ((t1 > 0 || t2 > 0) && GIVEOUTPUT)
        {
            printf("Matching invalid after flipping: %llu vertices have two adjacent edges.\n Vertices indicated by vector T1 and T2:\n", t1 + t2);
            GxB_print(T1, GxB_SHORT);
            GxB_print(T2, GxB_SHORT);
        }

        GrB_Vector_free(&FullVector);
        GrB_Vector_free(&T1);
        GrB_Vector_free(&T2);
    }

    // --- DIAGNOSTICS ---

    clock_gettime(CLOCK_REALTIME, &runtimeEnd);
    double runtime = TimeDifference(&runtimeStart, &runtimeEnd);

    // We report the actions performed by the flipping algorithm.
    GrB_Index nEdgesAdded, nEdgesRemoved;
    GrB_Matrix_nvals(&nEdgesAdded, *Augmentation); nEdgesAdded /= 2; // Divide by two because of symmetric matrix with two entries per edge.
    GrB_Matrix_nvals(&nEdgesRemoved, R); nEdgesRemoved /= 2;
    
    double oldWeight, weightAdded, weightRemoved, weightImprovement; 
    SumAllMatrixValues_FP64(Augmentation, &weightAdded); weightAdded /= 2;
    SumAllMatrixValues_FP64(&R, &weightRemoved); weightRemoved /= 2;
    weightImprovement = weightAdded - weightRemoved;

    (*logEntry).iterationType       = 'F';
    (*logEntry).runtime             = runtime;
    (*logEntry).weightAdded         = weightAdded;
    (*logEntry).weightRemoved       = weightRemoved;
    (*logEntry).weightImprovement   = weightImprovement;
    (*logEntry).oldWeight           = *CurrentMatchingWeight;
    (*logEntry).newWeight           = *CurrentMatchingWeight + weightAdded - weightRemoved;
    *CurrentMatchingWeight          = *CurrentMatchingWeight + weightAdded - weightRemoved;

    // Print flipping results.
    if (GIVEOUTPUT == 2)
    {
        printf("Flipping algorithm results: \n");
        printf(" - %llu edges added to the matching\n", nEdgesAdded);
        printf(" - %llu edges removed from the matching\n", nEdgesRemoved);
        printf("\n - Old matching weight: %f\n - New matching weight: %f\n - Added: %f   Removed: %f   Improvement: %f\n", 
                (*logEntry).oldWeight, (*logEntry).newWeight, 
                (*logEntry).weightAdded, (*logEntry).weightRemoved, 
                (*logEntry).weightImprovement);
        printf("Runtime of the flipping phase: %f seconds\n", runtime);
    }
    
    // We clear the augmentation matrix of all values so that it can be used in a new searching phase.
    GrB_Matrix_clear(*Augmentation);

    GrB_Vector_free(&m);
    GrB_Vector_free(&a);
    GrB_Vector_free(&r);

    GrB_Matrix_free(&R);

    if (GIVEOUTPUT == 2)
        printf("\n---------------------\nEnded flipping phase\n---------------------\n");

    return GrB_SUCCESS;
}
