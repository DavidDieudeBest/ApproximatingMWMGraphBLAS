/* C code document by David D. de Best for the master thesis 
   'Approximating the maximum weight matching problem with 
   positive-gain k-augmentations using the GraphBLAS standard', 
   supervised by prof. Rob H. Bisseling, Utrecht University, 2023 */


#include "generalheader.h"

// Finds the best path 2-augmentation for each potential center i.e. each matched edge (i,j).
// Its gain value is stored in matrix G2p at (i,j) and explicit values are stored in D2p at 
// (i,k) and (j,l) to indicate the used 1-augmentations.
// Permanently changes G2p, D2p and the augmentation matrix.
// O(n+m) runtime.
GrB_Info SearchPath2Augs (GrB_Matrix *Adjacency, GrB_Matrix *Matching, GrB_Matrix *Augmentation, 
                              GrB_Matrix *G2p, GrB_Matrix *D2p, GrB_Index n)
{
    // --- Allocate space for data.
    GrB_Matrix M_u, B, C, E, F, AwoM, G1, D1;
    GrB_Matrix_new(&M_u, GrB_FP64, n, n);
    GrB_Matrix_new(&B, GrB_FP64, n, n);
    GrB_Matrix_new(&C, GrB_FP64, n, n);
    GrB_Matrix_new(&E, GrB_FP64, n, n);
    GrB_Matrix_new(&F, GrB_FP64, n, n);
    GrB_Matrix_new(&AwoM, GrB_FP64, n, n);
    GrB_Matrix_new(&G1, GrB_FP64, n, n);
    GrB_Matrix_new(&D1, GrB_FP64, n, n);

    GrB_Vector FullVector, m, m_w, m_u, a, a_g, a_gi, a_gj, v, c, c_u, c_l, g2p;
    GrB_Vector_new(&FullVector, GrB_FP64, n);
    GrB_Vector_new(&m, GrB_FP64, n);
    GrB_Vector_new(&m_w, GrB_FP64, n);
    GrB_Vector_new(&m_u, GrB_FP64, n);
    GrB_Vector_new(&a, GrB_FP64, n);
    GrB_Vector_new(&a_g, GrB_FP64, n);
    GrB_Vector_new(&a_gi, GrB_FP64, n);
    GrB_Vector_new(&a_gj, GrB_FP64, n);
    GrB_Vector_new(&v, GrB_FP64, n);
    GrB_Vector_new(&c, GrB_FP64, n);
    GrB_Vector_new(&c_u, GrB_FP64, n);
    GrB_Vector_new(&c_l, GrB_FP64, n);
    GrB_Vector_new(&g2p, GrB_FP64, n);

    GrB_Index g2p_nvals;

    // --- Set up basic data.

    // Full vector constant size in bytes: 224 | Column vector
    GrB_assign(FullVector, NULL, NULL, 1, GrB_ALL, n, NULL);

    // Create vector m with all indices of the match of vertices
    GrB_mxv(m, NULL, NULL, GxB_PLUS_SECONDI_INT64, *Matching, FullVector, NULL);

    // Create vector m_w with the weights of all matched edges.
    GrB_assign(m_w, NULL, NULL, 0, GrB_ALL, n, NULL);
    GrB_reduce(m_w, NULL, GrB_PLUS_FP64, GrB_PLUS_MONOID_FP64, *Matching, NULL);

    // Create a vector m_u with indicating all vertices i for matched edges (i,j) with i < j.
    // Thus having a vector indicating exactly one of each pair of matched vertices.
    GrB_select(M_u, NULL, NULL, GrB_TRIU, *Matching, 0, NULL);
    GrB_reduce(m_u, NULL, NULL, GrB_PLUS_MONOID_FP64, M_u, NULL);


    // --- Algorithm start.

    // An adapted version of the 1-augmentation algorithm is performed to obtain G1 and D1 matrices.
    ApplyMaskToMatrix(Adjacency, Matching, &AwoM, true);

    MxMForcedDotMethod((GrB_Matrix)m_w, (GrB_Matrix)m_w, AwoM, GxB_PLUS_PLUS_FP64, &C);
    GrB_eWiseAdd(G1, NULL, NULL, GrB_MINUS_FP64, AwoM, C, NULL);
    GrB_Matrix_clear(C);

    SelectHighestValuePerRow(&G1, &D1, n);

    // We want to detect invalid 2-augmentations (i,j,k,l) where k=l and replace them.
    // First we only need the 1-augmentations for matched vertices.
    SelectRowsInMatrixFromVector(&D1, &m, D2p, n);

    // Save the potential k and l vertices
    GrB_mxv(a, NULL, NULL, GxB_PLUS_SECONDI_INT64, *D2p, FullVector, NULL); 
    // Save the gain of the potential 1-augmentations used as k and l vertices
    GrB_reduce(a_g, NULL, NULL, GxB_ANY_FP64, *D2p, NULL);

    // Calculate the combined gain of the best 1-augmentations on each side of matched edges.
    // Those are all the centers. An explicit value falls away when there are not two 1-augmentations available,
    // and thus a path 2-augmentation is not possible.
    MxMForcedDotMethod((GrB_Matrix)a_g, (GrB_Matrix)a_g, *Matching, GxB_PLUS_PLUS_FP64, G2p);

    // We add the weight of the center matched edge for each 2-augmentation to obtain its true gain.
    GrB_eWiseMult(*G2p, NULL, NULL, GrB_PLUS_FP64, *G2p, *Matching, NULL);

    // We select only the 2-augmentations that have a positive gain value. Others are not of interest since
    // changing either 1-augmentation will not improve the gain.
    GrB_select(*G2p, NULL, NULL, GrB_VALUEGT_FP64, *G2p, 0, NULL);

    GrB_Matrix_nvals(&g2p_nvals, *G2p);

    // --- Early stop when no positive-gain augmentations exist.
    if (g2p_nvals > 0)
    {  
        // Keep only the vertices k and l in vector a where we can make a positive-gain 2-augmentation.
        GrB_reduce(v, NULL, NULL, GrB_PLUS_FP64, *G2p, NULL);
        ApplyMaskToVector(&a, &v, NULL, false);
        GrB_Vector_clear(v);

        // We attempt to change both 1-augmentations of invalid 2-augmentations separately and choose the best option.
        // Save in v the permuted k and l vertices of each 2-augmentation
        PermuteVector(&a, Matching, NULL, &v);
        // We find the vertices where k=l occurs when performing a eq v.
        // We save the i vertices in c_u and j vertices in c_l
        GrB_eWiseMult(c_l, NULL, NULL, GrB_EQ_FP64, a, v, NULL);

        GrB_select(c_l, NULL, NULL, GrB_VALUENE_FP64, c_l, 0, NULL);
        ApplyMaskToVector(&c_l, &m_u, &c_u, false);
        ApplyMaskToVector(&c_l, &m_u, &c_l, true);
        

        // We remove the best 1-augmentation for the i vertex for each matched edge by storing them in B
        // and using its complement as a mask for G1, which contains all possible 1-augmentations.
        // Then we rerun the 1-augmentation algorithm and save the result in E and the obtained gain for each
        // new best 1-augmentation in a_gi.
        SelectRowsInMatrixFromVector(&D1, &c_u, &B, n);
        ApplyMaskToMatrix(&G1, &B, NULL, true);
        GrB_Matrix_clear(D1);
        SelectHighestValuePerRow(&G1, &D1, n);
        SelectRowsInMatrixFromVector(&D1, &m, &E, n);
        GrB_reduce(a_gi, NULL, NULL, GrB_PLUS_FP64, E, NULL);

        // We add the removed 1-augmentations back to G1.
        GrB_eWiseAdd(G1, NULL, NULL, GrB_PLUS_FP64, G1, B, NULL);
        GrB_Matrix_clear(B);
        GrB_Matrix_clear(C);

        // We perform the same steps by now removing the best 1-augmentations for the j vertex of matched edges,
        // and running the 1-augmentation algorithm, saving results in C and a_gj.
        SelectRowsInMatrixFromVector(&D1, &c_l, &B, n);
        ApplyMaskToMatrix(&G1, &B, NULL, true);
        GrB_Matrix_clear(D1);
        SelectHighestValuePerRow(&G1, &D1, n);
        SelectRowsInMatrixFromVector(&D1, &m, &C, n);
        GrB_reduce(a_gj, NULL, NULL, GrB_PLUS_FP64, C, NULL);

        // We recalculate the combined gain of the best 1-augmentations for each side of matched edges,
        // now using the results from the previous steps stored in a_gi and a_gj. Only needed at locations
        // where a positive-gain 2-augmentation was possible in the first place, stored in G2p.
        GrB_Matrix_clear(B);
        MxMForcedDotMethod((GrB_Matrix)a_gj, (GrB_Matrix)a_gj, *G2p, GxB_PLUS_PLUS_FP64, &B);
        MxMForcedDotMethod((GrB_Matrix)a_gi, (GrB_Matrix)a_gi, *G2p, GxB_PLUS_PLUS_FP64, G2p);

        // We take the best 2-augmentation found at each location from G2p and B.
        GrB_assign(F, B, NULL, 1, GrB_ALL, n, GrB_ALL, n, NULL);
        ApplyMaskToMatrix(&F, G2p, NULL, true);
        GrB_eWiseMult(F, NULL, GrB_PLUS_FP64, GrB_GT_FP64, B, *G2p, NULL);
        GrB_select(F, NULL, NULL, GrB_VALUEEQ_FP64, F, 1, NULL);

        GrB_eWiseAdd(*G2p, F, NULL, GrB_FIRST_FP64, B, *G2p, NULL); 

        // We save in B the locations where the 1-augmentations from a_gj gave the best 2-augmentation, and in v the vertices.
        GrB_Vector_clear(v);
        GrB_reduce(v, NULL, NULL, GrB_PLUS_FP64, F, NULL);

        // We add the weight of the center matched edge to each 2-augmentation.
        GrB_eWiseMult(*G2p, *G2p, NULL, GrB_PLUS_FP64, *G2p, *Matching, GrB_DESC_S);
        RoundOffAllValuesInMatrix(G2p, 6, NULL); // HOW MANY DECIMALS?
        // We only require the positive-gain 2-augmentations in G2p and D2p, so save these indices in g2p.
        GrB_select(*G2p, NULL, NULL, GrB_VALUEGT_FP64, *G2p, 0, NULL);

        GrB_reduce(g2p, NULL, NULL, GrB_PLUS_FP64, *G2p, NULL);
        // We want to use the rows from C indicated by v for making D2p
        SelectRowsInMatrixFromVector(&C, &v, NULL, n);

        // We want to use the rows from E indicated by a but not v for making D2p, and add the results of E and C.
        ApplyMaskToVector(&a, &v, &c, true);
        SelectRowsInMatrixFromVector(&E, &c, NULL, n);
        GrB_eWiseAdd(*D2p, NULL, NULL, GrB_PLUS_FP64, E, C, NULL);
        // Only leave the 2-augmentations that have positive gain.
        SelectRowsInMatrixFromVector(D2p, &g2p, NULL, n);
    }

    // --- Diagnostics.
    GrB_Index nPos2PathAugs;
    GrB_Matrix_nvals(&nPos2PathAugs, *G2p);
    if (GIVEOUTPUT == 2)
        printf("Searching path 2-augmentations results:\n - %llu positive-gain path 2-augmentations found, of which %llu unique\n", nPos2PathAugs, nPos2PathAugs / 2);
    

    
    // --- Free data.
    GrB_Vector_free(&FullVector);
    GrB_Vector_free(&m);
    GrB_Vector_free(&m_w);
    GrB_Vector_free(&m_u);
    GrB_Vector_free(&a);
    GrB_Vector_free(&a_g);
    GrB_Vector_free(&a_gi);
    GrB_Vector_free(&a_gj);
    GrB_Vector_free(&v);
    GrB_Vector_free(&c);
    GrB_Vector_free(&c_u);
    GrB_Vector_free(&c_l);
    GrB_Vector_free(&g2p);


    GrB_Matrix_free(&M_u);
    GrB_Matrix_free(&B);
    GrB_Matrix_free(&C);
    GrB_Matrix_free(&E);
    GrB_Matrix_free(&F);
    GrB_Matrix_free(&AwoM);
    GrB_Matrix_free(&G1);
    GrB_Matrix_free(&D1);


    return GrB_SUCCESS;
}

// Finds the best cycle 2-augmentation for each potential center i.e. each matched edge (i,j).
// Its gain value is stored in matrix G2c at (i,j) and explicit values are stored in D2c at 
// (i,k) and (j,l) to indicate the used 1-augmentations.
// Permanently changes G2c, D2c and the augmentation matrix.
// O(n+m) runtime.
GrB_Info SearchCycle2Augs (GrB_Matrix *Adjacency, GrB_Matrix *Matching, GrB_Matrix *Augmentation, 
                              GrB_Matrix *G2c, GrB_Matrix *D2c, GrB_Index n)
{
    // --- Allocate space for data.
    GrB_Vector FullVector, r, m, m_w, m_u;
    GrB_Vector_new(&FullVector, GrB_FP64, n);
    GrB_Vector_new(&r, GrB_FP64, n);
    GrB_Vector_new(&m, GrB_FP64, n);
    GrB_Vector_new(&m_w, GrB_FP64, n);
    GrB_Vector_new(&m_u, GrB_FP64, n);

    GrB_Matrix M_u, A_m, A_mM, D2cprime, C;
    GrB_Matrix_new(&M_u, GrB_FP64, n, n);
    GrB_Matrix_new(&A_m, GrB_FP64, n, n);
    GrB_Matrix_new(&A_mM, GrB_FP64, n, n);
    GrB_Matrix_new(&D2cprime, GrB_FP64, n, n);
    GrB_Matrix_new(&C, GrB_FP64, n, n);

    GrB_Index d2c_nvals;


    // --- Set up basic data.

    // Full vector constant size in bytes: 224 | column vector
    GrB_assign(FullVector, NULL, NULL, 1, GrB_ALL, n, NULL);
    
    // Create vector m with all indices of the match of vertices
    GrB_mxv(m, NULL, NULL, GxB_PLUS_SECONDI_INT64, *Matching, FullVector, NULL);

    // Create vector m_w with the weights of all matched edges.
    GrB_assign(m_w, NULL, NULL, 0, GrB_ALL, n, NULL);
    GrB_reduce(m_w, NULL, GrB_PLUS_FP64, GrB_PLUS_MONOID_FP64, *Matching, NULL);
    
    // Create vector m_u with explicit values for each matched pair (i,j) at i where i < j.
    GrB_select(M_u, NULL, NULL, GrB_TRIU, *Matching, 0, NULL);
    GrB_reduce(m_u, NULL, NULL, GrB_PLUS_MONOID_FP64, M_u, NULL);


    // --- Algorithm start.

    // Cycle augmentations consist fully of matched vertices so we select these rows and columns.
    SelectRowsInMatrixFromVector(Adjacency, &m, &A_m, n);
    SelectColumnsInMatrixFromVector(&A_m, &m, NULL, n);

    // To check if two unmatched edges connected to i and j are connected by another matched edge,
    // set up matrix Cstar with explicit values indicating the unmatched edges connected to j in row i
    PermuteMatrixRows(&A_m, Matching, NULL, &A_mM);
    PermuteMatrixColumns(&A_mM, Matching, NULL, NULL);

    // Element-wise wultiplying with the plus operator gives the sum of both values if they are explicit,
    // but an implicit value if any is implicit. The result signifies if an alternating cycle exists.
    GrB_eWiseMult(A_mM, NULL, NULL, GrB_PLUS_FP64, A_m, A_mM, GrB_DESC_R); 
    // Subtract the sum of both matched edges for each found cycle to obtain the gain of the cycle augmentation.
    MxMForcedDotMethod((GrB_Matrix)m_w, (GrB_Matrix)m_w, A_mM, GxB_PLUS_PLUS_FP64, &C);
    GrB_eWiseAdd(C, NULL, NULL, GrB_MINUS_FP64, A_mM, C, NULL);

    // Round-off errors can occur in the addition + subtraction.
    RoundOffAllValuesInMatrix(&C, 6, NULL);

    // Select the highest gain augmentation per center (i,j), but only in row i. 
    // Then only positive-gain augmentations are selected.
    SelectRowsInMatrixFromVector(&C, &m_u, D2c, n);
    SelectHighestValuePerRow(D2c, NULL, n);
    GrB_select(*D2c, NULL, NULL, GrB_VALUEGT_FP64, *D2c, 0, NULL);

    GrB_Matrix_nvals(&d2c_nvals, *D2c);

    // --- Early stop when no positive-gain augmentations exist.
    if (d2c_nvals > 0)
    {   
        // Permute the values in D2c by row and column with the matching and add to the original to obtain
        // the desired matrix with (i,k) and (j,l) explicit per cycle augmentation.
        PermuteMatrixRows(D2c, Matching, NULL, &D2cprime);
        PermuteMatrixColumns(&D2cprime, Matching, NULL, NULL);
        GrB_eWiseAdd(*D2c, NULL, NULL, GrB_PLUS_FP64, *D2c, D2cprime, NULL);

        // Create the corresponding values of positive-gain augmentations in G2c, with the gain of the augmentation at (i,j).
        GrB_reduce(r, NULL, NULL, GrB_PLUS_MONOID_FP64, *D2c, NULL);
        MxMForcedDotMethod((GrB_Matrix)r, (GrB_Matrix)FullVector, *Matching, GrB_MAX_FIRST_SEMIRING_FP64, G2c);
    }
    
    // --- Diagnostics.
    GrB_Index nPos2CycleAugs;
    GrB_Matrix_nvals(&nPos2CycleAugs, *G2c);
    if (GIVEOUTPUT == 2)
        printf("Searching cycle 2-augmentations results:\n - %llu positive-gain cycle 2-augmentations found, of which %llu unique\n", nPos2CycleAugs, nPos2CycleAugs / 2);
    

    // --- Free data.
    GrB_Vector_free(&FullVector);
    GrB_Vector_free(&r);
    GrB_Vector_free(&m);
    GrB_Vector_free(&m_w);
    GrB_Vector_free(&m_u);

    GrB_Matrix_free(&M_u);
    GrB_Matrix_free(&A_m);
    GrB_Matrix_free(&A_mM);
    GrB_Matrix_free(&D2cprime);
    GrB_Matrix_free(&C);

    return GrB_SUCCESS;
}


// The actual procedure for calculating a set of positive-gain 2-augmentations.
// We compare the best path and cycle augmentations found by the respective methods,
// choose the best ones and remove any conflicts.
// Permanently changes augmentation matrix.
// O(n+m) runtime.
GrB_Info Search2AugsProcedure(GrB_Matrix *Adjacency, GrB_Matrix *Matching, GrB_Matrix *Augmentation, 
                              LogEntry *logEntry, GrB_Index n)
{
    // --- Allocate space for data.
    GrB_Vector FullVector, m, cc, cp, r, r_comp, v, a_g;
    GrB_Vector_new(&FullVector, GrB_FP64, n);
    GrB_Vector_new(&m, GrB_FP64, n);
    GrB_Vector_new(&cc, GrB_FP64, n);
    GrB_Vector_new(&cp, GrB_FP64, n);
    GrB_Vector_new(&r, GrB_FP64, n);
    GrB_Vector_new(&r_comp, GrB_FP64, n);
    GrB_Vector_new(&v, GrB_FP64, n);
    GrB_Vector_new(&a_g, GrB_FP64, n);

    GrB_Matrix mdiag, G2, G2c, G2p, D2, D2c, D2p, Gamma_mask, Gamma, C;
    GrB_Matrix_new(&mdiag, GrB_FP64, n, n);
    GrB_Matrix_new(&G2, GrB_FP64, n, n);
    GrB_Matrix_new(&G2c, GrB_FP64, n, n);
    GrB_Matrix_new(&G2p, GrB_FP64, n, n);
    GrB_Matrix_new(&D2, GrB_FP64, n, n);
    GrB_Matrix_new(&D2c, GrB_FP64, n, n);
    GrB_Matrix_new(&D2p, GrB_FP64, n, n);
    GrB_Matrix_new(&Gamma_mask, GrB_FP64, n, n);
    GrB_Matrix_new(&Gamma, GrB_FP64, n, n);
    GrB_Matrix_new(&C, GrB_FP64, n, n);

    GrB_Index g2_nvals;


    // --- Set up basic data.

    // Full vector constant size in bytes: 224 | Column vector
    GrB_assign(FullVector, NULL, NULL, 1, GrB_ALL, n, NULL);

    // Create vector m with all indices of the match of vertices.
    GrB_mxv(m, NULL, NULL, GxB_PLUS_SECONDI_INT64, *Matching, FullVector, NULL);


    // --- Algorithm start.

    // Obtain the best cycle and path 2-augmentations (i,j,k,l) for each center (i,j).
    SearchCycle2Augs(Adjacency, Matching, Augmentation, &G2c, &D2c, n);
    SearchPath2Augs(Adjacency, Matching, Augmentation, &G2p, &D2p, n);
        

    // We choose from the path and cycle augmentation the one with the highest gain per center edge.
    // We make a mask C that indicates where we select cycle augmentations.
    GrB_assign(C, G2c, NULL, 1, GrB_ALL, n, GrB_ALL, n, NULL);
    ApplyMaskToMatrix(&C, &G2p, NULL, true);
    GrB_eWiseMult(C, NULL, GrB_PLUS_FP64, GrB_GT_FP64, G2c, G2p, NULL);
    GrB_select(C, NULL, NULL, GrB_VALUEEQ_FP64, C, 1, NULL);

    GrB_Matrix_dup(&G2, G2p);
    GrB_eWiseAdd(G2, C, NULL, GrB_FIRST_FP64, G2c, G2, NULL);    

    GrB_Matrix_nvals(&g2_nvals, G2);

    // --- Early stop when no positive-gain augmentations exist.
    if (g2_nvals > 0)
    {   
        // We save the gain of each best 2-augmentation in a_g
        GrB_reduce(a_g, NULL, NULL, GrB_PLUS_FP64, G2, NULL);

        // The locations of cycle augmentations are saved in vector r, and used to find the correct values in D2c.
        // We also create the complement of r, and thus finding all locations where we select path augmentations.
        GrB_reduce(r, NULL, NULL, GrB_PLUS_FP64, C, NULL);
        ApplyMaskToVector(&FullVector, &r, &r_comp, true);

        // We select the augmentations from D2c and D2p.
        SelectRowsInMatrixFromVector(&D2c, &r, NULL, n);
        SelectRowsInMatrixFromVector(&D2p, &r_comp, NULL, n);
        // Then we build matrix D2 combining both the best path and cycle augmentations.
        GrB_eWiseAdd(D2, NULL, NULL, GrB_PLUS_FP64, D2c, D2p, NULL);
        GrB_reduce(r, NULL, NULL, GrB_PLUS_FP64, D2, NULL);


        // --- Solve conflicts.

        // Make matrix Gamma with elements (i,i), (i,j), (i,k) and (i,l) per positive-gain 3-augmentation. See thesis.
        PermuteMatrixRows(&D2, Matching, NULL, &Gamma_mask); // (i,l)
        GrB_eWiseAdd(Gamma_mask, NULL, GrB_PLUS_FP64, GrB_PLUS_FP64, *Matching, D2, NULL); // (i,j), (i,k)
        GrB_Matrix_diag(&mdiag, m, 0); // (i,i)
        GrB_eWiseAdd(Gamma_mask, NULL, NULL, GrB_PLUS_FP64, Gamma_mask, mdiag, NULL); // (i,i)
        SelectRowsInMatrixFromVector(&Gamma_mask, &r, NULL, n);

        MxMForcedDotMethod((GrB_Matrix)a_g, (GrB_Matrix)FullVector, Gamma_mask, GxB_PLUS_FIRST_FP64, &Gamma);
        SelectHighestValuePerColumn(&Gamma, NULL, n);


        // Count the number of elements per row remaining and find where there are four left.
        GrB_mxv(v, NULL, NULL, GxB_PLUS_SECOND_FP64, Gamma, FullVector, NULL);
        GrB_select(v, NULL, NULL, GrB_VALUEEQ_FP64, v, 4, NULL);
        GrB_Vector_clear(r);
        PermuteVector(&v, Matching, NULL, &r);
        GrB_eWiseAdd(v, NULL, NULL, GrB_PLUS_FP64, v, r, NULL);

        // Add exactly those augmentations.
        SelectRowsInMatrixFromVector(&D2, &v, Augmentation, n);
        GrB_eWiseAdd(*Augmentation, NULL, NULL, GrB_PLUS_FP64, *Augmentation, *Augmentation, GrB_DESC_T1);
    }
    else
        return GrB_NO_VALUE;


    // --- Diagnostics.
    GrB_Index nPos2Augs, n2AugsNoConflicts;
    GrB_Matrix_nvals(&nPos2Augs, G2); 
    GrB_Matrix_nvals(&n2AugsNoConflicts, *Augmentation); n2AugsNoConflicts /= 4; // Unique 2-augmentations have 4 explicit values in the augmentation matrix, 2 for each unmatched edge.
    (*logEntry).nAugsFound = nPos2Augs;
    (*logEntry).nAugsNoConflicts = n2AugsNoConflicts;


    // --- Free data.
    GrB_Vector_free(&FullVector);
    GrB_Vector_free(&m);
    GrB_Vector_free(&cc);
    GrB_Vector_free(&cp);
    GrB_Vector_free(&r);
    GrB_Vector_free(&r_comp);
    GrB_Vector_free(&v);
    GrB_Vector_free(&a_g);

    GrB_Matrix_free(&mdiag);
    GrB_Matrix_free(&G2);
    GrB_Matrix_free(&G2c);
    GrB_Matrix_free(&G2p);
    GrB_Matrix_free(&D2);
    GrB_Matrix_free(&D2c);
    GrB_Matrix_free(&D2p);
    GrB_Matrix_free(&Gamma_mask);
    GrB_Matrix_free(&Gamma);
    GrB_Matrix_free(&C);

    return GrB_SUCCESS;
}
