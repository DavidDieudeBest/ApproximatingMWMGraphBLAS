/* C code document by David D. de Best for the master thesis 
   'Approximating the maximum weight matching problem with 
   positive-gain k-augmentations using the GraphBLAS standard', 
   supervised by prof. Rob H. Bisseling, Utrecht University, 2023 */


#include "generalheader.h"

// The procedure for calculating a set of positive-gain 3-augmentations.
// Permanently changes augmentation matrix.
// O(n+m) runtime.
GrB_Info Search3AugsProcedure (GrB_Matrix *Adjacency, GrB_Matrix *Matching, GrB_Matrix *Augmentation, LogEntry *logEntry, GrB_Index n)
{
    // --- Allocate space for data.
    GrB_Type ValueX1X2;
    GrB_Type_new(&ValueX1X2, sizeof(VX1X2Struct));

    GrB_BinaryOp MaxVX1X2;
    GrB_BinaryOp_new(&MaxVX1X2, (GxB_binary_function) MaxVX1X2Method, ValueX1X2, ValueX1X2, ValueX1X2);

    GrB_Vector FullVector, ZeroVector, m, m_w, kl_i, c, v, a_g;
    GrB_Vector_new(&FullVector, GrB_FP64, n);
    GrB_Vector_new(&ZeroVector, GrB_FP64, n);
    GrB_Vector_new(&m, GrB_FP64, n);
    GrB_Vector_new(&m_w, GrB_FP64, n);
    GrB_Vector_new(&kl_i, GrB_FP64, n);
    GrB_Vector_new(&c, GrB_FP64, n);
    GrB_Vector_new(&v, GrB_FP64, n);
    GrB_Vector_new(&a_g, GrB_FP64, n);

    GrB_Matrix AwoM_star, G1, C, K, L, Cstar, F, G3, G3_V, G3_P, G3_Q, D3, Gamma, Gamma_mask, mdiag, centers;
    GrB_Matrix_new(&AwoM_star, GrB_FP64, n, n);
    GrB_Matrix_new(&G1, GrB_FP64, n, n);
    GrB_Matrix_new(&C, GrB_FP64, n, n);
    GrB_Matrix_new(&Cstar, GrB_FP64, n, n);
    GrB_Matrix_new(&F, ValueX1X2, n, n);
    GrB_Matrix_new(&K, GrB_FP64, n, n);
    GrB_Matrix_new(&L, GrB_FP64, n, n);
    GrB_Matrix_new(&G3, ValueX1X2, n, n);
    GrB_Matrix_new(&G3_V, GrB_FP64, n, n);
    GrB_Matrix_new(&G3_P, GrB_INT64, n, n);
    GrB_Matrix_new(&G3_Q, GrB_INT64, n, n);
    GrB_Matrix_new(&D3, GrB_FP64, n, n);
    GrB_Matrix_new(&Gamma, GrB_FP64, n, n);
    GrB_Matrix_new(&Gamma_mask, GrB_FP64, n, n);
    GrB_Matrix_new(&mdiag, GrB_FP64, n, n);
    GrB_Matrix_new(&centers, GrB_FP64, n, n);

    GrB_Vector *k, *l, *kl_v;
    k = malloc(4 * sizeof(GrB_Vector));
    l = malloc(4 * sizeof(GrB_Vector));
    kl_v = malloc(4 * sizeof(GrB_Vector));

    GrB_Matrix *G_max, *Lprime, *Kprime;
    G_max = malloc(4 * sizeof(GrB_Matrix));
    Lprime = malloc(4 * sizeof(GrB_Matrix));
    Kprime = malloc(4 * sizeof(GrB_Matrix));

    for (size_t i = 0; i < 4; i++)
    {
        GrB_Vector_new(&k[i], ValueX1X2, n);
        GrB_Vector_new(&l[i], ValueX1X2, n);
        GrB_Vector_new(&kl_v[i], GrB_FP64, n);

        GrB_Matrix_new(&G_max[i], GrB_FP64, n, n);
        GrB_Matrix_new(&Lprime[i], GrB_FP64, n, n);
        GrB_Matrix_new(&Kprime[i], GrB_FP64, n, n);
    }

    GrB_Index g3_nvals;


    // --- Set up basic data.

    // Full vector constant size in bytes: 224 | Column vector
    GrB_assign(FullVector, NULL, NULL, 1, GrB_ALL, n, NULL);
    GrB_assign(ZeroVector, NULL, NULL, 0, GrB_ALL, n, NULL);

    // Create vector m with all indices of the match of vertices.
    GrB_mxv(m, NULL, NULL, GxB_PLUS_SECONDI_INT64, *Matching, FullVector, NULL);
    
    // We create full vector m_w to contain the weights of the matched edges.
    GrB_assign(m_w, NULL, NULL, 0, GrB_ALL, n, NULL);
    GrB_reduce(m_w, NULL, GrB_PLUS_FP64, GrB_PLUS_MONOID_FP64, *Matching, NULL);


    // --- Algorithm start.

    // We produce matrix [A \ M] containing all unmatched edges.
    ApplyMaskToMatrix(Adjacency, Matching, &AwoM_star, true);

    // Each explicit value of AwoM_star should point to a potential center (i,j) of a 3-augmentation,
    // so both vertices of the unmatched edge should be matched. We select these values.
    SelectRowsInMatrixFromVector(&AwoM_star, &m, NULL, n);

    // See 1-augmentation algorithm, we calculate all 1-augmentations for matched vertices.
    // Performed in between selection for matched vertices, since we want the 1-augmentations where
    // k in center (k,p) is matched, but not necessarily p.
    MxMForcedDotMethod((GrB_Matrix)m_w, (GrB_Matrix)m_w, AwoM_star, GxB_PLUS_PLUS_FP64, &C);
    GrB_eWiseAdd(G1, NULL, NULL, GrB_MINUS_FP64, AwoM_star, C, NULL);
    GrB_Matrix_clear(C);

    // We complete selection of the centers (i,j) where j is matched as well.
    SelectColumnsInMatrixFromVector(&AwoM_star, &m, NULL, n);

    // We need the top four highest-gain 1-augmentations per vertex in G_max, found by
    // by repeating the algorithm for selecting the maximum in rows.
    for (size_t i = 0; i < 4; i++)
    {
        SelectHighestValuePerRow(&G1, &G_max[i], n);
        ApplyMaskToMatrix(&G1, &G_max[i], NULL, true);
    }

    // We construct two matrices K and L. In K we indicate per center (i,j) the vertices of one of the top four k-arms that cannot be part of any l-arm in row i.
    // These are vertices i, k and p, and so (i,i), (i,k) and (i,p) are set explicit.
    // In L we indicate the vertices of one of the top four l-arms that cannot be part of any k-arm in column j.
    // These are vertices j, l and q, and so (j,j), (l,j) and (q,j) are set explicit.

    // (i,i) and (j,j) are found by taking all rows from AwoM_star indicating the centers and diagonalizing. 
    GrB_reduce(c, NULL, NULL, GxB_ANY_FP64_MONOID, AwoM_star, NULL);
    GrB_Matrix_diag(&C, c, 0);
    // (i,k) is found by permuting the columns of this diagonal matrix with the matching.
    PermuteMatrixColumns(&C, Matching, NULL, &K);
    GrB_eWiseAdd(K, NULL, NULL, GrB_PLUS_FP64, C, K, NULL);
    // (l,j) is found by permuting the rows of this diagonal matrix with the matching.
    PermuteMatrixRows(&C, Matching, NULL, &L);
    GrB_eWiseAdd(L, NULL, NULL, GrB_PLUS_FP64, C, L, NULL);

    // In a double for-loop we combine all options for top four k- and l-arms, so we re-assign elements (i,p) and (q,j) each iteration, corresponding to some arm.
    // These elements are retrieved from the G_max matrices. These contain (k,p) and (l,q), so these need to be permuted by the matching to obtain (i,p) and (j,q)
    // First select all the correct rows, indicated by permuting vector c.
    PermuteVector(&c, Matching, NULL, NULL);
    
    // We perform the work only on the upper triangular part of the symmetric matrix containing centers.
    GrB_select(centers, NULL, NULL, GrB_TRIU, AwoM_star, 0, NULL);

    // The four versions of K and L are finished here, with the top four p and q vertices.
    // The gain value and indices are saved in the user-defined type vectors k and l
    for (size_t x_1 = 0; x_1 < 4; x_1++)
    {
        GrB_Matrix_clear(C);
        GrB_Vector_clear(kl_i);
        SelectRowsInMatrixFromVector(&G_max[x_1], &c, &C, n);
        PermuteMatrixRows(&C, Matching, NULL, NULL);
        GrB_reduce(kl_v[x_1], NULL, NULL, GxB_ANY_FP64_MONOID, C, NULL);
        GrB_mxv(kl_i, NULL, NULL, GxB_PLUS_SECONDI_INT64, C, FullVector, NULL);
        ConstructVX1X2Vector(&kl_v[x_1], &kl_i, &ZeroVector, &k[x_1], ValueX1X2, n);
        ConstructVX1X2Vector(&kl_v[x_1], &ZeroVector, &kl_i, &l[x_1], ValueX1X2, n);
        GrB_eWiseAdd(Kprime[x_1], NULL, NULL, GrB_PLUS_FP64, K, C, NULL);
        GrB_eWiseAdd(Lprime[x_1], NULL, NULL, GrB_PLUS_FP64, L, C, GrB_DESC_T1);
        SelectRowsInMatrixFromVector(&Kprime[x_1], &k[x_1], NULL, n);
        SelectColumnsInMatrixFromVector(&Lprime[x_1], &l[x_1], NULL, n);
    }

    // Here we combine each k-arm with each l-arm. Multiplication of K and L sees where combination leads to conflict
    // These locations where no conflicts occur are saved in Cstar, and here we combine the arms via k and l
    // that is saved in F. In G3 we save the current best 3-augmentation found for each center.
    for (size_t x_1 = 0; x_1 < 4; x_1++)
    {
        for (size_t x_2 = 0; x_2 < 4; x_2++)
        {
            // The first iteration we found the locations where even the best arms do not combine into a positive-
            // gain 3-augmentation. These centers can be taken away for the rest of the algorithm.
            if (x_1 == 0 && x_2 == 0)
            {
                MxMForcedDotMethod((GrB_Matrix)kl_v[x_1], (GrB_Matrix)kl_v[x_2], centers, GxB_PLUS_PLUS_FP64, &centers);
                GrB_eWiseAdd(centers, centers, NULL, GrB_PLUS_FP64, centers, AwoM_star, NULL);
                GrB_select(centers, NULL, NULL, GrB_VALUEGT_FP64, centers, 0, NULL);
            }
            if (x_2 == 0)
                GrB_Matrix_dup(&Cstar, centers);
            else
                GrB_Matrix_dup(&Cstar, C);

            MxMForcedDotMethod2(Kprime[x_1], Lprime[x_2], Cstar, GxB_ANY_PAIR_FP64, &C);
            ApplyMaskToMatrix(&Cstar, &C, NULL, true);
            MultiplyVX1X2VecsToMatr(&k[x_1], &l[x_2], &Cstar, ValueX1X2, &F);


        }
        // We save the best l-arm combination with the k-arm of current x_1 iteration where it 
        // improves the best combination of previous x_1 combination.
        GrB_eWiseAdd(G3, NULL, NULL, MaxVX1X2, G3, F, NULL);    

        GrB_Matrix_clear(F);  
    }

    // We split the tuple types into individual matrices.
    SplitVX1X2Matrix(&G3, ValueX1X2, &G3_V, &G3_P, &G3_Q);

    // To obtain the true gain value per 3-augmentation, we need to add the weight of center edge (i,j)
    GrB_eWiseAdd(G3_V, G3_V, NULL, GrB_PLUS_FP64, G3_V, AwoM_star, NULL);
    // We select all positive-gain 3-augmentations and since we can add at most one per starting vertex i,
    // we select the highest-gain 3-augmentation per row.
    RoundOffAllValuesInMatrix(&G3_V, 6, NULL);
    GrB_select(G3_V, NULL, NULL, GrB_VALUEGT_FP64, G3_V, 0, NULL);
    GrB_Matrix_nvals(&g3_nvals, G3_V);

    GrB_Index nPos3Augs, n3AugsNoConflicts;

    // --- Early stop when no positive-gain augmentations exist.
    if (g3_nvals > 0)
    {   
        // Since we worked only on the upper triangular part, we need to copy over the transpose.
        // Note: G3_P^T should be added to G3_Q and G3_Q^T to G3_P
        GrB_eWiseAdd(G3_V, NULL, NULL, GrB_PLUS_FP64, G3_V, G3_V, GrB_DESC_T1);
        GrB_eWiseAdd(G3_P, NULL, NULL, GrB_PLUS_FP64, G3_P, G3_Q, GrB_DESC_T1);
        GrB_transpose(G3_Q, NULL, NULL, G3_P, NULL);


        SelectHighestValuePerRow(&G3_V, NULL, n);
        // We require the corresponding values from G3_P and G3_Q.
        ApplyMaskToMatrix(&G3_P, &G3_V, NULL, false);
        ApplyMaskToMatrix(&G3_Q, &G3_V, NULL, false);

        GrB_Vector_clear(c);
        GrB_Vector_clear(v);  
        // Here we transfrom G3_P and G3_Q into the D3 matrix, with explicit values per remaining 3-augmentation at (i,p) and (i,q).
        GrB_reduce(c, NULL, NULL, GxB_ANY_FP64, G3_P, NULL);
        GrB_reduce(v, NULL, NULL, GxB_ANY_FP64, G3_Q, NULL); 
        GrB_Matrix_clear(G3_P); 
        GrB_Matrix_clear(G3_Q);
        CreateMatrixFromColumnIndices(&c, &G3_P);
        CreateMatrixFromColumnIndices(&v, &G3_Q);

        // We save the gain of each best 3-augmentation in a_g
        GrB_reduce(a_g, NULL, NULL, GxB_ANY_FP64, G3_V, NULL);

        // --- Solve conflicts.

        GrB_Vector_clear(c);
        GrB_reduce(c, NULL, NULL, GxB_ANY_FP64, G3_V, NULL);

        // We make matrix Gamma with elements (i,i), (i,j), (i,k), (i,l), (i,p) and (i,q) per positive-gain 3-augmentation. See thesis.
        GrB_Matrix_diag(&mdiag, m, 0); // i,i
        PermuteMatrixColumns(&G3_V, Matching, NULL, &Gamma_mask);  // i,l
        GrB_eWiseAdd(Gamma_mask, NULL, NULL, GrB_PLUS_FP64, Gamma_mask, G3_Q, NULL); // i,l , i,q
        GrB_eWiseAdd(Gamma_mask, NULL, GrB_PLUS_FP64, GrB_PLUS_FP64, *Matching, G3_P, NULL); // i,k , i,p
        GrB_eWiseAdd(Gamma_mask, NULL, GrB_PLUS_FP64, GrB_PLUS_FP64, G3_V, mdiag, NULL); // i,j , i,i

        // Only keep the rows corresponding to a 3-augmentation, now with six elements per row
        SelectRowsInMatrixFromVector(&Gamma_mask, &c, NULL, n);

        // Create matrix Gamma with the gain value of each 3-augmentation in the six locations per row
        MxMForcedDotMethod((GrB_Matrix)a_g, (GrB_Matrix)FullVector, Gamma_mask, GxB_PLUS_FIRST_FP64, &Gamma);
        SelectHighestValuePerColumn(&Gamma, NULL, n);

        // Count the number of elements per row remaining and find where there are six left.
        GrB_mxv(v, NULL, NULL, GxB_PLUS_SECOND_FP64, Gamma, FullVector, NULL);
        GrB_select(v, NULL, NULL, GrB_VALUEEQ_FP64, v, 6, NULL);

        GrB_Matrix_nvals(&nPos3Augs, G3_V);
    
        SelectRowsInMatrixFromVector(&G3_P, &v, NULL, n);
        SelectRowsInMatrixFromVector(&G3_Q, &v, NULL, n);
        SelectRowsInMatrixFromVector(&G3_V, &v, NULL, n);
        GrB_Matrix_clear(C);
        GrB_transpose(C, NULL, NULL, G3_V, NULL);

        // We create matrix D3 with explicit values (k,p) and (l,q) per found 3-augmentation.
        PermuteMatrixRows(&G3_P, Matching, NULL, NULL);
        PermuteMatrixRows(&G3_Q, &C, NULL, NULL);
        PermuteMatrixRows(&G3_Q, Matching, NULL, NULL);
        
        GrB_eWiseAdd(D3, NULL, NULL, GrB_PLUS_FP64, G3_P, G3_Q, NULL);


        // We need the augmentation matrix to have (k,p), (l,q) and (i,j).
        // So here add G3 to the augmentation for edges (i,j).
        GrB_eWiseAdd(*Augmentation, NULL, NULL, GrB_PLUS_FP64, D3, G3_V, NULL); 
        // Add its own transpose to the augmentation for the (p,k), (q,l) and (j,i) elements.
        GrB_eWiseAdd(*Augmentation, NULL, NULL, GrB_PLUS_FP64, *Augmentation, *Augmentation, GrB_DESC_T1);
    }
    else
        return GrB_NO_VALUE;


    // --- Diagnostics.
    GrB_Matrix_nvals(&n3AugsNoConflicts, *Augmentation); n3AugsNoConflicts /= 6; // A unique 3-augmentation has six values in the augmentation matrix, two for each unmatched edge.
    (*logEntry).nAugsFound = nPos3Augs;
    (*logEntry).nAugsNoConflicts = n3AugsNoConflicts;


    // --- Free data.
    GrB_Type_free(&ValueX1X2);

    GrB_Vector_free(&FullVector); 
    GrB_Vector_free(&ZeroVector);   
    GrB_Vector_free(&m);   
    GrB_Vector_free(&m_w);
    GrB_Vector_free(&kl_i);
    GrB_Vector_free(&c);
    GrB_Vector_free(&v);
    GrB_Vector_free(&a_g);

    GrB_Matrix_free(&AwoM_star);
    GrB_Matrix_free(&G1);
    GrB_Matrix_free(&C);
    GrB_Matrix_free(&Cstar);
    GrB_Matrix_free(&F);
    GrB_Matrix_free(&K);
    GrB_Matrix_free(&L);
    GrB_Matrix_free(&G3);
    GrB_Matrix_free(&G3_V);
    GrB_Matrix_free(&G3_P);
    GrB_Matrix_free(&G3_Q);
    GrB_Matrix_free(&D3);
    GrB_Matrix_free(&Gamma);
    GrB_Matrix_free(&Gamma_mask);
    GrB_Matrix_free(&centers);

    for (size_t i = 0; i < 4; i++)
    {
        GrB_Vector_free(&k[i]);
        GrB_Vector_free(&l[i]);
        GrB_Vector_free(&kl_v[i]);

        GrB_Matrix_free(&G_max[i]);
        GrB_Matrix_free(&Kprime[i]);
        GrB_Matrix_free(&Lprime[i]);
    }

    free(k);
    free(l);
    free(kl_v);

    free(G_max);
    free(Kprime);
    free(Lprime);

    return GrB_SUCCESS;
}

