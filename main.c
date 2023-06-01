/* C code document by David D. de Best for the master thesis 
   'Approximating the maximum weight matching problem with 
   positive-gain k-augmentations using the GraphBLAS standard', 
   supervised by prof. Rob H. Bisseling, Utrecht University, 2023 */


#include "Code files/generalheader.h"
#include "Code files/helperfunctions.c"
#include "Code files/customtype.c"
#include "Code files/IO.c"
#include "Code files/1augmentations.c"
#include "Code files/2augmentations.c"
#include "Code files/3augmentations.c"
#include "Code files/flip&search.c"


// LAGRAPH_CATCH is required by LAGRAPH_TRY.  If an error occurs, this macro
// catches it and takes corrective action, then terminates this program.
#define LAGRAPH_CATCH(info)                 \
{                                           \
    GrB_free (&Adjacency);                  \
    GrB_free (&Matching);                   \
    GrB_free (&Augmentation);               \
    return (info) ;                         \
}

// Print diagnostics info on the program, GraphBLAS version, number of threads used.
GrB_Info PrintProgramInfo()
{
    if (GIVEOUTPUT) 
        printf ("Program data:\n");
    int libraryVersion[3];
    bool openMPEnabled;
    GxB_get(GxB_LIBRARY_OPENMP, &openMPEnabled);
    GxB_get(GxB_LIBRARY_VERSION, &libraryVersion);

    if (!openMPEnabled && GIVEOUTPUT) 
        printf ("GraphBLAS not compiled with OpenMP\n");
    if (openMPEnabled && GIVEOUTPUT) 
        printf ("GraphBLAS compiled with OpenMP\n");
    if (GIVEOUTPUT)
        printf ("GraphBLAS version %d.%d.%d\n", libraryVersion[0], libraryVersion[1], libraryVersion[2]);
    
    int maxNThreadsDetected, nThreadsUsed;
    GxB_get (GxB_NTHREADS, &maxNThreadsDetected);
    if (GIVEOUTPUT)
        printf ("Maximum number of threads detected for use: %d\n", maxNThreadsDetected);
    
    GxB_set (GxB_NTHREADS, NTHREADSTOUSE);
    GxB_get (GxB_NTHREADS, &nThreadsUsed);
    if (GIVEOUTPUT)
        printf ("Number of threads used set to %d\n", nThreadsUsed);

    double chunkSize;
    GxB_get (GxB_CHUNK, &chunkSize);
    if (GIVEOUTPUT)
        printf ("Chunk size: %g\n", chunkSize);

    return GrB_SUCCESS;
}

/// The start of the program.
int main (int argc, char **argv) 
{
    // ----- Initialization -----

    if (GIVEOUTPUT)
    {
        printf("\n----------------\nStarted Program\n----------------\n");
        fflush(stdout);
    }

    char msg [LAGRAPH_MSG_LEN];
    GrB_Matrix Adjacency, Matching, Augmentation;
    LAGRAPH_TRY (LAGraph_Init (msg));

    PrintProgramInfo();

    // ----- I/O -----

    // Find input data.
    char inputFolderName[256];
    int inputMatrixCount;
    if (SetUpInputFolder(argc, argv, inputFolderName, &inputMatrixCount) == GrB_INVALID_VALUE)
    {
        printf("Given input folder not available...\n\n----------------\nExiting program\n----------------\n");
        LAGRAPH_TRY (LAGraph_Finalize (msg));
        return 0;
    }

    char matrixNames[inputMatrixCount][256];
    ReadInputMatrixNames(inputFolderName, matrixNames, inputMatrixCount);

    // Set up output location.
    char logFolderName[256];
    char logFileName[256];
    if (SetUpOutputFolder(argc, argv, logFolderName) == GrB_INVALID_VALUE)
    {
        printf("Given output folder not available...\n\n----------------\nExiting program\n----------------\n");
        LAGRAPH_TRY (LAGraph_Finalize (msg));
        return 0;
    }
    if (CreateLog(logFolderName, logFileName) == GrB_OUTPUT_NOT_EMPTY)
    {
        printf("No override of output file...\n\n----------------\nExiting program\n----------------\n");
        LAGRAPH_TRY (LAGraph_Finalize (msg));
        return 0;
    }

    char logFullLocation[513];
    strcpy(logFullLocation, logFolderName);
    strcat(logFullLocation, "/");
    strcat(logFullLocation, logFileName);



    // ----- Main loop over graphs in input folder -----

    for (size_t g = 0; g < inputMatrixCount; g++)
    {
        char currentMatrixLocation[513] = "";
        strcat(currentMatrixLocation, inputFolderName);
        strcat(currentMatrixLocation, "/");
        strcat(currentMatrixLocation, matrixNames[g]);

        if (GIVEOUTPUT)
            printf("\n------------------------------------------------------------\nPerforming Maximum-Weight Matching algorithm on %s\n------------------------------------------------------------\n", matrixNames[g]);
        

        // ----- Data structure setup -----

        FILE *f = fopen(currentMatrixLocation, "r");
        LAGRAPH_TRY(LAGraph_MMRead(&Adjacency, f, msg));
        
        RemoveDiagonalFromMatrix(&Adjacency, NULL);

        double matchingWeight = 0;
        GrB_Index n, m;
        GrB_Matrix_nrows(&n, Adjacency);
        GrB_Matrix_nvals(&m, Adjacency);
        GrB_Matrix_new(&Matching, GrB_FP64, n, n);
        GrB_Matrix_new(&Augmentation, GrB_FP64, n, n);

        if (GIVEOUTPUT)
            printf("Graph has %llu vertices and %llu edges\n", n, m/2);

        LogEntry logEntry = {matrixNames[g], 0, 'U', -1, -1, -1, -1, -1, -1, -1, -1};


        // ----- Performing algorithm -----

        bool finish[4] = {true, false, false, false};
        int k = 1;

        while (!(finish[1] &&  finish[2] && finish[3]))
        {   
            if (SearchKAugmentations(&Adjacency, &Matching, &Augmentation, k, &logEntry, n) == GrB_NO_VALUE)
            {
                finish[k] = true;

                k++;
                if (k == 4)
                    k = 1;
            }   
            else
            {
                AddNewLogEntry(logEntry, logFullLocation);
                NewIterationLogEntry(&logEntry);
                FlipAugmentations(&Adjacency, &Matching, &Augmentation, &matchingWeight, &logEntry, n);
                AddNewLogEntry(logEntry, logFullLocation);
                NewIterationLogEntry(&logEntry);

                if (STRATEGY == ALTERNATING)
                {
                    k++;
                    if (k == 4)
                        k = 1;
                }

                if (STRATEGY == ONEAUG_PREFERENCE)
                {
                    k = 1;
                }
                
                finish[1] = false;
                finish[2] = false;
                finish[3] = false;
            }
        }

        // Calculate a basic MWM upper bound. See thesis.
        double mwmUpperBound = 0;
        MWMUpperBound(Adjacency, &mwmUpperBound);
        if (GIVEOUTPUT)
            printf("\nMatching weight upper bound: %f\n", mwmUpperBound);

        // Free the graphs matrices for use for next graph.
        GrB_Matrix_free(&Adjacency);
        GrB_Matrix_free(&Matching);
        GrB_Matrix_free(&Augmentation);
    }
    
    
    // ----- Finish -----

    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;

    if (GIVEOUTPUT)
    {
        printf("\n-----------------\nFinished Program\n-----------------\n\n");
        fflush(stdout);
    }

   return 0;
}