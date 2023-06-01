/* C code document by David D. de Best for the master thesis 
   'Approximating the maximum weight matching problem with 
   positive-gain k-augmentations using the GraphBLAS standard', 
   supervised by prof. Rob H. Bisseling, Utrecht University, 2023 */

#ifndef G_H
#define G_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <dirent.h>

// Optionally update location of these headers when the user cannot install them in the usual locations
#include "GraphBLAS.h"
#include "LAGraph.h"

// ----- USER SETTINGS -----

// Give the number of threads to use for program execution.
#define NTHREADSTOUSE 8

// GIVEOUTPUT: 0 for no output to console; 1 for basic information; 2 for extensive information per iteration.
#define GIVEOUTPUT 1
// DEBUGMODE provides extra checks per iteration on the validity of the augmentations and matching,
// turn of for best performance.
#define DEBUGMODE 0


#define BASIC 0
#define ALTERNATING 1
#define ONEAUG_PREFERENCE 2

// Select a strategy out of BASIC, ALTERNATING and ONEAUG_PREFERENCE. See thesis for definition.
#define STRATEGY BASIC


// Custom struct for saving data to csv file.
typedef struct LogEntry
{
    char *graphName;
    int iterationNumber;
    char iterationType;
    double runtime;
    int kAugSearch;
    int nAugsFound;
    int nAugsNoConflicts;
    double weightAdded;
    double weightRemoved;
    double weightImprovement;
    double oldWeight;
    double newWeight;
} LogEntry;


// Custom types for 3-augmentation algorithm.
typedef struct VX1X2Struct
{
    double value;
    int64_t x1;
    int64_t x2;
} VX1X2Struct;

typedef struct X1X2Struct
{
    int64_t x1;
    int64_t x2;
} X1X2Struct;

#endif

