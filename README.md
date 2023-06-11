# Maximum weight Matching Approximation Algorithm Using The GraphBLAS Standard

David D. de Best, 2023.

This C program implements a method of approximating the maximum weight matching (MWM) for general graphs
via positive-gain k-augmentations using the GraphBLAS standard.
It is the companion for the master thesis 'A 3/4-approximation for the maximum weight matching problem using the 
GraphBLAS standard', David D. de Best, supervised by prof. dr. Rob H. Bisseling, Utrecht University, 2023.

For an in-depth explanation of the procedures implemented, please see the thesis and provided comments 
in de code files.

----------------------------------------------------------------------------------------------------
## Software Dependencies

This program relies on:
 -  The SuiteSparse:GraphBLAS package by Timothy A. Davis, which is a full implementation
    of the GraphBLAS standard, see http://graphblas.org and https://people.engr.tamu.edu/davis/GraphBLAS.html.
 -  The LAGraph library, a collection of high level graph algorithms based on the GraphBLAS C API,
    see https://lagraph.readthedocs.io/en/latest/ and https://github.com/GraphBLAS/LAGraph.

To compile this C code, one needs to compile both SuiteSparse:GraphBLAS and LAGraph and add their folders to
the folder containing main.c.

----------------------------------------------------------------------------------------------------
## Compilation/Use

Requires compiled SuiteSparse:GraphBLAS and LAGraph folders in the same folder as main.c.
Compilation using some C compiler using the -llagraph and -lgraphblas flags, for example on MacOS:

    clang -o program main.c -llagraph -lgraphblas

Then either run the program from terminal straight as:

    ./program

Or provide input and output folders as first and second arguments, respectively:

    ./program /SomePath/Input /SomePath/Output

When not providing these folders, the program will ask for them in the terminal.
Next, the program will search for all graphs in Matrix Market format (.mtx) in the input folder
and perform the MWM algorithm one-by-one on each graph and write all output metadata to a single
.csv file named 'MWM_Log.csv' in the output folder.

---
Use of parallelism with OpenMP might cause problems during compilation on MacOS with the clang
compiler. Try using the GCC compiler and use the -fopenmp flag, for example:

    gcc main.c -o program -llagraph -lgraphblas -fopenmp

One might need to link the libraries via the terminal, something like:

    gcc main.c -o program -fopenmp -L /usr/local/lib -lgraphblas -llagraph -I /usr/local/include

Parallelism only works when SuiteSparse:GraphBLAS is compiled with OpenMP. For more on this and
other options on OpenMP parallelism, see the GraphBLAS User Guide.

---
Only compilation on MacOS has been tested, there are no guarantees for Linux/Windows systems but
with GraphBLAS, LAGraph and OpenMP installed compilation should likely be possible.

---
Calculation of the MWM approximation is done via strategy of selecting iterations.
The standard strategy is set to BASIC, but can be set to ALTERNATING or ONEAUG_PREFERENCE
in the generalheader.h file by changing the line:

    #define STRATEGY BASIC

to either one of the following:

    #define STRATEGY ALTERNATING
    #define STRATEGY ONEAUG_PREFERENCE

See the thesis for a description of each strategy.

---
Other options on what output is given in the terminal, debugmode and number of threads used can be 
changed in the generalheader.h file.

----------------------------------------------------------------------------------------------------
## Error messages

The program throws warnings when using GraphBLAS 8.0.0 with the latest version of LAGraph about the 
GraphBLAS JIT Kernel and user defined types like: "error: unknown type name 'VX1X2Struct'". These 
can be ignored, have no effect on the program outcome and are hopefully solved with the coming version 
of LAGraph.

----------------------------------------------------------------------------------------------------
