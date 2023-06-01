/* C code document by David D. de Best for the master thesis 
   'Approximating the maximum weight matching problem with 
   positive-gain k-augmentations using the GraphBLAS standard', 
   supervised by prof. Rob H. Bisseling, Utrecht University, 2023 */


#include "generalheader.h"


// ----- LOG MANAGEMENT -----

// Add the data of a finished iteration to the log with the name and path in logFullLocation.
GrB_Info AddNewLogEntry (LogEntry newEntry, char *logFullLocation)
{
    FILE *log;
    log = fopen(logFullLocation, "a");

    if (GIVEOUTPUT == 2)
        printf("\nAdding new log entry to file: \"%s\"\n", logFullLocation);

    fprintf(log, "%s;%d;%c;%f;%d;%d;%d;%f;%f;%f;%f;%f\n", newEntry.graphName, 
                                                            newEntry.iterationNumber, 
                                                            newEntry.iterationType, 
                                                            newEntry.runtime, 
                                                            newEntry.kAugSearch,
                                                            newEntry.nAugsFound,
                                                            newEntry.nAugsNoConflicts,
                                                            newEntry.weightAdded, 
                                                            newEntry.weightRemoved, 
                                                            newEntry.weightImprovement, 
                                                            newEntry.oldWeight, 
                                                            newEntry.newWeight);

    fclose(log);

    return GrB_SUCCESS;
}

// Reset a given LogEntry struct to default values to prepare for a new iteration.
GrB_Info NewIterationLogEntry (LogEntry *logEntry)
{
    (*logEntry).iterationNumber++;
    (*logEntry).iterationType = 'U';
    (*logEntry).runtime = -1;
    (*logEntry).kAugSearch = -1;
    (*logEntry).nAugsFound = -1;
    (*logEntry).nAugsNoConflicts = -1;
    (*logEntry).weightAdded = -1;
    (*logEntry).weightRemoved = -1;
    (*logEntry).weightImprovement = -1;
    (*logEntry).oldWeight = -1;
    (*logEntry).newWeight = -1;

    return GrB_SUCCESS;
}

// Create a new log file in the given location of logFolderName, and return the logs name as logName.
GrB_Info CreateLog (char *logFolderName, char *logName)
{
    // Create the full location of the new log.
    FILE *log;
    char logFileName[] = "MWM_Log.csv";
    char logFullLocation[513] = "";
    strcat(logFullLocation, logFolderName);
    strcat(logFullLocation, "/");
    strcat(logFullLocation, logFileName);

    // Check if a log with the same name already exist, if so ask if the user wants
    // to override the file.
    if (fopen(logFullLocation, "r") != NULL)
    {
        bool c = true;
        while (c)
        {
            char overwrite;
            printf("\nFilename is already in use. Do you want to overwrite the data in the file? [Y/N]: ");
            scanf("%c", &overwrite);
            if (overwrite == 'N')
                return GrB_OUTPUT_NOT_EMPTY;
            else if (overwrite == 'Y')
                c = false;
        }
            
    }

    // Create the new log and write the column names on the first line as a header.
    if (GIVEOUTPUT)
        printf("\nSaving log data in file: \"%s\"\n", logFullLocation);
    log = fopen(logFullLocation, "w");
    fprintf(log, "graph name;iteration number;iteration type;runtime;k-augs searched;#k-augs found;#k-augs no conflicts;weight added;weight removed;weight improvement;old matching weight;new matching weight\n");

    fclose(log);
    
    // Return the new log name as the logName variable.
    strcpy(logName, logFileName);

    return GrB_SUCCESS;
}


// ----- FILE MANAGEMENT -----

// Check if a given folder exists and can be opened.
GrB_Info CheckFolder (char *folderName)
{
    DIR *directory;
    struct dirent *directoryStream;

    if ((directory = opendir(folderName)) == NULL) 
    {
        printf("Problem opening folder '%s'\n", folderName);
        return GrB_INVALID_VALUE;
    }

    closedir(directory);
    
    return GrB_SUCCESS;
}

// Check if a given file has the Matrix Market extension name.
bool CheckFileMTXExtension (char *fileName)
{
    int fileNameLength = strlen(fileName);
    char mtxExt[] = ".mtx"; 
    char extension[5]; 
    memcpy(extension, &fileName[fileNameLength - 4], 4 * sizeof(char));
    
    if (strcmp(extension, ".mtx") == 0)
        return true;
    
    return false;
}

// Find and save the input folder given via the program arguments or ask the user.
// Check also the number of graphs in the input folder
GrB_Info SetUpInputFolder (int argc, char **argv, char *folderName, int *inputMatrixCount)
{
    // If an input folder is given in the program arguments, then check if it is correct.
    // If it is not given, ask the user for the folder name.
    if (argc > 1)
    {
        if (CheckFolder(argv[1]) == GrB_SUCCESS) 
            strcpy(folderName, argv[1]);
    }
    else
    {
        bool c = true;
        while (c)
        {
            char inputName[256];
            printf("\nPlease give the name of the folder containing all input matrices: ");
            scanf("%s", inputName);
            if (CheckFolder(inputName) == GrB_SUCCESS) 
            {
                strcpy(folderName, inputName);
                printf("Now working from input folder '%s'\n", folderName);
                c = false;
            }
        }
    }

    DIR *directory;
    struct dirent *currentDirectoryEntry;

    // Open the given input folder and count the number of files with the .mtx extension.
    directory = opendir(folderName);
    *inputMatrixCount = 0;
    do 
    {
        if ((currentDirectoryEntry = readdir(directory)) != NULL) 
            if (currentDirectoryEntry->d_type == DT_REG && CheckFileMTXExtension(currentDirectoryEntry->d_name))
                (*inputMatrixCount)++;
    } 
    while (currentDirectoryEntry != NULL);
    closedir(directory);

    if (GIVEOUTPUT)
        printf("\nPerforming Maximum-Weight Matching algorithm on %d matrices in input folder\n", *inputMatrixCount);

    return GrB_SUCCESS;
}

// Read and save all names of graphs in the given input folder.
GrB_Info ReadInputMatrixNames (char *folderName, char matrixNames[][256], int inputMatrixCount)
{
    DIR *directory;
    struct dirent *currentDirectoryEntry;

    directory = opendir(folderName);

    int c = 0;
    do 
    {
        if ((currentDirectoryEntry = readdir(directory)) != NULL) 
            if (currentDirectoryEntry->d_type == DT_REG && CheckFileMTXExtension(currentDirectoryEntry->d_name))
            {
                strcpy(matrixNames[c], currentDirectoryEntry->d_name);
                c++;
            }
    }
    while (currentDirectoryEntry != NULL);
    closedir(directory);
    
    return GrB_SUCCESS;
}

// Find and save the output folder given via the program arguments or ask the user.
GrB_Info SetUpOutputFolder (int argc, char **argv, char *logFolderName)
{    
    // If an output folder is given in the program arguments, then check if it is correct.
    // If it is not given, ask the user for the folder name.
    if (argc > 2)
    {
        strcpy(logFolderName, argv[2]);
        return(CheckFolder(logFolderName));
    }
    else
    {
        bool c = true;
        while (c)
        {
            char inputName[256];
            printf("\nPlease give the name of the folder where the log should be stored: ");
            scanf("%s", inputName);
            if (CheckFolder(inputName) == GrB_SUCCESS) 
            {
                strcpy(logFolderName, inputName);
                c = false;
            }
        }
    }

    return GrB_SUCCESS;
}   