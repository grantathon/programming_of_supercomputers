// Simple tool that creates a binary file from a dat file

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
    // Check for valid input arguement
    if(argc != 2)
    {
        printf("ERROR: Pass file URI as arguement\n");
        return 1;
    }

    // Program variables
    FILE *inputFile, *outputFile;
    char newFileName[128];

    // Open input file
    printf("Reading input file...\n");
    inputFile = fopen(argv[1], "r");
    printf("Completed reading input file!\n\n");

    // Construct name for binary file to be created
    strcpy(newFileName, "test.bin");
    // strcpy(newFileName, argv[1]);
    // strcat(newFileName, ".bin");

    // Create copy of binary file in directory of input file
    printf("Creating output file...\n");
    outputFile = fopen(newFileName, "wb");

    // Copy input file contents to output file in binary format
    int element = fgetc(inputFile);
    while(element != EOF)
    {
        for(int i = 0; i < 8; i++)
        {
            fputc(((element >> i) & 1), outputFile);
        }

        element = fgetc(inputFile);
    }

    printf("Completed creating output file!\n\n");

    // Finalize program
    fclose(inputFile);
    fclose(outputFile);

    return 0;
}