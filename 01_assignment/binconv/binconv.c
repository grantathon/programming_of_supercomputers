// Simple tool that creates a binary file from a dat file

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
    // Check for valid input arguement
    if(argc != 3)
    {
        printf("ERROR: Pass both arguements <input file> <output file> as URIs\n");
        return 1;
    }

    // File variables
    FILE *inputFile, *outputFile;

    // Open input file
    printf("Reading input file...\n");
    inputFile = fopen(argv[1], "r");
    printf("Completed reading input file!\n\n");

    // Create copy of binary file in directory of input file
    printf("Creating output file...\n");
    outputFile = fopen(argv[2], "wb");

    // Copy input file contents to output file in binary format
    char element = fgetc(inputFile);
    while(element != EOF)
    {
        fputc(element, outputFile);
        element = fgetc(inputFile);
    }

    printf("Completed creating output file!\n");

    // Finalize program
    fclose(inputFile);
    fclose(outputFile);

    return 0;
}