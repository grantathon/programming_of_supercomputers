// Simple tool that converts a text file to a binary file

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INT_ROW_SIZE 6
#define DOUBLE_ROW_SIZE 8

int main(int argc, char *argv[])
{
    // Check for valid input arguement
    if(argc != 3)
    {
        printf("ERROR: Pass both arguements <input file> <output file> as URIs\n");
        return 1;
    }

    // File input variables
    FILE *inputFile, *outputFile;
    int startInners, endInners, startOuters, endOuters;
    int intInput[INT_ROW_SIZE];
    double doubleInput[DOUBLE_ROW_SIZE];

    // Open text input file
    printf("Begin converting text file to binary...\n");
    inputFile = fopen(argv[1], "r");

    // Create binary output file
    outputFile = fopen(argv[2], "wb");

    // Read the initial row parameters that describe file construction
    fscanf(inputFile, "%d", &startInners);
    fscanf(inputFile, "%d", &endInners);
    fscanf(inputFile, "%d", &startOuters);
    fscanf(inputFile, "%d", &endOuters);

    // Read and write the integer-based rows
    for(int i = startInners; i < startOuters; i++)
    {
        // Read the row elements
        for(int j = 0; j < INT_ROW_SIZE; j++)
        {
            fscanf(inputFile, "%d", &intInput[j]);
        }

        // Output the row elements
        fwrite(intInput, INT_ROW_SIZE, sizeof(*intInput), outputFile);
    }

    // Read and write the double-based rows
    for(int i = startInners; i < startOuters; i++)
    {
        // Read the row elements
        for(int j = 0; j < DOUBLE_ROW_SIZE; j++)
        {
            fscanf(inputFile, "%lf", &doubleInput[j]);
        }

        // Output the row elements
        fwrite(doubleInput, DOUBLE_ROW_SIZE, sizeof(*doubleInput), outputFile);
    }

    // Close files
    fclose(inputFile);
    fclose(outputFile);

    printf("Completed converting text file to binary!\n");

    return 0;
}