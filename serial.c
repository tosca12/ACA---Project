#include <string.h>
#include <stdio.h>
#include <stdlib.h>


void binThreshold(int** imageMatrix, int value, int width, int height){
/*
    int rows=height;
    int cols=width;
*/
    for(int i=0; i<height; i++){
        for(int j=0; j<width; j++){

            if (imageMatrix[i][j] < value)
                imageMatrix[i][j] = 1;
            else
                imageMatrix[i][j] = 0;
        }
    }
}

void binComplement(int** imageMatrix, int width, int height){
    /*
    int rows=height;
    int cols=width;
*/

    for(int i=0; i<height; i++){
        for(int j=0; j<width; j++){
                imageMatrix[i][j] = (imageMatrix[i][j] -1) * (-1);
        }
    }
}


void binErosion(int** eroded, int width, int height){

        int** image = malloc(height * sizeof(int*));
        
        for (int i = 0; i < height; i++) {
            image[i] = malloc(width * sizeof(int));
            for (int j = 0; j < width; j++) {
                image[i][j] = eroded[i][j];
            }
        }
	
        for(int i = 1; i < height -1 ; i++){
            for(int j = 1; j < width -1; j++){
		
		if(image[i][j] * image[i+1][j] * image[i+1][j+1] *
		   image[i][j+1] * image[i-1][j+1] * image[i-1][j] *
		   image[i-1][j-1] * image[i][j-1] * image[i+1][j-1] != 0);
		   
		   else{
		            eroded[i][j] = 0;
		            eroded[i+1][j]= 0;
		            eroded[i+1][j+1]= 0;
		            eroded[i][j+1]= 0;
		            eroded[i-1][j+1]= 0;
		            eroded[i-1][j]= 0;
		            eroded[i-1][j-1]= 0;
		            eroded[i][j-1]= 0;
		            eroded[i+1][j-1]= 0; 
		    }
            }
        }    

        for (int i = 0; i < height; i++) {
            free(image[i]);
        }
        
        free(image);
 
}


void binDilation(int** dilated, int width, int height){

        int** image = malloc(height * sizeof(int*));
        
            for (int i = 0; i < height; i++) {
              image[i] = malloc(width * sizeof(int));
             for (int j = 0; j < width; j++) {
                    image[i][j] = dilated[i][j];
              }
          }

        for(int i = 1; i < height -1 ; i++){
            for(int j = 1; j < width -1; j++){
	
		if((image[i][j] + image[i+1][j] + image[i+1][j+1] +
		   image[i][j+1] + image[i-1][j+1] + image[i-1][j] +
		   image[i-1][j-1] + image[i][j-1] + image[i+1][j-1]) == 0);
		   
		   else{
		            dilated[i][j] = 1;
		            dilated[i+1][j]= 1;
		            dilated[i+1][j+1]= 1;
		            dilated[i][j+1]= 1;
		            dilated[i-1][j+1]= 1;
		            dilated[i-1][j]= 1;
		            dilated[i-1][j-1]= 1;
		            dilated[i][j-1]= 1;
		            dilated[i+1][j-1]= 1; 
		   }
	     }
          }    

        for (int i = 0; i < height; i++) {
            free(image[i]);
        }
        
        free(image);

}


void binOpening(int** imageMatrix, int width, int height){

    binErosion(imageMatrix, width, height);

    binDilation(imageMatrix, width, height);

}

void identifyBorders(int** imageMatrix, int width, int height) {

    int** image = malloc(height * sizeof(int*));
        for (int i = 0; i < height; i++) {
            image[i] = malloc(width * sizeof(int));
             for (int j = 0; j < width; j++) {
                image[i][j] = imageMatrix[i][j];
              }
        }

    binErosion(image, width, height);

    // Calculate the difference between the original and eroded images
    for(int i=0; i<height; i++){
        for(int j=0; j<width; j++){
                imageMatrix[i][j] = imageMatrix[i][j] - image[i][j];
        }
    }

    for (int i = 0; i < height; i++) {
            free(image[i]);
    }
    
        free(image);

}


int **saveImage (FILE *inputImage, int width, int height) {

    int **matrix = (int **)malloc(height * sizeof(int*));
    
    for (int i = 0; i < height; i++) {
        matrix[i] = (int *)malloc(width * sizeof(int));
    }

  
    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {
          matrix[i][j]=fgetc(inputImage);
      }
    }
  
  return matrix;
}


void writeImage(int **matrix, int width, int height, int maxVal, const char *outputFileName) {

    FILE *outputImage = fopen(outputFileName, "wb");
    
    if (outputImage == NULL) {
        printf("Error: Could not open %s for writing.\n", outputFileName);
        exit(1);
    }

    // Write PGM header
    fprintf(outputImage, "P2\n");
    fprintf(outputImage, "%d %d\n", width, height);
    fprintf(outputImage, "%d\n", maxVal); 

    // Write pixel values
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            fprintf(outputImage, "%d ", matrix[i][j]);
        }
        
        fprintf(outputImage, "\n");
    }

    fclose(outputImage);
}
    
  
  int main(int argc, char*argv[]) {

    FILE *inputImage;
    char magicNumber[3];
    int width, height, maxVal;
    int **imageMatrix;

    // The program requires the name of the image file to be elaborated and the threshold value that will be used; if the number of provided arguments is wrong, an error message is printed  
    if (argc != 3) {
        fprintf(stderr, "Error: Wrong number of arguments\nUsage: ./executable_name input_image.pgm <threshold>\n");
        return 1;
    }

    // Opening the input image to be elaborated; if the fopen returns a NULL value, the opening of the file was not successful and an error message is printed
    
    inputImage = fopen(argv[1], "rb");
    
    if (inputImage == NULL) {
        printf("Error: %s can't be opened\n", argv[1]);
        return 1;
    }

    // Reading the header of the input image; if the magic number of the header is different from P5 (magic number for the PGM file format), an error message is printed 
    
    fscanf(inputImage, "%s", magicNumber);
    
    if (magicNumber[0] != 'P' || magicNumber[1] != '5') {
        printf("Error: file format not supported\n");
        return 1;
    }

    // Reading the width, height and maximum value contained in the header of the image
    fscanf(inputImage, "%d %d %d", &width, &height, &maxVal);

    // Calling the saveImage function to obtain the matrix containint image data to be elaborated
    imageMatrix = saveImage(inputImage, width, height);

    // Closing the file
    fclose(inputImage);

    int treshold=atoi(argv[2]);
    
    binThreshold(imageMatrix, treshold, width, height);

    binComplement(imageMatrix, width, height);

    binOpening(imageMatrix, width, height);
    
    identifyBorders(imageMatrix, width, height);

    for(int i=0; i<height; i++){
        for(int j=0; j<width; j++){
             imageMatrix[i][j] = imageMatrix[i][j]*255;
        }
    }

    writeImage(imageMatrix, width, height, maxVal, "./out/borders.pgm");

 return 0;

}



/*

 gcc serial.c -o serial

 ./serial coins.pgm 80

*/
