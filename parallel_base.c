#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int size, my_rank;
MPI_Status status;

void binThreshold(int rows, int cols, int imageMatrix[rows][cols], int value){

    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){

            if (imageMatrix[i][j] < value)
                imageMatrix[i][j] = 1;
            else
                imageMatrix[i][j] = 0;
        }
    }
}

void binComplement(int rows, int cols, int imageMatrix[rows][cols]){

    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
                imageMatrix[i][j] = (imageMatrix[i][j] -1) * (-1);
        }
    }
}

void binErosion(int rows, int cols, int eroded[rows][cols]){

        
        int image[rows][cols];

        for (int i=0; i<rows; i++){
            for (int j=0; j<cols; j++){
                image[i][j]= eroded[i][j];
            }

        }
        
            for(int i = 1; i < rows -1 ; i++){
            for(int j = 1; j < cols -1; j++){
		
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

}


void binDilation(int rows, int cols, int dilated[rows][cols]){


        int image[rows][cols];

        for (int i=0; i<rows; i++){
            for (int j=0; j<cols; j++){
                image[i][j]= dilated[i][j];
            }
        }

        for(int i = 1; i < rows -1 ; i++){
            for(int j = 1; j < cols -1; j++){
	
		if((image[i][j] + image[i+1][j] + image[i+1][j+1] +
		   image[i][j+1] + image[i-1][j+1] + image[i-1][j] +
		   image[i-1][j-1] + image[i][j-1] + image[i+1][j-1])==0);
		   
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

}


void binOpening(int rows, int cols, int imageMatrix[rows][cols]){

    binErosion(rows, cols, imageMatrix);

    binDilation(rows, cols, imageMatrix);

}


void identifyBorders(int rows, int cols, int imageMatrix[rows][cols]) {

    int image[rows][cols];

        for (int i=0; i<rows; i++){
            for (int j=0; j<cols; j++){
                image[i][j]= imageMatrix[i][j];
            }

        }

    binErosion(rows, cols, imageMatrix);

    // Calculate the difference between the original and eroded images
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
                imageMatrix[i][j] = image[i][j] - imageMatrix[i][j];
        }
    }

}


void writeImage(int rows, int cols, int maxVal, int matrix[rows][cols], const char *outputFileName) {
    FILE *outputImage = fopen(outputFileName, "wb");
    if (outputImage == NULL) {
        printf("Error: Could not open %s for writing.\n", outputFileName);
        exit(1);
    }

    // Write PGM header
    fprintf(outputImage, "P2\n");
    fprintf(outputImage, "%d %d\n", cols, rows);
    fprintf(outputImage, "%d\n", maxVal); 

    // Write pixel values
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            fprintf(outputImage, "%d ", matrix[i][j]);
        }
        fprintf(outputImage, "\n");
    }

    fclose(outputImage);
}
    
  
int main(int argc, char*argv[]) {

    FILE *inputImage;
    char magicNumber[3];
    int width, height, maxVal, newRows, newHeight;
    double tot_time = 0.0;
    double start_time = 0.0;
    double stop_time = 0.0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Status status;


    if(my_rank==0){
    /* The program requires as an input the name of the image file to be elaborated; if the number of provided arguments is 
     * strictly less than 2, an error message is printed 
     */

    start_time = MPI_Wtime();

    if (argc != 3) {
        fprintf(stderr, "Error: Wrong number of arguments\nUsage: mpirun -n <size> input_image.pgm <threshold>\n");
        return 1;
    }

    /* Opening the input image to be elaborated; if the fopen returns a NULL value, the opening of the file was not successful 
    *  and an error message is printed
    */
  
    inputImage = fopen(argv[1], "rb");
    if (inputImage == NULL) {
        printf("Error: %s can't be opened\n", argv[1]);
        return 1;
    }

    /* Reading the header of the input image; if the magic number of the header is different from P2 (magic number for the PGM 
     * file format), an error message is printed 
     */
    fscanf(inputImage, "%s", magicNumber);


    if (magicNumber[0] != 'P' || magicNumber[1] != '2') {
        printf("Error: file format not supported\n");
        return 1;
    }

    // Reading the width, height and maximum value contained in the header of the image
    fscanf(inputImage, "%d %d %d", &width, &height, &maxVal);

    // We are computing the extra rows because the MPI_Scatter gives the remaining rows to the processes in a randomic way
    // we add extra rows to let the division to remain in order so each submatrix does not have the randomic rows

    newHeight = (height+size-1)/size;
    newRows = newHeight * size - height;

    }   

    MPI_Bcast(&height,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&newRows,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&width,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&maxVal,1,MPI_INT,0,MPI_COMM_WORLD);

    height+=newRows;

    int treshold=atoi(argv[2]);

	int recvMatrix[height/size][width];
    int result[height][width];
    int imageMatrix[height][width];


if(my_rank==0){

        for (int i = 0; i < height-newRows; i++) {
            for (int j = 0; j < width; j++) {
                if (fscanf(inputImage, "%d", &imageMatrix[i][j]) != 1) {
                        // Handle error or unexpected end of file
                        fprintf(stderr, "Error reading pixel values from the file.\n");
                        exit(1);
                    }
            }
        }

        for (int i = height-newRows; i < height; i++) {
            for (int j = 0; j < width; j++) {
                imageMatrix[i][j]= 0;
            }
        }


     // Closing the file
    fclose(inputImage);

}

  	MPI_Scatter(imageMatrix, height/size*width, MPI_INT, 
		    recvMatrix, height*width/size, MPI_INT,
                    0, MPI_COMM_WORLD);


    binThreshold(height/size, width, recvMatrix, treshold);
    binComplement(height/size, width, recvMatrix);

    binOpening(height/size, width, recvMatrix);
 

      for(int i=0; i<height/size; i++){
        for(int j=0; j<width; j++){
             recvMatrix[i][j] = recvMatrix[i][j]*255;
        }
    }


    MPI_Gather(recvMatrix, height*width/size, MPI_INT, result, height*width/size, MPI_INT, 0, MPI_COMM_WORLD);


    if (my_rank==0){
        writeImage(height, width, maxVal, result, "./out/opened.pgm");
    }

    identifyBorders(height, width, result);


     if (my_rank==0){
        writeImage(height, width, maxVal, result, "./out/borders.pgm");
        stop_time = MPI_Wtime();
        tot_time= stop_time - start_time;
        printf("Total time: %f", tot_time);
    }


    MPI_Finalize();
    return 0;

}


/*

 mpicc parallel_base.c -o parallel_base

 mpirun -n 3 parallel_base coins.pgm 80


*/