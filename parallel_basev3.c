#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int size, my_rank;
MPI_Status status;

void binThreshold(int rows, int cols, int** imageMatrix, int value){

    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){

            if (imageMatrix[i][j] < value)
                imageMatrix[i][j] = 1;
            else
                imageMatrix[i][j] = 0;
        }
    }
}

void binComplement(int rows, int cols, int** imageMatrix){

    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
                imageMatrix[i][j] = (imageMatrix[i][j] -1) * (-1);
        }
    }
}

void binErosion(int rows, int cols, int** eroded){

        
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


void binDilation(int rows, int cols, int** dilated){


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


void binOpening(int rows, int cols, int** imageMatrix){

    binErosion(rows, cols, imageMatrix);

    binDilation(rows, cols, imageMatrix);

}


void identifyBorders(int rows, int cols, int** imageMatrix) {

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
                imageMatrix[i][j] = imageMatrix[i][j] - image[i][j];
        }
    }

}


void writeImage(int rows, int cols, int maxVal, int** matrix, const char *outputFileName) {
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

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Status status;

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

    /* Reading the header of the input image; if the magic number of the header is different from P5 (magic number for the PGM 
     * file format), an error message is printed 
     */
    fscanf(inputImage, "%s", magicNumber);


    if (magicNumber[0] != 'P' || magicNumber[1] != '2') {
        printf("Error: file format not supported\n");
        return 1;
    }

    // Reading the width, height and maximum value contained in the header of the image
    fscanf(inputImage, "%d %d %d", &width, &height, &maxVal);

    int treshold=atoi(argv[2]);


    int **result = (int **)malloc(height * sizeof(int*));

    for( int i = 0; i < height; i++) {
        result[i] = (int *)malloc(width * sizeof(int));
    }


    //int result[height][width];


    int imageMatrix[height][width];


        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                if (fscanf(inputImage, "%d", &imageMatrix[i][j]) != 1) {
                        // Handle error or unexpected end of file
                        fprintf(stderr, "Error reading pixel values from the file.\n");
                        exit(1);
                    }
            }
        }


     // Closing the file
    fclose(inputImage);

  // Dividing the rows of the input image matrix among processes
    int rowsPerProcess = width / size;
    int remainingRows = width % size; // Calculate the remaining rows
    int startRow = my_rank * rowsPerProcess + (my_rank < remainingRows ? my_rank : remainingRows);
    int endRow = startRow + rowsPerProcess + (my_rank < remainingRows ? 1 : 0); // Distribute remaining rows

    // Printing information about the range of rows each process is working on
    printf("Process %d working on rows %d to %d.\n", my_rank, startRow, endRow - 1);


    int **recvMatrix = (int **)malloc((endRow-startRow) * sizeof(int*));

    for( int i = 0; i < height; i++) {
        recvMatrix[i] = (int *)malloc(width * sizeof(int));
    }
  

    for (int j = 0; j < width; j++) {
        for (int i = 0; i < rowsPerProcess; i++) {  
            recvMatrix[i][j] = imageMatrix[i][j + startRow];      
        } 
    }


    binThreshold(endRow-startRow, width, recvMatrix, treshold);
    binComplement(endRow-startRow, width, recvMatrix);

    binOpening(endRow-startRow, width, recvMatrix);
 

      for(int i=0; i<endRow-startRow; i++){
        for(int j=0; j<width; j++){
             recvMatrix[i][j] = recvMatrix[i][j]*255;
        }
    }

      if (my_rank==0){
         writeImage(endRow-startRow, width, maxVal, recvMatrix, "./out/opened0.pgm");
   }

      if (my_rank==1){
         writeImage(endRow-startRow, width, maxVal, recvMatrix, "./out/opened1.pgm");
   }

      if (my_rank==2){
         writeImage(endRow-startRow, width, maxVal, recvMatrix, "./out/opened2.pgm");
   }

      if (my_rank==3){
         writeImage(endRow-startRow, width, maxVal, recvMatrix, "./out/opened3.pgm");
   } 

       if (my_rank==4){
         writeImage(endRow-startRow, width, maxVal, recvMatrix, "./out/opened4.pgm");
   }

      if (my_rank==5){
         writeImage(endRow-startRow, width, maxVal, recvMatrix, "./out/opened5.pgm");
   } 

      if (my_rank==6){
         writeImage(endRow-startRow, width, maxVal, recvMatrix, "./out/opened6.pgm");
   } 



    if (my_rank == 0) {
            
              for(int i=0; i<endRow-startRow; i++){
                for(int j=0; j<width; j++){
                    result[i][j] = recvMatrix[i][j];
                 }
                }

        int *processStartRow = (int *)malloc(sizeof(int));
        int *processEndRow = (int *)malloc(sizeof(int));
            // Master process receives transformed columns from other processes
            for (int process = 1; process < size; process++) {
                MPI_Recv((processStartRow), 1, MPI_INT, process, 2, MPI_COMM_WORLD, &status);
                MPI_Recv((processEndRow), 1, MPI_INT, process, 3, MPI_COMM_WORLD, &status);
                //printf("ricevo dal processo: %d le righe: %d - %d\n", process, processStartRow, processEndRow);
                for (int i = *processStartRow; i < *processEndRow; i++) {
                    MPI_Recv(recvMatrix[i], width, MPI_INT, process, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            }

        } else {
            MPI_Send(&startRow, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
            MPI_Send(&endRow, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
            // Slave processes send their transformed rows to the master process
            for (int i = startRow; i < endRow; i++) {
                    MPI_Send(recvMatrix[i], width, MPI_INT, 0, i, MPI_COMM_WORLD);
            
            }
        }

        for (int i = 0; i < height; i++){
            free(recvMatrix[i]);
        }
        free(recvMatrix);
    
    if (my_rank==0){
        writeImage(height, width, maxVal, result, "./out/opened.pgm");
    }

     identifyBorders(height, width, result);


     if (my_rank==0){
        writeImage(height, width, maxVal, result, "./out/borders.pgm");
    }


    for (int i = 0; i < height; i++){
            free(result[i]);
        }
    free(result);

    MPI_Finalize();
    return 0;

}



/*

 mpicc parallel_basev3.c -o parallel_basev3

 mpirun -n 3 parallel_basev3 coins.pgm 80


*/