#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int size, my_rank;
MPI_Status status;

void binThreshold(int rows, int cols, int **imageMatrix, int value){

    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){

            if (imageMatrix[i][j] < value)
                imageMatrix[i][j] = 1;
            else
                imageMatrix[i][j] = 0;
        }
    }
}

int **imageToMatrix(FILE *inputImage, int width, int height)
{
    int **matrix = (int **)malloc(height * sizeof(int *));
    for (int i = 0; i < height; i++)
    {
        matrix[i] = (int *)malloc(width * sizeof(int));
    }

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            matrix[i][j]=fgetc(inputImage);
        }
    }
    
    return matrix;
}


void binComplement(int rows, int cols, int **imageMatrix){

    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
                imageMatrix[i][j] = (imageMatrix[i][j] -1) * (-1);
        }
    }
}

void binErosion(int rows, int cols, int **eroded){

        
        int **image = (int **)malloc(rows * sizeof(int *));
        for (int i = 0; i < rows; i++)
        {
            image[i] = (int *)malloc(cols * sizeof(int));
        }

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


void binDilation(int rows, int cols, int **dilated){


        int **image = (int **)malloc(rows * sizeof(int *));
        for (int i = 0; i < rows; i++)
        {
            image[i] = (int *)malloc(cols * sizeof(int));
        }

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


void binOpening(int rows, int cols, int **imageMatrix){

    binErosion(rows, cols, imageMatrix);

    binDilation(rows, cols, imageMatrix);

}


void identifyBorders(int rows, int cols, int **imageMatrix) {

        int **image = (int **)malloc(rows * sizeof(int *));
        for (int i = 0; i < rows; i++)
        {
            image[i] = (int *)malloc(cols * sizeof(int));
        }

        for (int i=0; i<rows; i++){
            for (int j=0; j<cols; j++){
                image[i][j]= imageMatrix[i][j];
            }
        }

    binErosion(rows, cols, imageMatrix);

    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
                imageMatrix[i][j] = image[i][j] - imageMatrix[i][j];
        }
    }

}


void writeImage(int rows, int cols, int maxVal, int **matrix, const char *outputFileName) {
    FILE *outputImage = fopen(outputFileName, "wb");
    if (outputImage == NULL) {
        printf("Error: Could not open %s for writing.\n", outputFileName);
        exit(1);
    }

    fprintf(outputImage, "P2\n");
    fprintf(outputImage, "%d %d\n", cols, rows);
    fprintf(outputImage, "%d\n", maxVal); 

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            fprintf(outputImage, "%d ", matrix[i][j]);
        }
        fprintf(outputImage, "\n");
    }

    fclose(outputImage);
}
    
  
int main(int argc, char *argv[]) {

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

    start_time = MPI_Wtime();

    if (argc != 3) {
        fprintf(stderr, "Error: Wrong number of arguments\nUsage: mpirun -n <size> input_image.pgm <threshold>\n");
        return 1;
    }
  
    inputImage = fopen(argv[1], "rb");
    if (inputImage == NULL) {
        printf("Error: %s can't be opened\n", argv[1]);
        return 1;
    }

    fscanf(inputImage, "%s", magicNumber);


    if (magicNumber[0] != 'P' || magicNumber[1] != '5') {
        printf("Error: file format not supported\n");
        return 1;
    }

    fscanf(inputImage, "%d %d %d", &width, &height, &maxVal);

    newHeight = (height+size-1)/size;
    newRows = newHeight * size - height;

    }   

    MPI_Bcast(&height,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&newRows,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&width,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&maxVal,1,MPI_INT,0,MPI_COMM_WORLD);

    height+=newRows;

    int treshold=atoi(argv[2]);

        int **recvMatrix = (int **)malloc(height/size * sizeof(int *));
        for (int i = 0; i < height/size; i++)
        {
            recvMatrix[i] = (int *)malloc(width * sizeof(int));
        }
	
    int **result = (int **)malloc(height * sizeof(int *));
        for (int i = 0; i < height; i++)
        {
            result[i] = (int *)malloc(width * sizeof(int));
        }

    int **imageMatrix;

if(my_rank==0){

	imageMatrix = imageToMatrix(inputImage, width, height);

    fclose(inputImage);

}

  MPI_Scatter(imageMatrix, height/size*width, MPI_INT, 
		    recvMatrix, height*width/size, MPI_INT,
                    0, MPI_COMM_WORLD);


    binThreshold(height/size, width, recvMatrix, treshold);
    binComplement(height/size, width, recvMatrix);

    binOpening(height/size, width, recvMatrix);

    identifyBorders(height/size, width, recvMatrix);
 

      for(int i=0; i<height/size; i++){
        for(int j=0; j<width; j++){
             recvMatrix[i][j] = recvMatrix[i][j]*255;
        }
    }

    MPI_Gather(recvMatrix, height*width/size, MPI_INT, result, height*width/size, MPI_INT, 0, MPI_COMM_WORLD);


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
