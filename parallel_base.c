#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int size, my_rank;
MPI_Status status;


int* matrixToArray(int rows, int cols, int** matrix) {
    int* array = (int*)malloc(rows * cols * sizeof(int));
    if (array == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    int index = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            array[index++] = matrix[i][j];
        }
    }

    return array;
}

int** arrayToMatrix(int rows, int cols, int* array) {
    int** matrix = (int**)malloc(rows * sizeof(int*));
    if (matrix == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < rows; i++) {
        matrix[i] = (int*)malloc(cols * sizeof(int));
        if (matrix[i] == NULL) {
            printf("Memory allocation failed.\n");
            exit(1);
        }
    }

    int index = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = array[index++];
        }
    }

    return matrix;
}


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

int **imageToMatrix(FILE *inputImage, int width, int height, int newRows)
{
    int **matrix = (int **)malloc(height * sizeof(int *));
    for (int i = 0; i < height; i++)
    {
        matrix[i] = (int *)malloc(width * sizeof(int));
    }

    for (int i = 0; i < height-newRows; i++)
    {
        for (int j = 0; j < width; j++)
        {
            matrix[i][j]=fgetc(inputImage);
        }
    } 

    for (int i = height-newRows; i < height; i++) {
        for (int j = 0; j < width; j++) {
                matrix[i][j]= 0;
            }
        }
    
    fclose(inputImage);
    return matrix;

}


void binComplement(int rows, int cols, int **imageMatrix){

    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
                imageMatrix[i][j] = (imageMatrix[i][j] -1) * (-1);
        }
    }
}

void binErosion(int rows, int cols, int **eroded)
{

    int **image = (int **)malloc(rows * sizeof(int *));
        for (int i = 0; i < rows; i++)
        {
            image[i] = (int *)malloc(cols * sizeof(int));
        }

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            image[i][j] = eroded[i][j];
        }
    }

    for (int i = 1; i < rows - 1; i++)
    {
        for (int j = 1; j < cols - 1; j++)
        {

            if (image[i][j] * image[i + 1][j] * image[i + 1][j + 1] *
                    image[i][j + 1] * image[i - 1][j + 1] * image[i - 1][j] *
                    image[i - 1][j - 1] * image[i][j - 1] * image[i + 1][j - 1] !=0)
                ;

            else
            {
                eroded[i][j] = 0;
                eroded[i + 1][j] = 0;
                eroded[i + 1][j + 1] = 0;
                eroded[i][j + 1] = 0;
                eroded[i - 1][j + 1] = 0;
                eroded[i - 1][j] = 0;
                eroded[i - 1][j - 1] = 0;
                eroded[i][j - 1] = 0;
                eroded[i + 1][j - 1] = 0;
            }
        }
    }
}

void binDilation(int rows, int cols, int **dilated)
{

        int **image = (int **)malloc(rows * sizeof(int *));
        for (int i = 0; i < rows; i++)
        {
            image[i] = (int *)malloc(cols * sizeof(int));
        }
    

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            image[i][j] = dilated[i][j];
        }
    }

    for (int i = 1; i < rows - 1; i++)
    {
        for (int j = 1; j < cols - 1; j++)
        {

            if ((image[i][j] + image[i + 1][j] + image[i + 1][j + 1] +
                 image[i][j + 1] + image[i - 1][j + 1] + image[i - 1][j] +
                 image[i - 1][j - 1] + image[i][j - 1] + image[i + 1][j - 1]) == 0)
                ;

            else
            {
                dilated[i][j] = 1;
                dilated[i + 1][j] = 1;
                dilated[i + 1][j + 1] = 1;
                dilated[i][j + 1] = 1;
                dilated[i - 1][j + 1] = 1;
                dilated[i - 1][j] = 1;
                dilated[i - 1][j - 1] = 1;
                dilated[i][j - 1] = 1;
                dilated[i + 1][j - 1] = 1;
            }
        }
    }
}

void binOpening(int rows, int cols, int **imageMatrix)
{

    binErosion(rows, cols, imageMatrix);

    binDilation(rows, cols, imageMatrix);
}

void identifyBorders(int rows, int cols, int **imageMatrix)
{

        int **image = (int **)malloc(rows * sizeof(int *));
        for (int i = 0; i < rows; i++)
        {
            image[i] = (int *)malloc(cols * sizeof(int));
        }

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            image[i][j] = imageMatrix[i][j];
        }
    }

    binErosion(rows, cols, imageMatrix);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            imageMatrix[i][j] = (imageMatrix[i][j] - image[i][j]) * (-1);
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

        int treshold=atoi(argv[2]); 

    newHeight = (height+size-1)/size;
    newRows = newHeight * size - height;  
    height= height + newRows;

    int **imageMatrix= imageToMatrix(inputImage, width, height, newRows);

    MPI_Barrier(MPI_COMM_WORLD);

    int *result = (int*)malloc(height * width * sizeof(int));

    printf("\n PRE Scatter, rank: %d  and value: %d\n\n", my_rank, imageMatrix[0][0]);

    int *sendArray = (int*)malloc(height * width * sizeof(int));
    sendArray = matrixToArray(height, width, imageMatrix);

    int *recvArray = (int*)malloc(height/size * width * sizeof(int));

    MPI_Scatter(sendArray, ((height/size)*width), MPI_INT, 
                recvArray, ((height/size)*width), MPI_INT,
                0, MPI_COMM_WORLD);

    int **recvMatrix = arrayToMatrix(height/size, width, recvArray);
  
    printf("\n POST Scatter, rank: %d  and value: %d %d\n\n", my_rank, recvMatrix[0][0], imageMatrix[0][0]);

     if (my_rank==0){

              writeImage(height/size, width, maxVal, recvMatrix, "./out/original_0.pgm");
    }

        if (my_rank==1){

              writeImage(height/size, width, maxVal, recvMatrix, "./out/original_1.pgm");
    }

        if (my_rank==2){

              writeImage(height/size, width, maxVal, recvMatrix, "./out/original_2.pgm");
    }

    binThreshold(height/size, width, recvMatrix, treshold);
    
    binComplement(height/size, width, recvMatrix);
   
    binOpening(height/size, width, recvMatrix);

    identifyBorders(height/size, width, recvMatrix);
 
      for(int i=0; i<height/size; i++){
        for(int j=0; j<width; j++){
             recvMatrix[i][j] = recvMatrix[i][j]*255;
        }
    }

         if (my_rank==0){

              writeImage(height/size, width, maxVal, recvMatrix, "./out/borders_0.pgm");
    }

        if (my_rank==1){

              writeImage(height/size, width, maxVal, recvMatrix, "./out/borders_1.pgm");
    }

        if (my_rank==2){

              writeImage(height/size, width, maxVal, recvMatrix, "./out/borders_2.pgm");
    }

    recvArray = matrixToArray(height/size, width, recvMatrix);

    MPI_Gather(recvArray, height*width/size, MPI_INT, result, height*width/size, MPI_INT, 0, MPI_COMM_WORLD);

    int **resultMatrix = arrayToMatrix(height, width, result);

    printf("\n post Gather, rank: %d \n", my_rank);

    for (int i = 0; i < height / size; i++) {
                free(recvMatrix[i]);
    }
    free(recvMatrix);

     
     if (my_rank==0){

        writeImage(height, width, maxVal, resultMatrix, "./out/borders.pgm");
        
        stop_time = MPI_Wtime();
        tot_time= stop_time - start_time;
        printf("Total time: %f", tot_time);

            for (int i = 0; i < height; i++) {
                free(imageMatrix[i]);
            }
        free(imageMatrix);
    }

    MPI_Finalize();
    return 0;

}



/*

 mpicc parallel_base.c -o parallel_base

 mpirun -n 3 parallel_base coins.pgm 80

*/
