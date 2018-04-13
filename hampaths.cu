#include <iostream>
#include <stdio.h>
#include <stdbool.h> /*helps with bool data type*/
#include <string.h> /* memset */
#include <unistd.h> /* close */
#include <emmintrin.h>
#include <sys/time.h> /*allows system type*/





struct timeval start, end;
void starttime() {
  gettimeofday( &start, 0 );
}

void endtime(const char* c) {
   gettimeofday( &end, 0 );
   double elapsed = ( end.tv_sec - start.tv_sec ) * 1000.0 + ( end.tv_usec - start.tv_usec ) / 1000.0;
   printf("%s: %f ms\n", c, elapsed); 
}

//initializes matrix values to 0,
//AS OF RIGHT NOW , ONLY USES THREADS!!!
__global__ void initTo(int *matx, int vs, int set){
	int index = threadIdx.x;
  	int stride = blockDim.x;
	for(int i = index ; i < vs; i += stride ){
		matx[i] = set;
	}

}
//initializes truth table to true
////AS OF RIGHT NOW , ONLY USES THREADS!!!
__global__ void initTo(bool *t_table, int vs, bool set){
	int index = threadIdx.x;
  	int stride = blockDim.x;
	for(int i = index ; i < vs; i += stride ){
		t_table[i] = set;
	}

}

__global__ void checkPathsTH(int *matrix, int vertices, int pathsize, int *paths, bool *t_table){
	//int index = threadIdx.x;
  	int stride = blockDim.x;
	int index = blockIdx.x*blockDim.x + threadIdx.x;

	//printf("%i\n",stride);
	//index is the actual index!!!
	//so 0-95
	//stops trash from going out of bounds
/*
	if(index<(pathsize)){
		for(int j = 0;j<vertices;j++){
			printf("%i's own %i called %i ",index,j,paths[(index*vertices)+j]);
			//THIS GETS YOU WHAT YOU NEED paths[(index*vertices)+j]
		}
		//printf("%i, %i, %i\n",index,stride,paths[index]);
	}
*/


	for(int i = index ; i < pathsize; i += stride){

		for(int j = 0;j<vertices;j++){

			if(j==0){
				if(matrix[(paths[(index*vertices)+j] * vertices) + paths[(index*vertices)+j+1]] != 1 &&
					matrix[(paths[(index*vertices)+j+1] * vertices) + paths[(index*vertices)+j]] != 1 )
				{
					t_table[index] = false;
					break;
				}
			}else if(j > 0 && j < vertices-1){
				if( (matrix[(paths[(index*vertices)+j] * vertices) + paths[(index*vertices)+j+1]] != 1 &&
					matrix[(paths[(index*vertices)+j+1] * vertices) + paths[(index*vertices)+j]] != 1 )
					||
				    (matrix[(paths[(index*vertices)+j] * vertices) + paths[(index*vertices)+j-1]] != 1 &&
					matrix[(paths[(index*vertices)+j-1] * vertices) + paths[(index*vertices)+j]] != 1 )
					)
					
				{
					t_table[index] = false;
					break;
				}

			}else if(j== vertices -1){
				if(matrix[(paths[(index*vertices)+j] * vertices) + paths[(index*vertices)+j-1]] != 1 &&
					matrix[(paths[(index*vertices)+j-1] * vertices) + paths[(index*vertices)+j]] != 1 )
				{
					t_table[index] = false;
					break;
				}


			}
			
		}
	}	

	//__syncthreads();
}


void checkPathLocal(int *ptr_path,int vertices,int pathsize, int *matrix, bool *possiblePaths ){
	for(int outcount = 0; outcount < pathsize;outcount++){
        for (int i = 0; i < vertices; i++){

            //the first one only checks the next one, checks if the path is even possible
            if(i == 0){
                if(matrix[(ptr_path[(outcount*vertices)+i] * vertices) + ptr_path[(outcount*vertices)+i+1]] != 1 &&
                   matrix[(ptr_path[(outcount*vertices)+i+1] * vertices) + ptr_path[(outcount*vertices)+i]] != 1)
                {
                    possiblePaths[outcount] = false;
                    break;
                }
            }
            else if(i > 0 && i < vertices - 1){
                if( (matrix[(ptr_path[(outcount*vertices)+i] * vertices) + ptr_path[(outcount*vertices)+i+1]] != 1 &&
                     matrix[(ptr_path[(outcount*vertices)+i+1] * vertices) + ptr_path[(outcount*vertices)+i]] != 1) 
                    ||
                    (matrix[(ptr_path[(outcount*vertices)+i] * vertices) + ptr_path[(outcount*vertices)+i-1]] != 1  &&
                     matrix[(ptr_path[(outcount*vertices)+i-1] * vertices) + ptr_path[(outcount*vertices)+i]] != 1)
                    )
                {
                    possiblePaths[outcount] = false;
                    break;
                }
            }
            else if(i == vertices - 1){
                if( matrix[(ptr_path[(outcount*vertices)+i] * vertices) + ptr_path[(outcount*vertices)+i-1]] != 1  &&
                    matrix[(ptr_path[(outcount*vertices)+i-1] * vertices) + ptr_path[(outcount*vertices)+i]] != 1)
                {
                    possiblePaths[outcount] = false;
                    break;
                } 
            }
        }   
    }
	endtime("NORMAL");
}




//just makes unweighted edges in a 1D array
void makeEdge(int to, int from, int *num, int vertices) {
        *(num + (((to) * (vertices)) + (from))) = 1;
	*(num + (((from) * (vertices)) + (to))) = 1;
	//matrix[((to) * (vertices)) + (from)] = 1;
        //matrix[((from) * (vertices)) + (to)] = 1;
        //matrix[from][to] = 1;
}

//swaps values, assists permute
    int swap(int a[], int i, int j)
    {
        int temp = a[i];
        a[i] = a[j];
        a[j] = temp;
        return a[j];
    }


//places in array
void getArray(int ha[] , int *ptr_path , int vertices, int pathsize){
        //finds the first path that is empty and places the array there!!!
        for(int i = 0; i < pathsize ;i++){
            //finds the first one that wasn't touched
            if( *(ptr_path + (i*vertices))  == -1){
                //maybe there could be a way to binary search this to make it faster!
                //i dunno :3
                //paths[i] = ha;
                for(int j =0; j < vertices ; j++){
                    //write to the original array
                    *(ptr_path + (i*vertices)+j) = ha[j];
                }
                //stops when found
                return;
            }
        }
}


void permute(int str[] , int l, int r, int *ptr_path, int vertices, int pathsize){
        //printf("Value of paths = %p\n",ptr_path);
        if (l == r){
            //array created to save each combination
            int raw [vertices]; 
            //prints the combinations and save them to the array.
            for(int j=0; j < vertices ;j++){
                //printf("%i", str[j]);
                raw[j] = str[j];
            }
           getArray(raw,ptr_path,vertices,pathsize);
           //printf("+++\n");
            //count = 0;
        }
        else
        {
            for (int i = l; i <= r; i++)
            {
                str[i] = swap(str,l,i);
                //permute(str, l+1, r);
                permute(str, l+1, r,ptr_path, vertices, pathsize);
                str[i] = swap(str,l,i);    
            }
        }
        
    }





int main(void) {

	const int vertices = 5;

	int *matrix;
	//int matrix[vertices*vertices] ;
 	//memset(matrix, 0, sizeof matrix);
	
	//1D MATRIX ARRAY INNITIALLIZED IN 
	//ALLOCATES
	cudaMallocManaged(&matrix, (vertices*vertices)*sizeof(int));
	//initiates to 0
	initTo<<<1,100>>>(matrix,vertices*vertices,0);
	cudaDeviceSynchronize();


	int pathsize = 1;
    	for (int i = 1; i < vertices + 1; i++)
    	{
        	pathsize = pathsize * i;
    	}
    	printf("%i\n",pathsize);

	bool *possiblePaths;
	cudaMallocManaged(&possiblePaths, (pathsize)*sizeof(bool));

	//bool possiblePaths[pathsize];
	//initialize to true
    	//memset(possiblePaths, true, sizeof possiblePaths);
	initTo<<<1,100>>>(possiblePaths,pathsize,true);
	cudaDeviceSynchronize();
	
	int *allPaths;
	cudaMallocManaged(&allPaths, (pathsize*vertices)*sizeof(int));
	initTo<<<1,100>>>(allPaths,pathsize*vertices,-1);
	cudaDeviceSynchronize();



//creating the graph 
	// ---
	//| / |
	// ---

	
	makeEdge(0,2,matrix,vertices);
	makeEdge(0,3,matrix,vertices);
   	makeEdge(1,4,matrix,vertices);
    	makeEdge(1,3,matrix,vertices);
	makeEdge(2,4,matrix,vertices);

	printf("The adjacency matrix for the given graph is: ");
   	 printf("\n  ");
    	for (int i = 0; i < vertices; i++)
       		printf("%i ",i+1);
    
    	for (int j = 0; j < vertices*vertices; j++) {
        	if(j % vertices == 0)
                	printf("\n%i ",(j/4)+1);
            printf("%i ", matrix[j]);
    	}
        printf("\n");

printf("%i\n",allPaths[0]);
printf("%i\n",allPaths[95]);	

int arr[vertices];
	for(int pattern = 0 ; pattern < vertices; pattern++){
        	arr[pattern] = pattern;
}

permute(arr,0,vertices - 1,allPaths,vertices,pathsize);

//printf("%i\n",allPaths[0]);
//printf("%i\n",allPaths[95]);
for (int i = 0; i < pathsize; i++)
        {
            if(i<10)
            printf(" %i |",i);
            else
                printf("%i |",i);

            for (int j = 0; j < vertices; j++)
            {
                printf("%d ", allPaths[(i*vertices)+j]);
            }
            //ptr_path++;
            printf("\n");
		
        }
/*
printf("%i\n",allPaths[92]);
printf("%i\n",allPaths[93]);
printf("%i\n",allPaths[94]);
printf("%i\n",allPaths[95]);
*/


starttime();
checkPathLocal(allPaths,vertices,pathsize,matrix,possiblePaths);
cudaDeviceSynchronize();

printf("POSSIBLE PATHS\n");
for(int wow = 0 ; wow<pathsize;wow++){
        //printf("");
        if(possiblePaths[wow] == true){
            for(int innerwow = 0 ; innerwow < vertices ; innerwow++){
                printf("%i", allPaths[(wow*vertices)+innerwow]);
            }
            printf("\n");
        }
        //printf("%d ",possiblePaths[wow]);
    }
printf("\n");

cudaDeviceSynchronize();
initTo<<<1,100>>>(possiblePaths,pathsize,true);
cudaDeviceSynchronize();

cudaDeviceSynchronize();
starttime();
checkPathsTH<<<2,75>>>(matrix,vertices,pathsize,allPaths,possiblePaths);
endtime("GPU THREADS");
cudaDeviceSynchronize();

cudaDeviceSynchronize();
printf("POSSIBLE PATHS\n");
for(int wow = 0 ; wow<pathsize;wow++){
        //printf("");
        if(possiblePaths[wow] == true){
            for(int innerwow = 0 ; innerwow < vertices ; innerwow++){
                printf("%i", allPaths[(wow*vertices)+innerwow]);
            }
            printf("\n");
        }
        //printf("%d ",possiblePaths[wow]);
    }
printf("\n");




	 // Free memory
  	cudaFree(matrix);
	cudaFree(possiblePaths);
	cudaFree(allPaths);
	//free(matrix);

	return 0;
}