#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <errno.h>
#include <mpi.h>

#define NAME_LEN 255

// Simple wrapper around a dense matrix. Can access elements either
// linearly through mat->all[i] or via 2d index through
// mat->data[i][j]
typedef struct {
  double *all;                  // Array of all elements, allows for easy sharing of entire matrix
  double **data;                // Pointers to individual rows to allow mat->data[r][c] access
  int nrows, ncols;             // Number of rows and columns
  int nnz;                      // Number of nonzeros in the matrix initially 
} densemat;
// Free a densemat
void free_densemat(densemat *mat){
  free(mat->all);
  free(mat->data);
  free(mat);
}
// Allocate a dense matrix of the given size. Initialize its elements
// to 0.0
densemat * allocate_densemat(int nrows, int ncols){
  densemat * mat = malloc(sizeof(densemat));
  mat->nrows = nrows;
  mat->ncols = ncols;
  mat->nnz = 0;
  int i,j;
  mat->all = malloc(nrows*ncols * sizeof(double *));
  for(i=0; i<nrows*ncols; i++){
    mat->all[i] = 0.0;
  }
  mat->data = malloc(nrows * sizeof(double**));
  for(i=0; i<nrows; i++){
    mat->data[i] = &mat->all[i*ncols];
  }
  return mat;
}    

// Read a dense matrix from the given file. Assumes a square matrix in
// the file.  Contents of the file are only row and column pairs whose
// value is assumed 1.0.  The format of the file starts with the
// number of rows and nonzeros (number of lines in file) on the first
// line. Each subsequent line has a row/col entry whose value is
// assumed to be 1.0
densemat * load_row_col_as_densemat(char *fname){
  FILE *f = fopen(fname,"r");
  if(f==NULL){
    perror(fname);
    exit(1);
  }
  int nrows, nnz;
  fscanf(f, "%d %d", &nrows, &nnz);
  densemat *mat = allocate_densemat(nrows,nrows);
  int i;
  for(i=0; i<nnz; i++){
    int row,col;
    fscanf(f,"%d %d",&row,&col);
    if(row >= nrows || col >= nrows){
      fprintf(stderr,"ERROR: line %d has row/col %d %d for matrix with rows/cols %d %d\n",
              i+1,row,col,nrows,nrows);
      exit(1);
    }
    mat->data[row][col] = 1.0;
    mat->nnz++;
  }
  assert(i==nnz);
  fclose(f);
  return mat;
}
int main(int argc, char **argv){
  int npes, proc_id, name_len;
  char proc_name[NAME_LEN];
  MPI_Init (&argc, &argv);                      /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &proc_id);     /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &npes);        /* get number of processes */
  MPI_Get_processor_name(proc_name, &name_len); /* get the symbolic host name */

  if(argc < 3){
    printf("usage: %s row_col.txt damping\n  0.0 < damping <= 1.0\n",argv[0]);
    return -1;
  }
   
  int n; //NxN broadcast to processors
  int root_proc = 0;
  int r,c,i;
  int MAX_ITER = 10000;
  int iter;
  int elements_per_proc;
  int surplus;
  int *counts;
  int *displs;
  int *displsmod;
  int *countsmod;
  double TOL = 1e-3;
  double change = TOL*10;
  double cur_norm;
  double *cur_ranks;
  double *old_ranks;
  double *indiv_cur_ranks;
  densemat *mat;
  densemat *gathermat;
  densemat *scattermat;
  double diff;
  double damping_factor;
  double *colsums;
  double zelem;

  if(proc_id == root_proc){//things only proc0 needs to do
    damping_factor = atof(argv[2]);
    mat = load_row_col_as_densemat(argv[1]);
    n = mat->nrows; //store N for NxN so the other processors can know the amount of columns/rows
    printf("Loaded %s: %d rows, %d nonzeros\n",argv[1],mat->nrows,mat->nnz);
    // Allocate space for the page ranks and a second array to track
    // page ranks from the last iterations
    cur_ranks = malloc(mat->nrows * sizeof(double));
    old_ranks = malloc(mat->nrows * sizeof(double));

    //only proc0 needs to do this since we broadcast old ranks
    for(c=0; c<(mat->nrows); c++){
      cur_ranks[c] = 1.0 / mat->nrows;
      old_ranks[c] = cur_ranks[c];
    }
  }//end pro0

  MPI_Bcast(&n, 1, MPI_INT, root_proc, MPI_COMM_WORLD); //send the amount of NxN to everyone
  MPI_Bcast(&damping_factor, 1, MPI_DOUBLE, root_proc, MPI_COMM_WORLD); //send the amount of damping to everyone
  //malloc the correct amounts
  colsums = malloc(n*sizeof(double));
  counts = malloc(npes * sizeof(int));
  displs = malloc(npes * sizeof(int));
  countsmod = malloc(npes * sizeof(int));
  displsmod = malloc(npes * sizeof(int));
  elements_per_proc = n/npes;
  surplus = n % npes;

  //for uneven processors
  for(i=0; i<npes; i++){
      counts[i] = (i<surplus) ? elements_per_proc+1 : elements_per_proc;
      displs[i] = (i==0) ? 0 : displs[i-1]+counts[i-1];
  }
  //this made scatterv easier, modify all values by the number of cols so we can use mat->all
  for(i=0; i<npes; i++){
      countsmod[i] = counts[i] * n;
      displsmod[i] = displs[i] * n;
  }
  if(proc_id != root_proc){ //allocate space for oldranks all but root
    old_ranks = malloc(n * sizeof(double));
    mat = allocate_densemat(n,n);//all other processors needs to have atleast a blank square
  }
  scattermat = allocate_densemat(counts[proc_id],n);
  indiv_cur_ranks = malloc(scattermat->nrows * sizeof(double));

  MPI_Scatterv(mat->all, countsmod,displsmod, MPI_DOUBLE, //scatter the mat values from p0 to all
   	       scattermat->all, countsmod[proc_id], MPI_DOUBLE,
  	       root_proc, MPI_COMM_WORLD);

  //start normalize 
  for(c=0; c<n; c++){
      colsums[c] = 0.0;
    }
  for(r=0; r<counts[proc_id]; r++){
    for(c=0; c<n; c++){
	colsums[c] += scattermat->data[r][c];
    }
  }
  //need to make every processor have the colsums for the division step
  MPI_Allreduce(MPI_IN_PLACE, colsums, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for(r=0; r<counts[proc_id]; r++){
    for(c=0; c<n; c++){
      scattermat->data[r][c] /= colsums[c];
    }
  }
  //end normalize

  //apply damping
  zelem = (1.0-damping_factor) / n;
  for(r=0; r<counts[proc_id]; r++){
    for(c=0; c<n; c++){
      if(scattermat->data[r][c] != 0.0){ // Scale down nonzero entries
	scattermat->data[r][c] *= damping_factor;
      }
      scattermat->data[r][c] += zelem; // Every entry increases a little
    }
  }//end damping

  if(proc_id == root_proc){
  printf("Beginning Computation\n\n%4s %8s %8s\n","ITER","DIFF","NORM");
  }
  for(iter=1; change > TOL && iter<=MAX_ITER; iter++){
    // old_ranks assigned to cur_ranks
    if(proc_id == root_proc){
      for(c=0; c<mat->ncols; c++){//sets old ranks to current ranks
	old_ranks[c] = cur_ranks[c];
      }
    }
    MPI_Bcast(old_ranks, n, MPI_DOUBLE, root_proc, MPI_COMM_WORLD);//broadcast old ranks
    change = 0.0;
    cur_norm = 0.0;

    // Compute matrix-vector product: cur_ranks = Matrix * old_ranks 
    for(r=0; r<scattermat->nrows; r++){
      // Dot product of matrix row with old_ranks column
      indiv_cur_ranks[r] = 0.0;
      for(c=0; c<scattermat->ncols; c++){
	indiv_cur_ranks[r] += scattermat->data[r][c] * old_ranks[c];
      }
    }
    //gather the smaller chunks back into current ranks on p0
    MPI_Gatherv(indiv_cur_ranks, counts[proc_id], MPI_DOUBLE,
    	       cur_ranks, counts, displs, MPI_DOUBLE,
    	       root_proc, MPI_COMM_WORLD);
  
    if(proc_id == root_proc){
      for(c=0; c<n; c++){
       	diff = cur_ranks[c] - old_ranks[c]; //compute the difference
       	change += diff>0 ? diff : -diff;
       	cur_norm += cur_ranks[c]; // Tracked to detect any errors
       }
      printf("%3d: %8.2e %8.2e\n",iter,change,cur_norm);
    }   
    MPI_Bcast(&change, 1, MPI_DOUBLE, root_proc, MPI_COMM_WORLD); //broadcast the loop condition
  }

  if(proc_id == root_proc){//proc0 printing
    if(change < TOL){
      printf("CONVERGED\n");
    }
    else{
      printf("MAX ITERATION REACHED\n");
    }
    printf("\nPAGE RANKS\n");
    for(r=0; r<mat->nrows; r++){
      printf("%.8f\n",cur_ranks[r]);
    }
  }//end proc0

  //free the structures
   if(proc_id == root_proc){
     free(cur_ranks);//only exists on p0
   }
   free(old_ranks);
   free(indiv_cur_ranks);
   free_densemat(mat);
   free_densemat(scattermat);
   free(counts);
   free(countsmod);
   free(displs);
   free(displsmod);
   free(colsums);
    
  MPI_Finalize();
  return 0;
}
