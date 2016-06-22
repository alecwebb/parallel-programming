#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define NAME_LEN 255

double calc_next(double posleft, double pos, double posright){
  double left_diff, right_diff, delta;
  double k = 0.5; //thermal conductivity constant
  left_diff = pos - posleft;
  right_diff = pos - posright;
  delta = -k*(left_diff + right_diff);
  return (pos + delta);
}
int main(int argc, char **argv){
  int npes, proc_id, name_len;
  char proc_name[NAME_LEN];
  MPI_Init (&argc, &argv);                      /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &proc_id);     /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &npes);        /* get number of processes */
  MPI_Get_processor_name(proc_name, &name_len); /* get the symbolic host name */

  if(argc < 4){
    printf("usage: %s max_time width print\n max_time: int\n width: int\n print: 1 print output, 0 no printing\n",
	    argv[0]);
    return 0;
  }
  
  int max_time = atoi(argv[1]); // Number of time steps to simulate
  int width = atoi(argv[2]);    // Number of cells in the rod
  int print = atoi(argv[3]);    // print option
  double initial_temp = 50.0;   // Initial temp of internal cells 
  double L_bound_temp = 20.0;   // Constant temp at Left end of rod
  double R_bound_temp = 10.0;   // Constant temp at Right end of rod
  double **H;                   // 2D array of temps at times/locations 
  double **root_data;           // To store the final data on proc0
  double left_val, right_val;   // used for communication of left and right node values
  int rootproc = 0;             // 0 is the root processor
  int indiv_width = width/npes; // To determine how many columns each extra processor gets
  int internal = indiv_width-2; // to determine how many cols dont deal with edge data
  int t,p;
  
  H = malloc(sizeof(double*)*max_time); 
  for(t=0;t<max_time;t++){
     H[t] = malloc(sizeof(double*)*indiv_width);
  }
  t = 0;
  for(p=0; p<indiv_width; p++){//we dont care about last columns, deal with that in calculation step
    H[t][p] = initial_temp;    //initialize to initial temperature
  }
  
  // Simulate the temperature changes for internal cells
  for(t=0; t<max_time-1; t++){
    //ESTABLISH COMMUNICATION AND TRANSFERS
    if(npes > 1){//communication only happens when we have more than one processor
      if(proc_id == rootproc){//even root
	MPI_Send(&H[t][indiv_width-1], 1, MPI_DOUBLE, proc_id+1, 1, MPI_COMM_WORLD);
	MPI_Recv(&right_val, 1, MPI_DOUBLE, proc_id+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      else if(proc_id == npes-1 && proc_id % 2 == 0){//even last
	MPI_Recv(&left_val, 1, MPI_DOUBLE, proc_id-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Send(&H[t][0], 1, MPI_DOUBLE, proc_id-1, 1, MPI_COMM_WORLD);
      }
      else if(proc_id == npes-1 && proc_id % 2 == 1){//odd last
	MPI_Recv(&left_val, 1, MPI_DOUBLE, proc_id-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Send(&H[t][0], 1, MPI_DOUBLE, proc_id-1, 1, MPI_COMM_WORLD);
      }
      else if(proc_id % 2 == 0){//even middle
	MPI_Send(&H[t][indiv_width-1], 1, MPI_DOUBLE, proc_id+1, 1, MPI_COMM_WORLD);
	MPI_Recv(&right_val, 1, MPI_DOUBLE, proc_id+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(&left_val, 1, MPI_DOUBLE, proc_id-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Send(&H[t][0], 1, MPI_DOUBLE, proc_id-1, 1, MPI_COMM_WORLD);
      }
      else if(proc_id % 2 == 1){//odd middle
	MPI_Recv(&left_val, 1, MPI_DOUBLE, proc_id-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Send(&H[t][0], 1, MPI_DOUBLE, proc_id-1, 1, MPI_COMM_WORLD);
	MPI_Send(&H[t][indiv_width-1], 1, MPI_DOUBLE, proc_id+1, 1, MPI_COMM_WORLD);
	MPI_Recv(&right_val, 1, MPI_DOUBLE, proc_id+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      else{ printf("couldnt find a processor for communications"); }
    }
 
    //PERFORM THE NECESSARY CALCULATIONS    
    if(proc_id == rootproc){//do root
      H[t][0] = L_bound_temp; //fill in left static col
      if(npes < 2){ //in the case we only had one processor working
	H[t][indiv_width-1] = R_bound_temp; //fill in right static col
      }
      //internal is up to indiv_width-2, so 3 cols then 3-2 = 1 internal col
      for(p=1; p<=internal; p++){//0 is the static temp col
	H[t+1][p] = calc_next( H[t][p-1], H[t][p], H[t][p+1] );
      }      
      H[t+1][indiv_width-1] = calc_next( H[t][indiv_width-2], H[t][indiv_width-1], right_val );//handle the right most column using communicated data
      if(t==max_time-2){//set the last static right col (bottom left)
	H[t+1][0] = L_bound_temp;
	if(npes < 2){ //in the case we only had one processor working (static col bottom right)
	  H[t+1][indiv_width-1] = R_bound_temp;
	}
      }
    }
    else if((proc_id > 0) && (proc_id < npes-1)){//do internal proc calculations
      H[t+1][0] = calc_next( left_val, H[t][0] , H[t][1]);//handle value communicated from the left  
      for(p=1; p<indiv_width-1; p++){//handle middle values requiring no externals
	H[t+1][p] = calc_next( H[t][p-1], H[t][p], H[t][p+1] );
      }
      H[t+1][indiv_width-1] = calc_next( H[t][indiv_width-2], H[t][indiv_width-1], right_val );//handle value communicated from the right
    }
    else if(proc_id == npes-1 && npes > 1){//do last, this is only needed when there are multiple processors.
      H[t][indiv_width-1] = R_bound_temp;//set last static column
      H[t+1][0] = calc_next( left_val, H[t][0], H[t][1] );//handle value communicated from the left
      for(p=1; p<=internal; p++){//handle values up to static column
	H[t+1][p] = calc_next( H[t][p-1], H[t][p], H[t][p+1] );
      } 
      if(t==max_time-2){//set the last static right col (bottom right)
	H[t+1][indiv_width-1] = R_bound_temp;
      }
    }
    else{ printf("couldnt find a processor for calculations"); }
  }
  //begin gather procedure
  if(proc_id == rootproc){//make space for root_data array
    root_data = malloc(sizeof(double*)*max_time);
    for(t=0;t<max_time;t++){
       root_data[t] = malloc(sizeof(double*)*width);
    }
  }
  for(t=0; t<max_time; t++){//gather the data to root
     MPI_Gather(H[t], indiv_width, MPI_DOUBLE,
		root_data[t],indiv_width, MPI_DOUBLE,
		rootproc, MPI_COMM_WORLD);
  }
  //end gather procedure
  free(H);//everyone free H
  if(print == 1){
    if(proc_id == rootproc){//start proc0 printing
      // Print results
      printf("Temperature results for 1D rod\n");
      printf("Time step increases going down rows\n");
      printf("Position on rod changes going accross columns\n");
      // Column headers
      printf("%3s| ","");
      for(p=0; p<width; p++){
	printf("%5d ",p);
      }
      printf("\n");
      printf("%3s+-","---");
      for(p=0; p<width; p++){
	printf("------");
      }
      printf("\n");
      // Row headers and data
      for(t=0; t<max_time; t++){
	printf("%3d| ",t);
	for(p=0; p<width; p++){
	  // printf("%5.1f ",H[t][p]);
	  printf("%5.1f ",root_data[t][p]);
	}
	printf("\n");
      }
      free(root_data);//free the root_data array
    }//end proc0 printing
  }
  MPI_Finalize();
  return 0;
}
      
