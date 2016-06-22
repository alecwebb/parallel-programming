#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <crack.h>
#include <omp.h>
#include <pthread.h>

volatile int pfound = 0;//for use in pthreads to determine if a solution has been found

int try_crackomp(char *target,
              dict_t **dicts, int dicts_len, int dict_pos,
              char *buf, int buflen, int bufpos)
{
  int found = 0;
  #pragma omp parallel shared(found)
  {
  //printf("thread %d\n",omp_get_thread_num());
  char *buff = malloc(buflen *sizeof(char));
  dict_t *cur_dict = dicts[dict_pos];
    #pragma omp for
    for(int i=0; i<dict_get_word_count(cur_dict); i++){
      if(found)
        continue;

      // Append a new word to the end of buf
      char *word = dict_get_word(cur_dict, i);
      int bp = bufpos + strlen(word);
      if(bp > buflen){
        fprintf(stderr,"WARNING: Buffer capacity exceeded: buflen= %d buflim= %d\n",
              buflen, bp);
      }
      strncpy((buff)+bufpos, word, (buflen-bufpos));

      // Descend another layer
      int success = try_crack(target,
                            dicts, dicts_len, dict_pos+1,
                            buff, buflen, bp);

      if(success==1){             // Check for success
        strcpy(buf,buff);
        found = 1;
      }
      // No match, next iteration overwrites previous word
    }
  }
  // Tried all dictionary words at this level with no luck
  return found;
}
struct thread_data {
    //STRUCT TO HOLD DATA FOR PASSING
    char *target;
    dict_t **dicts;
    int dicts_len;
    int dict_pos;
    char *buf;
    int buflen;
    int bufpos;
    long thread_id;
    int num_threads;
};
void *pcrack(void *arg){
    //GATHER DATA FED IN FROM STRUCT
    struct thread_data *my_data;
    my_data = (struct thread_data *) arg;
    char *target = my_data->target;
    dict_t **dicts = my_data->dicts;
    int dicts_len = my_data->dicts_len;
    int dict_pos = my_data->dict_pos;
    char *buf = my_data->buf;
    int buflen = my_data->buflen;
    int bufpos = my_data->bufpos;
    long thread_id = my_data->thread_id;
    int num_threads = my_data->num_threads;
    long junk = 0; //junkretval

    //CREATE A BUFFER FOR THE THREAD TO WORK WITH
    char *buff = malloc(buflen *sizeof(char));
    dict_t *cur_dict = dicts[dict_pos];

    //FIGURE OUT HOW MANY ITERATIONS PER THREAD
    long iterations = dict_get_word_count(cur_dict);
    long my_start = thread_id * iterations / num_threads;
    long my_end = (thread_id + 1) * iterations / num_threads;

    //printf("I am thread %ld, mystart: %ld  myend: %ld\n",thread_id,my_start,my_end);//DEBUG
    for(int n=my_start; n<my_end; ++n){//loop the right amount of times
      if(pfound)//stop doing work
        continue;
      // Append a new word to the end of buf
      char *word = dict_get_word(cur_dict, n);
      int bp = bufpos + strlen(word);
      if(bp > buflen){
        fprintf(stderr,"WARNING: Buffer capacity exceeded: buflen= %d buflim= %d\n",
              buflen, bp);
      }
      strncpy((buff)+bufpos, word, (buflen-bufpos));
      // Descend another layer
      int success = try_crack(target,
                            dicts, dicts_len, dict_pos+1,
                            buff, buflen, bp);

      if(success==1){// Check for success
        strcpy(buf,buff);//COPY BACK TO GLOBAL BUFFER
        pfound = 1;
      }
    }
    return (void *) junk;//WHO CARES..
}
int try_crackpthread(char *target,
              dict_t **dicts, int dicts_len, int dict_pos,
      char *buf, int buflen, int bufpos, int num_threads)
{
    //create numthreads amount of initial threads, like the omp parallel shared
    //then give a portion of the for loop to each thread
    //communicate a found
    //signal on the global pfound

    pthread_t threads[num_threads];
    struct thread_data thread_data_array[num_threads];
    pfound = 0;//RESET GVAR

    for(long p=0; p<num_threads; p++){
      //LOAD THE DATA
      thread_data_array[p].target = target;
      thread_data_array[p].dicts = dicts;
      thread_data_array[p].dicts_len = dicts_len;
      thread_data_array[p].dict_pos = dict_pos;
      thread_data_array[p].buf = buf;
      thread_data_array[p].buflen = buflen;
      thread_data_array[p].bufpos = bufpos;
      thread_data_array[p].thread_id = p;
      thread_data_array[p].num_threads = num_threads;

      //pthread_create...
      pthread_create(&threads[p],NULL,pcrack,(void *) &thread_data_array[p]);
    }

    for(long p=0; p<num_threads; p++){
      //pthread_join...
      long retval;
      pthread_join(threads[p],(void **) &retval);
    }
    return pfound;
}
