// Main entry point to crack passwords. Parallelized using PThreads

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <crack.h>
#include <ctype.h>
#include <pthread.h>

int main(int argc, char **argv) {
  if(argc < 3) {
    printf("usage: %s <encrypted_file> <dict1> [dict2] ...\n",argv[0]);
    printf("  <encrypted_file> : encrypted password file, one per line\n");
    printf("  <dict1>          : dictionary to try for passwords for word 1\n");
    printf("  [dict2]          : additional dictionary to try for word 2\n");
    printf("                   : further dictionaries may be specified for more words\n");
    return 0;
  }

  //check env variable for number of threads
  int nthreads = 4;
  char *nthreads_str = getenv("PASSCRACK_NUMTHREADS");
  if(nthreads_str != NULL){
    nthreads = atoi(nthreads_str);
  }

  // Load the passwords from a file
  dict_t *passwords = dict_load(argv[1]);
  printf("found %d passwords to crack\n",dict_get_word_count(passwords));

  // Load all dictionaries of words
  char **dict_files = &(argv[2]);
  int dicts_len = argc-2;
  dict_t **dicts = dict_load_dicts(dict_files, dicts_len);

  // Show information on loaded dictionaries and compute the maximum
  // length of the temporary buffer required to accomodate all
  // possible passwords.
  int buflen = 0;
  for(int i=0; i<dicts_len; i++){
    printf("dict %d: %d words\n",i,dict_get_word_count(dicts[i]));
    buflen += dict_get_longest_word_length(dicts[i]);
  }
  char *buf = malloc(buflen * sizeof(char));

  // Loop through each password in the password file and attempt to
  // crack it.
  int successes = 0;
  for(int i=0; i<dict_get_word_count(passwords); i++){
    // int i = 0;
      char *encrypted = dict_get_word(passwords,i);
      int success = try_crackpthread(encrypted, dicts, dicts_len, 0,
				       buf, buflen, 0,nthreads);
      if(success == 1){
	printf("%3d: SUCCES: %s <-- %s\n",i,encrypted,buf);
	successes++;
      }
      else{
	printf("%3d: FAILED: %s <-- %s\n",i,encrypted,"???");
      }
    }
  printf("%d / %d passwords cracked\n",successes,dict_get_word_count(passwords));

  // Free up memory and bail out
  dict_free(passwords);
  dict_free_dicts(dicts, dicts_len);
  free(buf);
  return 0;
}

