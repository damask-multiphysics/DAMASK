/* Unix */
#include <stdio.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>
#include <signal.h>
#include <pwd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "zlib.h"

#define PATHLEN 4096
#define STRLEN   256


int setcwd_c(const char *cwd){
  return chdir(cwd);
}


void getcwd_c(char cwd[], int *stat ){
  char cwd_tmp[PATHLEN+1];
  if(getcwd(cwd_tmp, sizeof(cwd_tmp))){
    strcpy(cwd,cwd_tmp); // getcwd guarantees a NULL-terminated string
    *stat = 0;
  }
  else{
    *stat = 1;
  }
}


void gethostname_c(char hostname[], int *stat){
  char hostname_tmp[STRLEN];
  if(gethostname(hostname_tmp, sizeof(hostname_tmp)) == 0){
    strncpy(hostname,hostname_tmp,sizeof(hostname_tmp)+1); // gethostname does not guarantee a NULL-terminated string
    *stat = 0;
  }
  else{
    *stat = 1;
  }
}


void getusername_c(char username[], int *stat){
  struct passwd *pw = getpwuid(getuid());
  if(pw && strlen(pw->pw_name) <= STRLEN){
    strncpy(username,pw->pw_name,STRLEN+1);
    *stat = 0;
  }
  else{
    *stat = 1;
  }
}


void signalterm_c(void (*handler)(int)){
  signal(SIGTERM, handler);
}

void signalusr1_c(void (*handler)(int)){
  signal(SIGUSR1, handler);
}

void signalusr2_c(void (*handler)(int)){
  signal(SIGUSR2, handler);
}


void inflate_c(const uLong *s_deflated, const uLong *s_inflated, const Byte deflated[], Byte inflated[]){
  /* make writable copy, uncompress will write to it */
  uLong s_inflated_,i;
  s_inflated_ = *s_inflated;

  if(uncompress((Bytef *)inflated, &s_inflated_, (Bytef *)deflated, *s_deflated) == Z_OK)
    return;
  else{
    for(i=0;i<*s_inflated;i++){
      inflated[i] = 0;
    }
  }
}
