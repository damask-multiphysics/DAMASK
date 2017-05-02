/* Unix */
#include <stdio.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>

/* http://stackoverflow.com/questions/30279228/is-there-an-alternative-to-getcwd-in-fortran-2003-2008 */

int isdirectory_c(const char *dir){
  struct stat statbuf;
  if(stat(dir, &statbuf) != 0)
    return 0;
  return S_ISDIR(statbuf.st_mode);
}


void getcurrentworkdir_c(char cwd[], int *stat ){
  char cwd_tmp[1024];
  if(getcwd(cwd_tmp, sizeof(cwd_tmp)) == cwd_tmp){
    strcpy(cwd,cwd_tmp);
    *stat = 0;
  }
  else{
    *stat = 1;
  }
}


void gethostname_c(char hostname[], int *stat ){
  char hostname_tmp[1024];
  if(gethostname(hostname_tmp, sizeof(hostname_tmp)) == 0){
    strcpy(hostname,hostname_tmp);
    *stat = 0;
  }
  else{
    *stat = 1;
  }
}
