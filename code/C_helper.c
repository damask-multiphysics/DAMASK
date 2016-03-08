/* Unix */
#include <unistd.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>

#define GETCWD getcwd


void getCurrentWorkDir(char *str,int *stat){
  if(GETCWD(str, sizeof(str)) == str){
    *stat = 0;
  } 
  else{
    *stat = 1;
  }
}

int isDirectory(const char *path){
  struct stat statbuf;
  if(stat(path, &statbuf) != 0)
    return 0;
  return S_ISDIR(statbuf.st_mode);
}
