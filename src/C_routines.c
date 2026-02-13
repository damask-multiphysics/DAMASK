// SPDX-License-Identifier: AGPL-3.0-or-later
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>
#include <string.h>
#include <pwd.h>
#include <limits.h>
#include <ISO_Fortran_binding.h>
#include <stdbool.h>

#ifdef GRID
#include <zlib.h>
#endif

#ifdef FYAML
#include <libfyaml.h>
#endif

#ifndef MARC_SOURCE
extern bool f_sigint, f_sigusr1, f_sigusr2;

static void signalHandler(int signum) {
  if      (signum == SIGINT)  f_sigint  = true;
  else if (signum == SIGUSR1) f_sigusr1 = true;
  else if (signum == SIGUSR2) f_sigusr2 = true;
}

void init_signal_c() {
  signal(SIGINT, signalHandler);
  signal(SIGUSR1, signalHandler);
  signal(SIGUSR2, signalHandler);
}
#endif

int set_cwd_c(const char *cwd){
  return chdir(cwd);
}

#ifndef OLD_STYLE_C_TO_FORTRAN_STRING
void get_cwd_c(CFI_cdesc_t *cwd, int *stat){
  char cwd_tmp[PATH_MAX]; // getcwd returns a string <= PATH_MAX (incl NULL)
  if (getcwd(cwd_tmp, sizeof(cwd_tmp))){
    if (CFI_allocate(cwd, (CFI_index_t *)0, (CFI_index_t *)0, strlen(cwd_tmp)) == 0){
      memcpy(cwd->base_addr, cwd_tmp, strlen(cwd_tmp));
      *stat = 0;
    }
    else *stat = 2;
  }
  else *stat = 1;
}

void get_hostname_c(CFI_cdesc_t *hostname, int *stat){
  char hostname_tmp[_POSIX_HOST_NAME_MAX+1]; // for a host name of size HOST_NAME_MAX, NULL is not included
  *stat = gethostname(hostname_tmp, sizeof(hostname_tmp));

  if (*stat == 0){
    if (CFI_allocate(hostname, (CFI_index_t *)0, (CFI_index_t *)0, strlen(hostname_tmp)) == 0)
      memcpy(hostname->base_addr, hostname_tmp, strlen(hostname_tmp));
    else *stat = 1;
  }
}


void get_username_c(CFI_cdesc_t *username, int *stat){
  struct passwd *pw = getpwuid(getuid());
  if (pw){
    if (CFI_allocate(username, (CFI_index_t *)0, (CFI_index_t *)0, strlen(pw->pw_name)) == 0){
      memcpy(username->base_addr, pw->pw_name, strlen(pw->pw_name));
      *stat = 0;
    }
    else *stat = 2;
  }
  else *stat = 1;
}
#else
void get_cwd_c(char cwd[], int *stat ){
  char cwd_tmp[4096+1];
  if (getcwd(cwd_tmp, sizeof(cwd_tmp))){
    strcpy(cwd, cwd_tmp); // getcwd guarantees a NULL-terminated string
    *stat = 0;
  }
  else *stat = 1;
}


void get_hostname_c(char hostname[], int *stat){
  char hostname_tmp[256];
  if (gethostname(hostname_tmp, sizeof(hostname_tmp)) == 0){
    strncpy(hostname, hostname_tmp, strlen(hostname_tmp)+1); // gethostname does not guarantee a NULL-terminated string
    *stat = 0;
  }
  else *stat = 1;
}


void get_username_c(char username[], int *stat){
  struct passwd *pw = getpwuid(getuid());
  if (pw && strlen(pw->pw_name) <= 256){
    strncpy(username, pw->pw_name,strlen(pw->pw_name)+1);
    *stat = 0;
  }
  else *stat = 1;
}
#endif

bool isatty_stdout_c(){
  return isatty(STDOUT_FILENO) == 1;
}

bool isatty_stderr_c(){
  return isatty(STDERR_FILENO) == 1;
}

bool isatty_stdin_c(){
  return isatty(STDIN_FILENO) == 1;
}


#ifdef GRID
void inflate_c(const uLong *s_deflated, const uLong *s_inflated,
               const Byte deflated[*s_deflated], Byte inflated[*s_inflated], int *stat){
  /* make writable copy, uncompress will write to it */
  uLong s_inflated_;
  s_inflated_ = *s_inflated;

  // https://stackoverflow.com/questions/51334741
  if (uncompress((Bytef *)inflated, &s_inflated_, (Bytef *)deflated, *s_deflated) == Z_OK){
    if (*s_inflated == s_inflated_) *stat = 0;
  }
  else{
    *stat = 1;
    memset(inflated, 0, (size_t)*s_inflated);
  }
}
#endif


#ifdef FYAML
#ifndef OLD_STYLE_C_TO_FORTRAN_STRING
void to_flow_c(CFI_cdesc_t *flow, const char *mixed, int *stat){
  struct fy_document *fyd = NULL;
  enum fy_emitter_cfg_flags emit_flags = FYECF_MODE_FLOW_ONELINE
                                       | FYECF_WIDTH_INF
                                       | FYECF_STRIP_LABELS
                                       | FYECF_STRIP_TAGS
                                       | FYECF_STRIP_DOC
                                       | FYECF_DOC_START_MARK_OFF;

  fyd = fy_document_build_from_string(NULL, mixed, -1);
  if (!fyd) {
    *stat = 1;
    return;
  }

  int err = fy_document_resolve(fyd);
  if (err) {
    *stat = 2;
    return;
  }

  char *flow_tmp = fy_emit_document_to_string(fyd, emit_flags); /* ends with new line */
  if (CFI_allocate(flow, (CFI_index_t *)0, (CFI_index_t *)0, strlen(flow_tmp)-1) == 0){
    memcpy(flow->base_addr, flow_tmp, strlen(flow_tmp)-1);
  }
  else {
    *stat = 3;
    return;
  }

  free(flow_tmp);
  fy_document_destroy(fyd);
  *stat = 0;
}
#else
void to_flow_c(char **flow, long* length_flow, const char *mixed){
  struct fy_document *fyd = NULL;
  enum fy_emitter_cfg_flags emit_flags = FYECF_MODE_FLOW_ONELINE
                                       | FYECF_WIDTH_INF
                                       | FYECF_STRIP_LABELS
                                       | FYECF_STRIP_TAGS
                                       | FYECF_STRIP_DOC
                                       | FYECF_DOC_START_MARK_OFF;

  fyd = fy_document_build_from_string(NULL, mixed, -1);
  if (!fyd) {
    *length_flow = -1;
    return;
  }
  int err = fy_document_resolve(fyd);
  if (err) {
    *length_flow = -1;
    return;
  }

  *flow = fy_emit_document_to_string(fyd, emit_flags);
  *length_flow = (long) strlen(*flow);

  fy_document_destroy(fyd);
}

#endif
#endif
