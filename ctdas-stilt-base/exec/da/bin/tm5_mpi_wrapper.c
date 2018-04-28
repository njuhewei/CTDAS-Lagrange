#include <stdio.h>
#include <errno.h>
#include "mpi.h"

#define CMDLENGTH 200

int safe_system (const char *command);

int main(int argc, char *argv[]) { 

  int ierr;
  int myrank;
  char cmd[CMDLENGTH];

  ierr = MPI_Init(&argc, &argv);
  if (argc != 2) {
    fprintf(stderr, "[tm5_mpi_wrapper] Expecting 1 argument, got %d.\n",argc);
    exit(-1);
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  snprintf(cmd,CMDLENGTH,"%s %03d",argv[1],myrank);
  //snprintf(cmd,CMDLENGTH,"%s %03d >> tm5.%03d.log",argv[1],myrank,myrank);//
  printf( "MPI rank %d about to execute command \"%s\".\n",myrank,cmd );
  ierr = safe_system(cmd);
  if(ierr != 0) {
     MPI_Abort(MPI_COMM_WORLD,ierr);
     exit(ierr);
  }
  ierr = MPI_Finalize(  );
  exit(ierr);
}





int safe_system (const char *command) {
  int pid, status;

  if (command == 0)
    return 1;
  pid = fork();
  if (pid == -1) // fork failed
    return -1;
  if (pid == 0) {  // then this is the child
    char *argv[4];
    argv[0] = "sh";
    argv[1] = "-c";
    argv[2] = (char*)command;
    argv[3] = 0;
    execv("/bin/sh", argv);
    _exit(127);
  }
  do {
    if (waitpid(pid, &status, 0) == -1) {
      if (errno != EINTR)
        return -1;
    } else
      return status;
  } while(1);
}
