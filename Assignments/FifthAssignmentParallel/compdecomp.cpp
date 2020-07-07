/*
 * miniz source code: https://github.com/richgel999/miniz
 * https://code.google.com/archive/p/miniz/
 * 
 * This is an extended version of the example3.c file distributed with the miniz.c.
 * --------------------
 * example3.c - Demonstrates how to use miniz.c's deflate() and inflate() functions for simple file compression.
 * Public domain, May 15 2011, Rich Geldreich, richgel99@gmail.com. See "unlicense" statement at the end of tinfl.c.
 * For simplicity, this example is limited to files smaller than 4GB, but this is not a limitation of miniz.c.
 * -------------------
 *
 */
/* Author: Massimo Torquati <massimo.torquati@unipi.it>
 * This code is a mix of POSIX C code and some C++ library call 
 * (mainly for strings manipulation).
 */

#include "utility.hpp"

static inline void usage(const char *argv0) {
    printf("--------------------\n");
    printf("Usage: %s c|d|C|D file-or-directory [file-or-directory]\n",argv0);
    printf("\nModes:\n");
    printf("c - Compresses file infile to a zlib stream into outfile\n");
    printf("d - Decompress a zlib stream from infile into outfile\n");
    printf("C - Like c but remove the input file\n");
    printf("D - Like d but remove the input file\n");
    printf("--------------------\n");
}


int main(int argc, char *argv[]) {
  if (argc < 3) {
      usage(argv[0]);
      return -1;
  }
  const char *pMode = argv[1];
  --argc;
  if (!strchr("cCdD", pMode[0])) {
      printf("Invalid option!\n\n");
      usage(argv[0]);
      return -1;
  }
  const bool compress = ((pMode[0] == 'c') || (pMode[0] == 'C'));
  REMOVE_ORIGIN = ((pMode[0] == 'C') || (pMode[0] == 'D'));

  bool success = true;
  while(argc>1) {
      struct stat statbuf;
      if (stat(argv[argc], &statbuf)==-1) {
	  perror("stat");
	  fprintf(stderr, "Error: stat %s\n", argv[argc]);
	  --argc;
	  continue;
      }
      bool dir=false;
      if (S_ISDIR(statbuf.st_mode))  {
	  success &= walkDir(argv[argc], statbuf.st_size, compress);
      } else {
	  success &= doWork(argv[argc], statbuf.st_size, compress);
      }
      --argc;
  }
  if (!success) {
      printf("Exiting with (some) Error(s)\n");
      return -1;
  }
  printf("Exiting with Success\n");
  return 0;
}