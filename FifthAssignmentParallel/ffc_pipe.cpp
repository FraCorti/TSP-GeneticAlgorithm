/*
 * Simple file compressor using miniz and the FastFlow pipeline.
 *
 * miniz source code: https://github.com/richgel999/miniz
 * https://code.google.com/archive/p/miniz/
 * 
 */
/* Author: Massimo Torquati <massimo.torquati@unipi.it>
 * This code is a mix of POSIX C code and some C++ library call 
 * (mainly for strings manipulation).
 */


#include <miniz.h>
#include <string>
#include <ff/pipeline.hpp>
#include <ff/ff.hpp>
#include "utility.hpp"
using namespace ff;

struct Task {
  Task(unsigned char *ptr, size_t size, const std::string &name) :
      ptr(ptr), size(size), cmp_size(0), filename(name) {}

  unsigned char *ptr;
  size_t size;      // input size
  size_t cmp_size;  // size compressed file
  const std::string filename;
};

// 1st stage
struct Read : ff_node_t<Task> {
  Read(const char **argv, int argc) : argv(argv), argc(argc){}

  // ------------------- utility functions
  // It memory maps the input file and then assigns a task to
  // one Worker
  bool doWork(const std::string &fname, size_t size) {
    unsigned char *ptr = nullptr;
    if (!mapFile(fname.c_str(), size, ptr)) return false;
    Task *t = new Task(ptr, size, fname);
    ff_send_out(t); // sending to the next stage
    return true;
  }
  // walks through the directory tree rooted in dname
  bool walkDir(const std::string &dname, size_t size) {
    DIR *dir;
    if ((dir = opendir(dname.c_str())) == NULL) {
      perror("opendir");
      fprintf(stderr, "Error: opendir %s\n", dname.c_str());
      return false;
    }
    struct dirent *file;
    bool error = false;
    while ((errno = 0, file = readdir(dir)) != NULL) {
      struct stat statbuf;
      std::string filename = dname + "/" + file->d_name;
      if (stat(filename.c_str(), &statbuf) == -1) {
        perror("stat");
        fprintf(stderr, "Error: stat %s\n", filename.c_str());
        return false;
      }
      if (S_ISDIR(statbuf.st_mode)) {
        if (!isdot(filename.c_str())) {
          if (!walkDir(filename, statbuf.st_size)) error = true;
        }
      } else {
        if (!doWork(filename, statbuf.st_size)) error = true;
      }
    }
    if (errno != 0) {
      perror("readdir");
      error = true;
    }
    closedir(dir);
    return !error;
  }
  // -------------------

  Task *svc(Task *) {
    for (long i = 0; i < argc; ++i) {
      struct stat statbuf;
      if (stat(argv[i], &statbuf) == -1) {
        perror("stat");
        fprintf(stderr, "Error: stat %s\n", argv[i]);
        continue;
      }
      bool dir = false;
      if (S_ISDIR(statbuf.st_mode)) {
        success &= walkDir(argv[i], statbuf.st_size);
      } else {
        success &= doWork(argv[i], statbuf.st_size);
      }
    }
    return EOS;
  }

  void svc_end() {
    if (!success) {
      printf("Read stage: Exiting with (some) Error(s)\n");
      return;
    }
  }

  const char **argv;
  const int argc;
  bool success = true;
};
// 2nd stage
struct Compressor : ff_node_t<Task> {
  Task *svc(Task *task) {
    unsigned char *inPtr = task->ptr;
    size_t inSize = task->size;
    // get an estimation of the maximum compression size
    unsigned long cmp_len = compressBound(inSize);
    // allocate memory to store compressed data in memory
    unsigned char *ptrOut = new unsigned char[cmp_len];
    if (compress(ptrOut, &cmp_len, (const unsigned char *) inPtr, inSize) != Z_OK) {
      printf("Failed to compress file in memory\n");
      success = false;
      delete[] ptrOut;
      return GO_ON;
    }
    task->ptr = ptrOut;
    task->cmp_size = cmp_len;
    ff_send_out(task);

    unmapFile(inPtr, inSize);
    return GO_ON;
  }
  void svc_end() {
    if (!success) {
      printf("Comprssor stage: Exiting with (some) Error(s)\n");
      return;
    }
  }
  bool success = true;
};
// 3rd stage
struct Write : ff_node_t<Task> {
  Task *svc(Task *task) {
    const std::string outfile = task->filename + ".zip";
    // write the compressed data into disk
    success &= writeFile(outfile, task->ptr, task->cmp_size);
    if (success && REMOVE_ORIGIN) {
      unlink(task->filename.c_str());
    }
    delete[] task->ptr;
    delete task;
    return GO_ON;
  }
  void svc_end() {
    if (!success) {
      printf("Write stage: Exiting with (some) Error(s)\n");
      return;
    }
  }
  bool success = true;
};

static inline void usage(const char *argv0) {
  printf("--------------------\n");
  printf("Usage: %s file-or-directory [file-or-directory]\n", argv0);
  printf("\nModes: COMPRESS ONLY\n");
  printf("--------------------\n");
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    usage(argv[0]);
    return -1;
  }
  argc--;
  Read reader(const_cast<const char**>(&argv[1]), argc);
  Compressor compress;
  Write writer;
  ff_Pipe pipe(reader, compress, writer);
  if (pipe.run_and_wait_end()<0) {
    error("running pipeline\n");
    return -1;
  }
  bool success = true;
  success &= reader.success;
  success &= compress.success;
  success &= writer.success;
  if (success) printf("Done.\n");

  return 0;
}