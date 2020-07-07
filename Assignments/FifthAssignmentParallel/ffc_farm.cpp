//
// Created by Francesco on 24/05/20.
//

/*
 * Simple file compressor using miniz and the FastFlow pipeline.
 *
 * miniz source code: https://github.com/richgel999/miniz
 * https://code.google.com/archive/p/miniz/
 *
 */

#include <miniz.h>
#include <string>
#include <iostream>
#include <ff/ff.hpp>
#include <ff/farm.hpp>
#include <utility.hpp>
#include <map>
#include <filesystem>
using namespace ff;

struct Task {
  Task(unsigned char *chunkPointer,
       unsigned char *beginFilePointer,
       const std::string &chunkName,
       const std::string &completeFileName,
       size_t size) :
      ptr(chunkPointer),
      beginFilePointer(beginFilePointer),
      size(size),
      filename(chunkName),
      completeFilename(completeFileName),
      cmp_size(0) {}
  unsigned char *ptr;
  unsigned char *beginFilePointer;
  const std::string filename;
  const std::string completeFilename;
  size_t size;
  size_t cmp_size;  // compression size, set to zero and initialized in Compress() method
};

//! 1st stage
struct Read : ff_node_t<Task> {
  Read(const char **argv, int argc, size_t threshold) : argv(argv), argc(argc), threshold(threshold) {}

  //! Create chunks and pass them to second stage (Compress)
  bool doWork(const std::string &fname, size_t size) {
    unsigned char *filePointer = nullptr;

    //! number of chunk needed to manage the file
    size_t chunksNumber = size / threshold;
    size_t remained = size % threshold;
    size_t extraChunk = 0;
    if (remained != 0) {
      extraChunk = 1;  // manage extra chunk to update chunksNumber
    }
    int currentChunk = 0; //! identify chunks
    for (; currentChunk < chunksNumber; currentChunk++) {
      if (!mapFile(fname.c_str(), size, filePointer)) return false;

      Task *t = new Task(filePointer + currentChunk * threshold,
                         filePointer,
                         "files/" + fname + "/" + fname + "_Chunk" + std::to_string(currentChunk),
                         fname,
                         threshold);
      ff_send_out(t);
    }

    //! manage extra (last) chunk
    if (extraChunk != 0) {
      currentChunk++;
      if (!mapFile(fname.c_str(), size, filePointer)) return false;
      Task *t = new Task(filePointer + currentChunk * threshold,
                         filePointer,
                         "files/" + fname + "/" + fname + "_Chunk" + std::to_string(currentChunk),
                         fname,
                         threshold);
      ff_send_out(t);
    }
    return true;
  }

  //! walks through the directory tree dname
  bool walkDir(const std::string &dname, size_t size) {
    DIR *dir;
    if ((dir = opendir(dname.c_str())) == NULL) {
      perror("opendir");
      fprintf(stderr, "Error: opendir %s\n", dname.c_str());
      return false;
    }
    struct dirent *file;
    bool error = false;

    //! read current directory (dir)
    while ((errno = 0, file = readdir(dir)) != NULL) {

      //! structure containing file information
      struct stat statbuf;
      std::string filename = dname + "/" + file->d_name;

      //! retrieve file information and storing it in statbuf
      if (stat(filename.c_str(), &statbuf) == -1) {
        perror("stat");
        fprintf(stderr, "Error: stat %s\n", filename.c_str());
        return false;
      }

      //! if is a dir recurse on walkDir
      if (S_ISDIR(statbuf.st_mode)) {
        if (!isdot(filename.c_str())) {
          if (!walkDir(filename, statbuf.st_size)) error = true;
        }
      }
        //! if is a file pass the job to doWork method with the size and path of the file
      else {
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
  size_t threshold;
};

//! 2nd stage
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
    unmapFile(task->beginFilePointer, inSize);
    return GO_ON;
  }
  void svc_end() {
    if (!success) {
      printf("Compressor stage: Exiting with (some) Error(s)\n");
      return;
    }
  }
  bool success = true;
};

//! 3rd stage
struct Write : ff_node_t<Task> {
  Task *svc(Task *task) {
    const std::string outfile = task->filename + ".zip";

    //! check if exists a folder for the zipped chunk inside files, if not create it
    if (!std::filesystem::exists("files/" + task->completeFilename)) {
      std::filesystem::create_directories("files/" + task->completeFilename);
    }

    //! write the compressed data into disk
    success &= writeFile(outfile, task->ptr, task->cmp_size);
    if (success && REMOVE_ORIGIN) {
      unlink(task->filename.c_str());
    }
    delete[] task->ptr;
    delete task;
    return GO_ON;
  }
  void svc_end() {

    //! remove old folder and create new one
    std::filesystem::remove_all("mergeZippedFiles");
    std::filesystem::create_directories("mergeZippedFiles");

    std::string enterFilesDirectory = "files/";
    int check = chdir(enterFilesDirectory.c_str());
    if (check !=0){
      success = false;
    }

    //! iterate in the files/ directory and compress all the folders with the chunk
    for (const auto &entry : std::filesystem::directory_iterator("../files")) {
      std::string currentDirToZip = std::string(entry.path().filename());
      std::string command_string = " tar -zcf ../mergeZippedFiles/"+ currentDirToZip+ ".tar.gz " + currentDirToZip;

      int r = system(command_string.c_str());
      if(r!= 0){
        success = false;
      }
    }

    //! delete files directory
    std::filesystem::remove_all("../files");

    if (!success) {
      printf("Write stage: Exiting with (some) Error(s)\n");
      return;
    }
  }
  //std::map<std::string, int> files;  // <std::string originalFilename, int chunksNumber> files
  bool success = true;
};

static inline void usage(const char *argv0) {
  printf("--------------------\n");
  printf("Usage: %s nworkers threshold [file-or-directory]\n", argv0);
  printf("\nModes: COMPRESS ONLY\n");
  printf("--------------------\n");
}

int main(int argc, char *argv[]) {
  if (argc < 4) {
    usage(argv[0]);
    return -1;
  }
  argc--;

  const size_t nworkers = std::stol(argv[1]);
  argc--;

  const size_t threshold = std::stol(argv[2]);
  argc--;

  Read reader(const_cast<const char **>(&argv[3]), argc, threshold);
  Write writer;

  //! Remove old "files" directory and create new one
  std::filesystem::remove_all("files");
  std::filesystem::create_directories("files");

  //! Farm creation
  ff_farm farm;
  std::vector<ff_node *> W;
  for (size_t i = 0; i < nworkers; ++i)
    W.push_back(new Compressor);
  farm.add_workers(W);
  farm.add_collector(nullptr);
  farm.cleanup_workers();

  ff_pipeline pipe;
  pipe.add_stage(reader);
  pipe.add_stage(farm);
  pipe.add_stage(writer);

  ff::ffTime(ff::START_TIME);
  //! start the pipeline
  if (pipe.run_and_wait_end() < 0) {
    error("running pipeline\n");
    return -1;
  }

  ff::ffTime(ff::STOP_TIME);
  std::cout << nworkers << " " << ff::ffTime(ff::GET_TIME) << std::endl;
  bool success = true;
  success &= reader.success;
  success &= writer.success;
  if (success) printf("Done.\n");

  return 0;
}