//
// miniz source code: https://github.com/richgel999/miniz
// https://code.google.com/archive/p/miniz/
// 
// This is a modified version of the example3.c file distributed with the miniz.c, it inflates/deflates 
// block by block. 
// --------------------
// example3.c - Demonstrates how to use miniz.c's deflate() and inflate() functions for simple file compression.
// Public domain, May 15 2011, Rich Geldreich, richgel99@gmail.com. See "unlicense" statement at the end of tinfl.c.
// For simplicity, this example is limited to files smaller than 4GB, but this is not a limitation of miniz.c.
// -------------------
//
// Author: Massimo Torquati <massimo.torquati@unipi.it>
// This code is a mix of POSIX C code and some C++ library call (mainly for the strings).
//

#if !defined _UTILITY_HPP
#define _UTILITY_HPP

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <dirent.h>
#include <sys/stat.h>

#include <algorithm>
#include <string>

#include <miniz.h>

#define FASTER_COMPRESSION
#define BUF_SIZE (1024 * 1024)

static bool REMOVE_ORIGIN=false;
static unsigned char s_inbuf[BUF_SIZE];
static unsigned char s_outbuf[BUF_SIZE];



// map the file pointed by filepath in memory
static inline bool mapFile(const char fname[], size_t size, unsigned char *&ptr) {
    // open input file.
    int fd = open(fname,O_RDONLY);
    if (fd<0) {
	printf("Failed opening file %s\n", fname);
	return false;
    }
#if 0    
    struct stat s;
    if (fstat (fd, &s)) {
	printf("Failed to stat file %s\n", fname);
	return false;
    }
    // checking for regular file type
    if (!S_ISREG(s.st_mode)) return false;
    // check the size of the file
    assert(size == s.st_size);
#endif     
    // map all the file in memory
    ptr = (unsigned char *) mmap (0, size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (ptr == MAP_FAILED) {
	printf("Failed to memory map file %s\n", fname);
	return false;
    }
    close(fd);
    return true;
}
// unmap the file from memory
static inline void unmapFile(unsigned char *ptr, size_t size) {
    if (munmap(ptr, size)<0) {
	printf("Failed to unmap file\n");
    }
}
// write size bytes starting from ptr into a file pointed by filename
static inline bool writeFile(const std::string &filename, unsigned char *ptr, size_t size) {
    FILE *pOutfile = fopen(filename.c_str(), "wb");
    if (!pOutfile) {
	printf("Failed opening output file %s!\n", filename.c_str());
	return false;
    }
    if (fwrite(ptr, 1, size, pOutfile) != size) {
	printf("Failed writing to output file %s\n", filename.c_str());
	return false;
    }
    if (fclose(pOutfile) != 0) return false;
    return true;
}


static inline bool isdot(const char dir[]) {
  int l = strlen(dir);  
  if ( (l>0 && dir[l-1] == '.') ) return true;
  return false;
}

// --------------------------------------------------------------------------
// compress the input file (fname) having size infile_size
// it returns 0 for success and -1 if something went wrong
// if removeOrigin is true, the source file will be removed if successfully compressed
static inline int compressFile(const char fname[], size_t infile_size,
			       const bool removeOrigin=REMOVE_ORIGIN) {
    // define the output file name 
    const std::string infilename(fname);
    std::string outfilename = std::string(fname) + ".zip";

#if defined(FASTER_COMPRESSION)

    unsigned char *ptr = nullptr;
    if (!mapFile(fname, infile_size, ptr)) return -1;
    // get an estimation of the maximum compression size
    unsigned long cmp_len = compressBound(infile_size);
    // allocate memory to store compressed data in memory
    unsigned char *ptrOut = new unsigned char[cmp_len];
    if (compress(ptrOut, &cmp_len, (const unsigned char *)ptr, infile_size) != Z_OK) {
	printf("Failed to compress file in memory\n");	   
	delete [] ptrOut;
	return -1;
    }
    // write the compressed data into disk 
    bool success = writeFile(outfilename, ptrOut, cmp_len);
    if (success && removeOrigin) {
	unlink(fname);
    }
    unmapFile(ptr, infile_size);
    delete [] ptrOut;

#else    
    FILE *pInfile;
    FILE *pOutfile;

    // Open input file.
    pInfile = fopen(fname, "rb");
    if (!pInfile) {
	printf("Failed opening input file!\n");
	return -1;
    }
    char *fnameOut = const_cast<char *>(outfilename.c_str());
    // Open output file.
    pOutfile = fopen(fnameOut, "wb");
    if (!pOutfile) {
	printf("Failed opening output file!\n");
	fclose(pInfile);
	return -1;
    }

    z_stream stream;
    // Init the z_stream
    memset(&stream, 0, sizeof(stream));
    stream.next_in = s_inbuf;
    stream.avail_in = 0;
    stream.next_out = s_outbuf;
    stream.avail_out = BUF_SIZE;
        
    size_t infile_remaining = infile_size;
    if (deflateInit(&stream, Z_BEST_COMPRESSION) != Z_OK) {
	printf("deflateInit() failed!\n");
	fclose(pInfile); fclose(pOutfile);
	return -1;
    }    
    for ( ; ; ) {
	if (!stream.avail_in)  {
	    // Input buffer is empty, so read more bytes from input file.
	    size_t n = std::min((size_t)BUF_SIZE, infile_remaining);
	    if (fread(s_inbuf, 1, n, pInfile) != n) {
		printf("Failed reading from input file!\n");
		fclose(pInfile); fclose(pOutfile);
		return -1;
	    }
	    stream.next_in    = s_inbuf;
	    stream.avail_in   = n;	  
	    infile_remaining -= n;
	}      
	int status = deflate(&stream, infile_remaining ? Z_NO_FLUSH : Z_FINISH);
	if ((status == Z_STREAM_END) || (!stream.avail_out))  {
	    // Output buffer is full, or compression is done, so write buffer to output file.
	    size_t n = BUF_SIZE - stream.avail_out;
	    if (fwrite(s_outbuf, 1, n, pOutfile) != n) {
		printf("Failed writing to output file!\n");
		fclose(pInfile); fclose(pOutfile);
		return -1;
	    }
	    stream.next_out  = s_outbuf;
	    stream.avail_out = BUF_SIZE;
	}
	
	if (status == Z_STREAM_END)  break; // done
	else if (status != Z_OK) {
	    printf("deflate() failed with status %i!\n", status);
	    fclose(pInfile); fclose(pOutfile);
	    return -1;
	}
    } // for
    if (deflateEnd(&stream) != Z_OK) {
	printf("deflateEnd() failed!\n");
	fclose(pInfile); fclose(pOutfile);
	return -1;
    }
    //printf("Total input bytes: %u\n", (mz_uint32)stream.total_in);
    //printf("Total output bytes: %u\n", (mz_uint32)stream.total_out);
    fclose(pInfile); fclose(pOutfile);
    if (removeOrigin) unlink(fname);
#endif 
    return 0;
}


// uncompress the input file (fname) having size infile_size
// it returns 0 for success and -1 if something went wrong
// if removeOrigin is true, the input compressed file will be removed if successfully uncompressed
static inline int decompressFile(const char fname[], size_t infile_size,
				 const bool removeOrigin=REMOVE_ORIGIN) {
    FILE *pInfile;
    FILE *pOutfile;
    // Open input file.
    pInfile = fopen(fname, "rb");
    if (!pInfile) {
	printf("Failed opening input file!\n");
	return -1;
    }
    // define the output file name
    // if the input file does not have a ".zip" extention
    // the output file name will terminate with "_decomp"
    // and the original file will not be removed even if
    // removeOrigin is true
    const std::string infilename(fname);
    std::string outfilename;
    int n = infilename.find(".zip");
    if (n>0) outfilename = infilename.substr(0,n);
    else     outfilename = infilename + "_decomp";
    char *fnameOut = const_cast<char *>(outfilename.c_str());

    // Open output file.
    pOutfile = fopen(fnameOut, "wb");
    if (!pOutfile) {
	printf("Failed opening output file!\n");
	fclose(pInfile);
	return -1;
    }

    z_stream stream;
    // Init the z_stream
    memset(&stream, 0, sizeof(stream));
    stream.next_in = s_inbuf;
    stream.avail_in = 0;
    stream.next_out = s_outbuf;
    stream.avail_out = BUF_SIZE;
    
    size_t infile_remaining = infile_size;
    if (inflateInit(&stream)) {
	printf("inflateInit() failed!\n");
	fclose(pInfile); fclose(pOutfile);
	return -1;
    }
    
    for ( ; ; ) {
      if (!stream.avail_in)   {
	  // Input buffer is empty, so read more bytes from input file.
	  size_t n = std::min((size_t)BUF_SIZE, infile_remaining);
	  if (fread(s_inbuf, 1, n, pInfile) != n) {
	      printf("Failed reading from input file!\n");
	      fclose(pInfile); fclose(pOutfile);
	      return -1;
	  }
	  stream.next_in    = s_inbuf;
	  stream.avail_in   = n;	  
	  infile_remaining -= n;
      }
      int status = inflate(&stream, Z_SYNC_FLUSH);
      if ((status == Z_STREAM_END) || (!stream.avail_out)) {
	  // Output buffer is full, or decompression is done, so write buffer to output file.
	  size_t n = BUF_SIZE - stream.avail_out;
	  if (fwrite(s_outbuf, 1, n, pOutfile) != n) {
	      printf("Failed writing to output file!\n");
	      fclose(pInfile); fclose(pOutfile);
	      return -1;
	  }
	  stream.next_out = s_outbuf;
	  stream.avail_out = BUF_SIZE;
      }      
      if (status == Z_STREAM_END)  break; // done
      else if (status != Z_OK) {
	  printf("inflate() failed with status %i!\n", status);
	  fclose(pInfile); fclose(pOutfile);
	  return -1;
      }
    }// for
    if (inflateEnd(&stream) != Z_OK) {
	printf("inflateEnd() failed!\n");
	fclose(pInfile); fclose(pOutfile);
	return -1;
    }
    //printf("Total input bytes: %u\n", (mz_uint32)stream.total_in);
    //printf("Total output bytes: %u\n", (mz_uint32)stream.total_out);
    fclose(pInfile); fclose(pOutfile);
    if (n>0 && removeOrigin) unlink(fname);
    return 0;
}
// --------------------------------------------------------------------------



// returns false in case of error
static inline bool doWork(const char fname[], size_t size, const bool comp) {
    if (comp) {
	if (compressFile(fname, size, REMOVE_ORIGIN)<0) return false;
    } else {
	if (decompressFile(fname, size, REMOVE_ORIGIN)<0) return false;
    }
    return true;
}

// returns false in case of error
static inline bool walkDir(const char dname[], size_t size, const bool comp) {

    if (chdir(dname) == -1) {
	perror("chdir");
	fprintf(stderr, "Error: chdir %s\n", dname);
	return false;
    }
    DIR *dir;	
    if ((dir=opendir(".")) == NULL) {
	perror("opendir");
	fprintf(stderr, "Error: opendir %s\n", dname);
	return false;
    }
    struct dirent *file;
    bool error=false;
    while((errno=0, file =readdir(dir)) != NULL) {
	struct stat statbuf;
	if (stat(file->d_name, &statbuf)==-1) {
	    perror("stat");
	    fprintf(stderr, "Error: stat %s\n", file->d_name);
	    return false;
	}
	if(S_ISDIR(statbuf.st_mode)) {
	    if ( !isdot(file->d_name) ) {
		if (walkDir(file->d_name,statbuf.st_size, comp)) chdir("..");
		else error  = true;
	    }
	} else {
	    if (!doWork(file->d_name, statbuf.st_size, comp)) error = true;
	}
    }
    if (errno != 0) { perror("readdir"); error=true; }
    closedir(dir);
    return !error;
}



#endif 
