//
// Created by checco on 23/05/20.
//

#ifndef FOURTHASSIGNMENTPARALLEL__IMAGEPROCESSING_H_
#define FOURTHASSIGNMENTPARALLEL__IMAGEPROCESSING_H_
#include <iostream>
#include <string>
#include <filesystem>
#include <utility>
#include <vector>
#include <functional>
#include <chrono>
#include <opencv2/core.hpp>
#include "opencv2/imgproc.hpp"
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <algorithm>
#include <random>
#include <thread>
#include <mkl.h>

class ImageProcessing {
 private:
  std::vector<std::string> filePaths;
  std::vector<cv::Mat> images;
  int numberFilesReplica;
  int numberWorkers;
  inline void loadImagesFromDisk(std::vector<cv::Mat> &images) {
    auto imagesIterator = images.begin();

    //! load images disk from their path
    for (std::string &currentFile: filePaths) {
      *imagesIterator = cv::imread(currentFile);
      imagesIterator++;
    }
  }
 public:
  ImageProcessing(const std::string directory_path,
                  const int number_workers,
                  const int number_files_replica
  )
      :
      numberFilesReplica(number_files_replica),
      numberWorkers(number_workers) {

    //! store the complete path of the file in fileNames
    for (const auto &entry : std::filesystem::directory_iterator(directory_path)) {
      //! replicate the work
      for (int i = 0; i < number_files_replica; i++) {
        filePaths.emplace_back(entry.path());
      }
    }

    //! randomly shuffle the paths of the files
    std::shuffle(std::begin(filePaths), std::end(filePaths), std::default_random_engine{});
  }

  void SequentialTwoStage() {

    //! declare image containers
    std::vector<cv::Mat> destination(filePaths.size());
    std::vector<cv::Mat> images(filePaths.size());
    loadImagesFromDisk(images);

    auto destinationIterator = destination.begin();
    auto startGaussian = std::chrono::high_resolution_clock::now();

    //! apply GaussianBlur
    for (cv::Mat &currentImage: images) {
      cv::GaussianBlur(currentImage, *destinationIterator, cv::Size(3, 3), 0, 0, cv::BORDER_DEFAULT);
      destinationIterator++;
    }
    auto endGaussian = std::chrono::system_clock::now();
    std::cout << "Sequential time Gaussian bluer: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(endGaussian - startGaussian).count() << "ms"
              << std::endl;

    auto startSobel = std::chrono::high_resolution_clock::now();

    //! apply Sobel filter
    for (cv::Mat &currentImage: destination) {
      cv::Sobel(currentImage, currentImage, -1, 1, 1);
    }
    auto endSobel = std::chrono::system_clock::now();

    std::cout << "Sequential time Sobel bluer: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(endSobel - startSobel).count() << "ms"
              << std::endl;

  }

  void SequentialOneStage() {

    //! declare image container
    std::vector<cv::Mat> images(filePaths.size());
    auto start = std::chrono::high_resolution_clock::now();

    loadImagesFromDisk(images);

    //! apply filters
    for (cv::Mat &currentImage: images) {
      cv::GaussianBlur(currentImage, currentImage, cv::Size(3, 3), 0, 0, cv::BORDER_DEFAULT);
      cv::Sobel(currentImage, currentImage, -1, 1, 1);
    }
    auto end = std::chrono::system_clock::now();
    std::cout << "Sequential time one stage took: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms"
              << std::endl;

  }

  static void WorkerComputation(const std::vector<std::string> *filePaths, const int indexStart, const int indexEnd) {

    std::vector<cv::Mat> threadImages((indexEnd - indexStart) + 1);
    auto threadImagesIterator = threadImages.begin();
    for (int currentPath = indexStart; currentPath <= indexEnd; currentPath++) {
      *threadImagesIterator = cv::imread((*filePaths)[currentPath]);
      cv::GaussianBlur(*threadImagesIterator, *threadImagesIterator, cv::Size(3, 3), 0, 0, cv::BORDER_DEFAULT);
      cv::Sobel(*threadImagesIterator, *threadImagesIterator, -1, 1, 1);
      threadImagesIterator++;
    }
  }

  void ParallelStandard() {

    //! declare image container
    std::vector<cv::Mat> images(filePaths.size());

    //! items to compute
    int workerCellsToCompute = std::floor((float) filePaths.size() / (float) numberWorkers);
    int remainedElements = (filePaths.size() % numberWorkers);
    std::vector<std::thread> computationThreads(numberWorkers);
    int currentIndex = 0;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < numberWorkers; i++) {

      int endIndex = currentIndex + (workerCellsToCompute - 1);

      //! consider extra elements
      if (remainedElements) {
        endIndex++;
        remainedElements--;
      }

      computationThreads.emplace_back(std::thread(
          &ImageProcessing::WorkerComputation,
          &filePaths,
          currentIndex,
          endIndex));
      currentIndex = endIndex + 1;
    }

    //! wait the threads
    for (std::thread &currentThread : computationThreads) {
      if (currentThread.joinable()) {
        currentThread.join();
      }
    }

    auto end = std::chrono::system_clock::now();
    std::cout << "Standard thread with: " << numberWorkers << " workers took " << " "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms" << std::endl;
  }

  void ParallelOpenMP() {

    //! declare image container
    std::vector<cv::Mat> images(filePaths.size());

    auto start = std::chrono::high_resolution_clock::now();

    //! apply filters
#pragma omp parallel for num_threads(numberWorkers) schedule(dynamic, 1)
    for (int i = 0; i < filePaths.size(); i++) {
      images[i] = cv::imread(filePaths[i]);
      cv::GaussianBlur(images[i], images[i], cv::Size(3, 3), 0, 0, cv::BORDER_DEFAULT);
      cv::Sobel(images[i], images[i], -1, 1, 1);
    }
    auto end = std::chrono::system_clock::now();
    std::cout << "OpenMP with " << numberWorkers << " workers took " <<
              std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms" << std::endl;
  }
};
#endif //FOURTHASSIGNMENTPARALLEL__IMAGEPROCESSING_H_