/**
  Mergesort: sort an array of N integer in parallel using the DAC pattern
  Author: Tiziano De Matteis <dematteis@di.unipi.it>
 */
#include <iostream>
#include <functional>
#include <vector>
#include <algorithm>
#include <cstring>
#include <chrono>
#include "divideAndConquer.cpp"

#define CUTOFF 125000000
#define BASE_CONDITION_CHECK 500


/* ------------------------------------------------------------------ */
//maximum value for arrays elements
const int MAX_NUM=99999;

static bool isArraySorted(int *a, int n) {
  for(int i=1;i<n;i++)
    if(a[i]<a[i-1])
      return false;
  return true;
}

static int *generateRandomArray(int n) {
  srand ((time(0)));
  int *numbers=new int[n];
  for(int i=0;i<n;i++)
    numbers[i]=(int) (rand()) / ((RAND_MAX/MAX_NUM));
  return numbers;
}
/* ------------------------------------------------------------------ */


// Operand and Results share the same format
struct ops{
  int *array=nullptr;		    //array (to sort/sorted)
  int left=0;                     //left index
  int right=0;                    //right index
};

typedef struct ops Operand;
typedef struct ops Result;


/*
 * The divide simply 'split' the array in two: the splitting is only logical.
 * The recursion occur on the left and on the right part
 * UPDATED
 */
std::vector<Operand> divide(Operand op)
{
  std::vector<Operand> result;
  int mid=(op.left+op.right)/2;
  Operand a;
  a.array=op.array;
  a.left=op.left;
  a.right=mid;
  result.push_back(a);

  Operand b;
  b.array=op.array;
  b.left=mid+1;
  b.right=op.right;
  result.push_back(b);
  return result;
}


/*
 * For the base case we resort to std::sort
 * UPDATED
 */
Operand seq(Operand op)
{
  std::sort(&(op.array[op.left]),&(op.array[op.right+1]));

  //the result is essentially the same of the operand
  return op;
}


/*
 * The Merge (Combine) function start from two ordered sub array and construct the original one
 * It uses additional memory
 * UPDATED
 */
Operand mergeMS(std::vector<Operand>ress)
{
  Operand ret;
  //compute what is needed: array pointer, mid, ...
  int *a=ress[0].array;                 //get the array
  int mid=ress[0].right;                //by construction
  int left=ress[0].left, right=ress[1].right;
  int size=right-left+1;
  int *tmp=new int[size];
  int i=left,j=mid+1;

  //merge in order
  for(int k=0;k<size;k++)
  {
    if(i<=mid && (j>right || a[i]<= a[j]))
    {
      tmp[k]=a[i];
      i++;
    }
    else
    {
      tmp[k]=a[j];
      j++;
    }
  }

  //copy back
  memcpy(a+left,tmp,size*sizeof(int));

  delete[] tmp;

  //build the result
  ret.array=a;
  ret.left=left;
  ret.right=right;
  return ret;
}


/*
 * Base case condition
 * UPDATED
 */
bool cond(Operand op)
{
  return (op.right-op.left<=BASE_CONDITION_CHECK);
}

bool checkCutOff(Operand op){
  return (op.right-op.left<=CUTOFF );
}


int main(int argc, char *argv[])
{
  if(argc<2)
  {
    std::cerr << "Usage: "<<argv[0]<< " <num_elements> <num_workers>"<<std::endl;
    exit(-1);
  }
  std::function<std::vector<Operand>(Operand)> div(divide);
  std::function <Operand(Operand)> sq(seq);
  std::function <Operand(std::vector<Operand>)> mergef(mergeMS);
  std::function<bool(Operand)> cf(cond);
  std::function<bool(Operand)> cs(checkCutOff);

  int num_elem=atoi(argv[1]);
  int nwork=atoi(argv[2]);
  //generate a random array
  int *numbers=generateRandomArray(num_elem);

  //build the operand
  Operand op;
  op.array=numbers;
  op.left=0;
  op.right=num_elem-1;
  Result res;

  auto start = std::chrono::high_resolution_clock::now();
  if(nwork == 1)
    // sequential version
    res = sequentialDc(op, cond ,seq ,divide ,mergeMS);
  else
  {
    res = parallelDc(op, cond, seq, divide, mergeMS, checkCutOff);
  }
  /*
   for(int i=0; i< num_elem; i++){
    std::cout << *(res.array+i)<< " ";
  }
  std::cout << std::endl;
   */

  auto elapsed = std::chrono::high_resolution_clock::now() - start;
  auto usec    = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

  if(!isArraySorted(numbers,num_elem))
  {
    fprintf(stderr,"Error: array is not sorted!!\n");
    exit(-1);
  }
  printf("Time (usecs): %ld\n",usec);

  return 0;
}