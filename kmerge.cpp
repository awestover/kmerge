/*
Author: Alek Westover
TODOs:
- performance testing
- cache locality
- tune

Description:
parallel k-merge

Say that you have k (a power of 2 for now) sorted lists
and want to merge them into a single sorted list using p (a power of 2)
processors.
We make a tree with p leaves, where each vertex in the tree corresponds to
merging two lists. Parallelism across these vertices is easy.
Within each of these pair merges we can also achieve parallelism via merging odd
/ even halves of the array first and then doing some swaps in parallel.
*/

#include <cassert>
#include <iostream>
#include <thread>
#include <vector>

using namespace std;

// helper function: print an array
void printarr(int *A, int n) {
  for (int i = 0; i < n; i++) std::cout << A[i] << " ";
  std::cout << endl;
}

// compare elements start and start+step
// and swap if they are out of order
void enforce_order(int *X, uint start, uint step) {
  if (X[start] > X[start + step]) std::swap(X[start], X[start + step]);
}

// takes an array, namely every step-th element of A starting from offset
// and applies enforce_order to every adjacent pair of elements
// i.e. enforce_order(A, offset+step, step), enforce_order(A, offset+2*step,
// step),
//  ... enforce_order(A, offset+(n-1)*step, step)
// n is the number of elements in the array (which is a sub-array of A)
void neighbor_swaps(int *A, int *B, uint n, uint step, uint offset, uint p) {
  const uint MINSIZE = 256;
  uint num_pairs = (n / p) / 2;
  if (p > 1 && num_pairs >= MINSIZE) {
    std::vector<std::thread> handles;
    for (int i = 0; i < p; ++i) {
      handles.push_back(thread([&, i]() {
        int is_final = i == p - 1;
        for (int j = 0; j < (n / p) / 2 - is_final; ++j) {
          enforce_order(A, offset + (2 * (j + i * num_pairs) + 1) * step, step);
          enforce_order(B, offset + (2 * (j + i * num_pairs) + 1) * step, step);
        }
      }));
    }
    for (std::thread &handle : handles) handle.join();
  } else {  // p=1
    for (int j = 0; j < n / 2 - 1; ++j) {
      enforce_order(A, offset + (j * 2 + 1) * step, step);
      enforce_order(B, offset + (j * 2 + 1) * step, step);
    }
  }
  uint last_a = offset + (n - 1) * step;
  uint first_b = offset;
  if (A[last_a] > B[first_b]) std::swap(A[last_a], B[first_b]);
}

// merge two sorted arrays in serial
// creates a large auxiliary array
// in particular the arrays are
// A[offset], A[offset +step]...A[offset + step*(n-1)]
// and same for B
// at the end this function writes the data back on to A, B
void serial_merge_pair(int *A, int *B, uint n, uint step, uint offset) {
  int *C = (int *)malloc(2 * n * sizeof(int));
  uint a = 0;
  uint b = 0;
  while (a < n && b < n) {
    if (A[offset + step * a] < B[offset + step * b]) {
      C[a + b] = A[offset + step * a];
      a += 1;
    } else {
      C[a + b] = B[offset + step * b];
      b += 1;
    }
  }
  while (a < n) {
    C[a + b] = A[offset + step * a];
    a += 1;
  }
  while (b < n) {
    C[a + b] = B[offset + step * b];
    b += 1;
  }
  for (int i = 0; i < n; ++i) {
    A[offset + step * i] = C[i];
    B[offset + step * i] = C[n + i];
  }
  free(C);
}

// merge two sorted arrays in parallel
// in particular the arrays are
// A[offset], A[offset +step]...A[offset + step*(n-1)]
// and same for B
// at the end this function writes the data back on to A, B
// the function works as follows:
// if enough parallelism available, split into even odd
// merge even odd (recursively, with this method!). then stitch them together
void parallel_merge_pair(int *A, int *B, uint n, uint step, uint offset,
                         uint p) {
  const int THRESHOLD = 256;
  if (p <= 1 || n < THRESHOLD) {
    serial_merge_pair(A, B, n, step, offset);
  } else {
    // merge even odd
    std::thread even_thread([&]() {
      parallel_merge_pair(A, B, n >> 1, step << 1, offset, p >> 1);
    });
    std::thread odd_thread([&]() {
      parallel_merge_pair(A, B, n >> 1, step << 1, offset + step, p >> 1);
    });
    even_thread.join();
    odd_thread.join();

    // stitch results together
    neighbor_swaps(A, B, n, step, offset, p);
  }
}

// assumes num_arrays and p (num_processors) are powers of 2
// Takes a list of sorted arrays, namely
// arrays, arrays+array_size, ... arrays+array_size*(num_arrays-1)
// each of size array_size and merges them
// into the array arrays
void kmerge_pll_tree(int *arrays, uint array_size, uint num_arrays, uint p) {
  const int THRESHOLD = 256;
  if (num_arrays > 2) {
    if (p > 1 || array_size * num_arrays <= THRESHOLD) {
      std::thread left_thread([&]() {
        kmerge_pll_tree(arrays, array_size, num_arrays >> 1, p >> 1);
      });
      std::thread right_thread([&]() {
        kmerge_pll_tree(arrays + array_size * num_arrays / 2, array_size,
                        num_arrays >> 1, p >> 1);
      });
      left_thread.join();
      right_thread.join();
    } else {
      kmerge_pll_tree(arrays, array_size, num_arrays >> 1, p);
      kmerge_pll_tree(arrays + array_size * num_arrays / 2, array_size,
                      num_arrays >> 1, p);
    }
  }
  parallel_merge_pair(arrays, arrays + array_size * num_arrays / 2,
                      array_size * num_arrays / 2, 1, 0, p);
}

// *******TESTS ***********//
bool is_equal(int *A, int *B, int n) {
  for (int i = 0; i < n; i++) {
    if (A[i] != B[i]) return false;
  }
  return true;
}

void test_enforce_order() {
  int X[10] = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
  enforce_order(X, 0, 1);
  int expected1[10] = {8, 9, 7, 6, 5, 4, 3, 2, 1, 0};
  assert(is_equal(X, expected1, 10));
  enforce_order(X, 0, 1);
  int expected2[10] = {8, 9, 7, 6, 5, 4, 3, 2, 1, 0};
  assert(is_equal(X, expected2, 10));
  enforce_order(X, 2, 4);
  int expected3[10] = {8, 9, 3, 6, 5, 4, 7, 2, 1, 0};
  assert(is_equal(X, expected3, 10));
}

void test_neighbor_swaps() {
  int A1[4] = {9, 8, 7, 6};
  int B1[4] = {1, 2, 3, 4};
  int A2[4] = {9, 7, 8, 1};
  int B2[4] = {6, 2, 3, 4};
  neighbor_swaps(A1, B1, 4, 1, 0, 1);
  assert(is_equal(A1, A2, 4));
  assert(is_equal(B1, B2, 4));

  const int N = 4096;
  int *A4 = new int[N];
  int *B4 = new int[N];
  for (int i = 0; i < N; i++) {
    A4[i] = i - (2 * ((i + 1) % 2)) + 1;
    B4[i] = 4 * N;
  }
  A4[0] = 0;
  A4[N - 1] = N - 1;
  int *A5 = new int[N];
  int *B5 = new int[N];
  for (int i = 0; i < N; i++) {
    A5[i] = i;
    B5[i] = 4 * N;
  }
  neighbor_swaps(A4, B4, N, 1, 0, 4);

  assert(is_equal(A4, A5, N));
  assert(is_equal(B4, B5, N));

  delete[] A4;
  delete[] A5;
  delete[] B4;
  delete[] B5;

  int A6[8] = {0, 4, 0, 6, 10, 12, 11, 21};
  int B6[8] = {19, 24, 23, 28, 27, 28, 36, 47};
  neighbor_swaps(A6, B6, 8, 1, 0, 2);

  for (int i = 0; i < 7; i++) {
    assert(A6[i] <= A6[i + 1]);
    assert(B6[i] <= B6[i + 1]);
  }
  assert(A6[7] <= B6[0]);
}

void test_serial_merge() {
  int A[4] = {1, 4, 5, 7};
  int B[4] = {2, 3, 6, 8};
  serial_merge_pair(A, B, 4, 1, 0);
  int A2[4] = {1, 2, 3, 4};
  int B2[4] = {5, 6, 7, 8};
  assert(is_equal(A, A2, 4));
  assert(is_equal(B, B2, 4));
}

void test_parallel_merge() {
  const int N = 1024;
  int *A = new int[N];
  int *B = new int[N];
  A[0] = 0;
  B[0] = 0;
  for (int i = 1; i < N; i++) {
    A[i] = A[i - 1] + 1 + (5 * i) % 7;
    B[i] = B[i - 1] + 1 + (3 * i) % 11;
  }
  parallel_merge_pair(A, B, N, 1, 0, 2);
  for (int i = 0; i < N - 1; i++) {
    assert(A[i] <= A[i + 1]);
    assert(B[i] <= B[i + 1]);
  }
  assert(A[N - 1] <= B[0]);
  delete[] A;
  delete[] B;
}

void test_parallel_merge2() {
  const int N = 1024;
  int *A = new int[N * 2];
  A[0] = 0;
  A[N] = 0;
  for (int i = 1; i < N; i++) {
    A[i] = A[i - 1] + 1 + (5 * i) % 7;
    A[i + N] = A[N + i - 1] + 1 + (3 * i) % 11;
  }
  parallel_merge_pair(A, A + N, N, 1, 0, 2);
  for (int i = 0; i < 2 * N - 1; i++) {
    assert(A[i] <= A[i + 1]);
  }
  delete[] A;
}

void test_kmerge_pll_tree() {
  const int N = 1024;
  int *A = new int[4 * N];
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < N; j++) {
      A[N * i + j] = N * (5 - i) + j;
    }
  }
  kmerge_pll_tree(A, N, 4, 1);
  for (int i = 0; i < 4 * N - 1; i++) {
    assert(A[i] <= A[i + 1]);
  }
  delete[] A;
}

int main() {
  std::cout << "RUNNING TESTS" << endl;
  test_enforce_order();
  test_neighbor_swaps();
  test_serial_merge();
  test_parallel_merge();
  test_parallel_merge2();
  test_kmerge_pll_tree();
  std::cout << "TESTS COMPLETED" << endl;
}
