#pragma once

#include <string_view>
#include <iterator>
#include <type_traits>
#include <fstream>
#include <iostream>

namespace pgo
{
namespace BasicIO
{
template<typename T>
concept IsBasicScalarType = std::is_arithmetic<std::decay_t<T>>::value;

template<class T>
struct output_iterator_traits : std::iterator_traits<T>
{
};

template<class Container>
struct output_iterator_traits<std::back_insert_iterator<Container>>
{
  using value_type = Container::value_type;
};

template<class Container>
struct output_iterator_traits<std::front_insert_iterator<Container>>
{
  using value_type = Container::value_type;
};

template<class Container>
struct output_iterator_traits<std::insert_iterator<Container>>
{
  using value_type = Container::value_type;
};

template<typename Iter>
  requires(std::forward_iterator<Iter> && IsBasicScalarType<typename std::iter_value_t<Iter>>)
int write1DText(const char *filename, Iter itBeg, Iter itEnd)
{
  std::ofstream outfile(filename);
  if (!outfile)
    return 1;

  for (auto it = itBeg; it != itEnd; ++it) {
    outfile << *it << ' ';
  }
  outfile.close();

  return 0;
}

template<typename Iter>
  requires(std::output_iterator<Iter, typename output_iterator_traits<Iter>::value_type> && IsBasicScalarType<typename output_iterator_traits<Iter>::value_type>)
int read1DText(const char *filename, Iter itBeg, size_t maxAllowedCount = 0ull)
{
  std::ifstream infile(filename);
  if (!infile)
    return 1;

  if (maxAllowedCount == 0ull)
    maxAllowedCount = ~0ull;

  typename output_iterator_traits<Iter>::value_type tempScalar;
  size_t count = 0ull;
  while (!infile.eof() && count < maxAllowedCount && infile >> tempScalar) {
    (*itBeg) = tempScalar;
    ++itBeg;
    ++count;
  }
  infile.close();

  return 0;
}

template<typename Iter>
  requires(std::forward_iterator<Iter> && IsBasicScalarType<typename std::iter_value_t<Iter>>)
int write1DBinary(const char *filename, Iter itBeg, Iter itEnd)
{
  std::ofstream outfile(filename, std::ios_base::binary);
  if (!outfile)
    return 1;

  for (auto it = itBeg; it != itEnd; ++it) {
    typename std::iter_value_t<Iter> &refVal = *it;
    outfile.write((const char *)&refVal, sizeof(refVal));
  }

  outfile.close();

  return 0;
}

template<typename Iter>
  requires(std::output_iterator<Iter, typename output_iterator_traits<Iter>::value_type> && IsBasicScalarType<typename output_iterator_traits<Iter>::value_type>)
int read1DBinary(const char *filename, Iter itBeg, size_t maxAllowedCount = 0ull)
{
  std::ifstream infile(filename, std::ios_base::binary);
  if (!infile)
    return 1;

  if (maxAllowedCount == 0ull)
    maxAllowedCount = ~0ull;

  typename output_iterator_traits<Iter>::value_type tempScalar;
  size_t count = 0ull;
  while (!infile.eof() && count < maxAllowedCount && infile.read((char *)tempScalar, sizeof(tempScalar))) {
    (*itBeg) = tempScalar;
    ++itBeg;
    ++count;
  }
  infile.close();

  return 0;
}

template<typename real>
  requires(IsBasicScalarType<real>)
int write2DMatrixToStream(std::ofstream &ofile, int m, int n, const real *matrix)
{
  ofile.write((const char *)matrix, sizeof(real) * m * n);
  if (!ofile)
    return 1;

  return 0;
}

// writes out the m x n matrix header onto the stream
int write2DMatrixHeaderToStream(std::ofstream &ofile, int m, int n)
{
  ofile.write((const char *)&m, sizeof(int));
  if (!ofile)
    return 1;

  ofile.write((const char *)&n, sizeof(int));
  if (!ofile)
    return 1;

  return 0;
}

template<typename real>
  requires(IsBasicScalarType<real>)
int write2DMatrix(const char *filename, int m, int n, const real *matrix)
{
  std::ofstream file(filename, std::ios::binary);
  if (!file) {
    std::cerr << "Error opening the file for writing: " << filename << std::endl;
    return 1;
  }

  if (write2DMatrixHeaderToStream(file, m, n) != 0) {
    std::cerr << "Error writing the matrix header to disk file: " << filename << std::endl;
    return 1;
  }

  if (write2DMatrixToStream(file, m, n, matrix) != 0) {
    std::cerr << "Error writing the matrix to disk file: " << filename << std::endl;
    return 1;
  }

  file.close();

  return 0;
}

// read the m x n matrix from the stream, in binary format
template<typename real>
  requires(IsBasicScalarType<real>)
int read2DMatrixFromStream(std::ifstream &ifile, int M, int N, real *matrix)
{
  ifile.read((char *)matrix, sizeof(real) * M * N);
  if (!ifile) {
    std::cerr << "Error reading the matrix data from stream." << std::endl;
    return 1;
  }

  return 0;
}

int read2DMatrixSizeFromStream(std::ifstream &ifile, int *m, int *n)
{
  ifile.read((char *)m, sizeof(int));
  if (!ifile) {
    std::cerr << "Error reading the matrix row size from stream." << std::endl;
    return 1;
  }

  ifile.read((char *)n, sizeof(int));
  if (!ifile) {
    std::cerr << "Error reading the matrix column size from stream." << std::endl;
    return 1;
  }

  return 0;
}

template<typename real>
  requires(IsBasicScalarType<real>)
int readMatrixFromDisk(const char *filename, int *m, int *n, real *matrix)
{
  std::ifstream infile(filename, std::ios::binary);
  if (!infile) {
    std::cerr << "Error opening the disk file: " << filename << std::endl;
    return 1;
  }

  if (read2DMatrixSizeFromStream(infile, m, n) != 0) {
    std::cerr << "Error reading matrix header from disk file: " << filename << std::endl;
    return 1;
  }

  if (matrix) {
    if (read2DMatrixFromStream(infile, *m, *n, matrix) != 0) {
      std::cerr << "Error reading matrix data from disk file: " << filename << std::endl;
      return 1;
    }
  }

  return 0;
}

int testIO();

}  // namespace BasicIO
}  // namespace pgo