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

int testIO();

}  // namespace BasicIO
}  // namespace pgo