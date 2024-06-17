#include "basicIO.h"

#include "pgoLogging.h"

#include <vector>

namespace pgo::BasicIO
{
template<typename T>
int testIO_T()
{
  std::vector<T> testAry{ 1, 2, 3, 4, 5, 0, 1, 2 };
  int ret = write1DText("1.txt", testAry.begin(), testAry.end());

  std::vector<T> testAry1;
  ret = read1DText("1.txt", std::back_inserter(testAry1));
  PGO_ALOG(ret == 0);

  PGO_ALOG(testAry1.size() == testAry.size());
  PGO_ALOG(std::memcmp(testAry.data(), testAry1.data(), testAry1.size() * sizeof(T)) == 0);

  ret = write1DText("2.txt", testAry.data(), testAry.data() + testAry.size());
  PGO_ALOG(ret == 0);

  std::vector<T> testAry2(testAry.size());
  ret = read1DText("2.txt", testAry2.data(), testAry2.size());
  PGO_ALOG(ret == 0);

  PGO_ALOG(std::memcmp(testAry.data(), testAry2.data(), testAry2.size() * sizeof(T)) == 0);
  return 0;
}
}  // namespace pgo::BasicIO

int pgo::BasicIO::testIO()
{
  testIO_T<int>();

  return 0;
}