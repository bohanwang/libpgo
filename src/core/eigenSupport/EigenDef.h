/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#if defined(PGO_HAS_MKL)
#  include <Eigen/PardisoSupport>
#endif

namespace pgo
{
namespace EigenSupport
{
typedef Eigen::Index IDX;

typedef Eigen::Triplet<double> TripletD;
typedef Eigen::Triplet<Eigen::Index> TripletI;

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMatD;
typedef Eigen::SparseMatrix<Eigen::Index, Eigen::RowMajor> SpMatI;

typedef Eigen::Matrix<double, 1, 1> M1d;

typedef Eigen::Vector2d V2d;
typedef Eigen::Vector2i V2i;

typedef Eigen::Vector3d V3d;
typedef Eigen::Vector3i V3i;

typedef Eigen::Vector4d V4d;
typedef Eigen::Vector4i V4i;

typedef Eigen::Matrix<double, 6, 1> V6d;
typedef Eigen::Matrix<int, 6, 1> V6i;

typedef Eigen::Matrix<double, 7, 1> V7d;
typedef Eigen::Matrix<int, 7, 1> V7i;

typedef Eigen::Matrix<double, 8, 1> V8d;
typedef Eigen::Matrix<int, 8, 1> V8i;

typedef Eigen::Matrix<double, 9, 1> V9d;
typedef Eigen::Matrix<int, 9, 1> V9i;

typedef Eigen::Matrix<double, 10, 1> V10d;
typedef Eigen::Matrix<int, 10, 1> V10i;

typedef Eigen::Matrix<double, 12, 1> V12d;
typedef Eigen::Matrix<int, 12, 1> V12i;

typedef Eigen::Matrix<double, 18, 1> V18d;
typedef Eigen::Matrix<int, 18, 1> V18i;

typedef Eigen::Matrix2d M2d;
typedef Eigen::Matrix2i M2i;

typedef Eigen::Matrix3d M3d;
typedef Eigen::Matrix3i M3i;

typedef Eigen::Matrix4d M4d;
typedef Eigen::Matrix4i M4i;

typedef Eigen::Matrix<double, 6, 6> M6d;
typedef Eigen::Matrix<int, 6, 6> M6i;

typedef Eigen::Matrix<double, 7, 7> M7d;
typedef Eigen::Matrix<int, 7, 7> M7i;

typedef Eigen::Matrix<double, 9, 9> M9d;
typedef Eigen::Matrix<int, 9, 9> M9i;

typedef Eigen::Matrix<double, 10, 10> M10d;
typedef Eigen::Matrix<int, 10, 10> M10i;

typedef Eigen::Matrix<double, 18, 18> M18d;
typedef Eigen::Matrix<double, 24, 24> M24d;
typedef Eigen::Matrix<int, 18, 18> M18i;

typedef Eigen::Matrix<double, 12, 12> M12d;
typedef Eigen::Matrix<double, 6, 9> M6x9d;

typedef Eigen::Matrix<double, 6, 18> M6x18d;
typedef Eigen::Matrix<double, 18, 6> M18x6d;
typedef Eigen::Matrix<double, 3, 18> M3x18d;

typedef Eigen::Matrix<double, 2, 3> M2x3d;
typedef Eigen::Matrix<double, 3, 2> M3x2d;
typedef Eigen::Matrix<double, 9, 12> M9x12d;
typedef Eigen::Matrix<double, 12, 9> M12x9d;
typedef Eigen::Matrix<double, 12, 12> M12x12d;
typedef Eigen::Matrix<int, 12, 12> M12i;

typedef Eigen::Matrix<double, 3, 2> M3x2d;
typedef Eigen::Matrix<double, 2, 3> M2x3d;

typedef Eigen::Matrix<double, 4, 3> M4x3d;
typedef Eigen::Matrix<double, 4, 9> M4x9d;
typedef Eigen::Matrix<double, 4, 18> M4x18d;

typedef Eigen::Matrix<double, 6, 9> M6x9d;
typedef Eigen::Matrix<double, 9, 6> M9x6d;

typedef Eigen::Matrix<double, 9, 12> M9x12d;
typedef Eigen::Matrix<double, 12, 9> M12x9d;

typedef Eigen::VectorXd VXd;
typedef Eigen::VectorXi VXi;

typedef Eigen::MatrixXd MXd;
typedef Eigen::MatrixXi MXi;

typedef Eigen::Matrix3Xd M3Xd;
typedef Eigen::MatrixX3d MX3d;

typedef Eigen::Affine2d A2d;
typedef Eigen::Affine3d A3d;

typedef Eigen::Quaterniond Qd;

// float
typedef Eigen::Vector2f V2f;
typedef Eigen::Vector3f V3f;
typedef Eigen::Vector4f V4f;
typedef Eigen::Matrix2f M2f;
typedef Eigen::Matrix3f M3f;
typedef Eigen::Matrix4f M4f;

template<typename T>
using EigenArray = std::vector<T, typename Eigen::aligned_allocator<T>>;

template<typename T, int MapOptions = Eigen::Unaligned, typename StrideType = Eigen::Stride<0, 0>>
using Mp = Eigen::Map<T, MapOptions, StrideType>;

#if defined(PGO_HAS_MKL)
typedef Eigen::PardisoLLT<SpMatD> SPDSolver;
typedef Eigen::PardisoLDLT<SpMatD> SymSolver;
typedef Eigen::PardisoLU<SpMatD> LUSolver;
#else
typedef Eigen::SimplicialLLT<SpMatD> SPDSolver;
typedef Eigen::SimplicialLDLT<SpMatD> SymSolver;
typedef Eigen::SparseLU<SpMatD> LUSolver;
#endif

typedef const Eigen::Ref<const Eigen::VectorXd> ConstRefVecXd;
typedef const Eigen::Ref<const Eigen::MatrixXd> ConstRefMatXd;

typedef Eigen::Ref<Eigen::VectorXd> RefVecXd;
typedef Eigen::Ref<Eigen::MatrixXd> RefMatXd;

typedef Eigen::Ref<Eigen::VectorXi> RefVecXi;

// from boost
inline void hashCombine(std::size_t &h, std::size_t k)
{
  const std::size_t m = UINT64_C(0xc6a4a7935bd1e995);
  const int r = 47;

  k *= m;
  k ^= k >> r;
  k *= m;

  h ^= k;
  h *= m;

  // Completely arbitrary number, to prevent 0's
  // from hashing to 0.
  h += 0xe6546b64;
}

struct IntPairHash
{
  std::size_t operator()(const std::pair<int, int> &pr) const
  {
    std::size_t h = 0;
    hashCombine(h, std::hash<int>()(pr.first));
    hashCombine(h, std::hash<int>()(pr.second));

    return h;
  }
};

struct IntPairEqual
{
  bool operator()(const std::pair<int, int> &pr1, const std::pair<int, int> &pr2) const
  {
    return pr1.first == pr2.first && pr1.second == pr2.second;
  }
};

template<int N>
struct IntArrayHash
{
  std::size_t operator()(const std::array<int, N> &val) const
  {
    std::size_t h = 0;
    for (int i = 0; i < N; i++) {
      hashCombine(h, std::hash<int>()(val[i]));
    }

    return h;
  }
};

template<int N>
struct IntArrayEqual
{
  bool operator()(const std::array<int, N> &v1, const std::array<int, N> &v2) const
  {
    return std::memcmp(v1.data(), v2.data(), sizeof(std::array<int, N>)) == 0;
  }
};

template<int N>
struct IntArrayLess
{
  bool operator()(const std::array<int, N> &v1, const std::array<int, N> &v2) const
  {
    return std::memcmp(v1.data(), v2.data(), sizeof(std::array<int, N>)) < 0;
  }
};

using EntryMap = std::unordered_map<std::pair<int, int>, std::ptrdiff_t, IntPairHash, IntPairEqual>;

struct SymbolicMmData;
}  // namespace EigenSupport
}  // namespace pgo
