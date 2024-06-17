#pragma once

#include "EigenDef.h"

#include <type_traits>

namespace pgo
{

using Vec2d = EigenSupport::V2d;
using Vec3d = EigenSupport::V3d;
using Vec4d = EigenSupport::V4d;

using Vec2i = EigenSupport::V2i;
using Vec3i = EigenSupport::V3i;
using Vec4i = EigenSupport::V4i;

using Mat3d = EigenSupport::M3d;
using Mat4d = EigenSupport::M4d;

static_assert(sizeof(Vec2d) == sizeof(double) * 2, "Size problem");
static_assert(sizeof(Vec3d) == sizeof(double) * 3, "Size problem");
static_assert(sizeof(Vec4d) == sizeof(double) * 4, "Size problem");

static_assert(sizeof(Vec2i) == sizeof(int) * 2, "Size problem");
static_assert(sizeof(Vec3i) == sizeof(int) * 3, "Size problem");
static_assert(sizeof(Vec4i) == sizeof(int) * 4, "Size problem");

static_assert(sizeof(Mat3d) == sizeof(double) * 9, "Size problem");
static_assert(sizeof(Mat4d) == sizeof(double) * 16, "Size problem");

inline Vec3d asVec3d(const double &v)
{
  return Vec3d(v, v, v);
}

inline Vec3d asVec3d(const double v[3])
{
  return Vec3d(v[0], v[1], v[2]);
}

inline Mat3d asMat3d(const double &v)
{
  Mat3d A;
  A.setZero();
  A(0, 0) = v;
  A(1, 1) = v;
  A(2, 2) = v;

  return A;
}

inline Mat3d asMat3d(const double &a0, const double &a1, const double &a2,
  const double &a3, const double &a4, const double &a5,
  const double &a6, const double &a7, const double &a8)
{
  Mat3d A;
  A << a0, a1, a2, a3, a4, a5, a6, a7, a8;

  return A;
}

inline Mat3d asMat3d(const Vec3d &r0, const Vec3d &r1, const Vec3d &r2)
{
  Mat3d A;
  A.row(0) = r0;
  A.row(1) = r1;
  A.row(2) = r2;

  return A;
}

inline Mat3d asMat3d(const double a[9])
{
  Mat3d A;
  A << a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8];

  return A;
}

inline int getInvertedIndex(const Vec3i &idx, int value)
{
  if (value == idx[0])
    return 0;
  if (value == idx[1])
    return 1;
  if (value == idx[2])
    return 2;
  return -1;
}

inline int getInvertedIndex(const Vec4i &idx, int value)
{
  if (value == idx[0])
    return 0;
  if (value == idx[1])
    return 1;
  if (value == idx[2])
    return 2;
  if (value == idx[3])
    return 3;
  return -1;
}

inline bool operator<(const Vec3i &vec1, const Vec3i &vec2)
{
  if (vec1[0] < vec2[0])
    return true;

  if (vec1[0] > vec2[0])
    return false;

  if (vec1[1] < vec2[1])
    return true;

  if (vec1[1] > vec2[1])
    return false;

  return vec1[2] < vec2[2];
}

inline bool operator<(const Vec3d &vec1, const Vec3d &vec2)
{
  if (vec1[0] < vec2[0])
    return true;
  if (vec1[0] > vec2[0])
    return false;
  if (vec1[1] < vec2[1])
    return true;
  if (vec1[1] > vec2[1])
    return false;
  return vec1[2] < vec2[2];
}

inline bool operator<(const Vec4i &vec1, const Vec4i &vec2)
{
  if (vec1[0] < vec2[0])
    return true;
  if (vec1[0] > vec2[0])
    return false;
  if (vec1[1] < vec2[1])
    return true;
  if (vec1[1] > vec2[1])
    return false;
  if (vec1[2] < vec2[2])
    return true;
  if (vec1[2] > vec2[2])
    return false;
  return vec1[3] < vec2[3];
}

inline bool overlap(const Vec3i &vec1, const Vec3i &vec2)
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      if (vec1[i] == vec2[j])
        return true;

  return false;
}

template<typename T, typename U>
constexpr bool is_decay_same = std::is_same_v<std::decay_t<T>, U>;

template<typename Vec1, typename Vec2>
  requires((is_decay_same<Vec1, Vec3i> || is_decay_same<Vec1, Vec4i>) &&
    (is_decay_same<Vec2, Vec3i> || is_decay_same<Vec2, Vec4i>))
inline int firstNonOverlapIndex(const Vec1 &vec1, const Vec2 &vec2)
{
  for (int i = 0; i < int(vec1.size()); i++) {
    bool found = false;
    for (int j = 0; j < int(vec2.size()); j++) {
      if (vec1(i) == vec2(j)) {
        found = true;
        break;
      }
    }
    if (!found)
      return i;
  }

  return -1;
}

inline void rotate(Vec4i &v, int newStartIndex)
{
  // rotate left one time, (1, 2, 3, 0)
  if (newStartIndex == 1) {
    int tmp = v[0];
    v[0] = v[1];
    v[1] = v[2];
    v[2] = v[3];
    v[3] = tmp;
  }
  // rotate left two times (2, 3, 0, 1)
  else if (newStartIndex == 2) {
    std::swap(v[0], v[2]);
    std::swap(v[1], v[3]);
  }
  // rotate right one time (3, 0, 1, 2)
  else if (newStartIndex == 3) {
    int tmp = v[3];
    v[3] = v[2];
    v[2] = v[1];
    v[1] = v[0];
    v[0] = tmp;
  }
}

}  // namespace pgo