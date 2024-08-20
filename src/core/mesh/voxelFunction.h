#ifndef GRIDFUNCTION_H
#define GRIDFUNCTION_H

#include "boundingBox.h"
#include "vec3d.h"
#include "vec3i.h"
#include <vector>

// a voxel grid
class VoxelGrid
{
public:
  using packedIndex = std::size_t;
  VoxelGrid() {}

  inline void initialize(const Vec3d & bmin, const Vec3d & bmax, const Vec3i & resolution);

  // get #voxels on each axis
  inline const Vec3i & getResolution() const { return resolution; }

  inline Vec3i voxelIndex(const Vec3d & position) const;
  inline packedIndex packedVoxelIndex(const Vec3d & position) const;

  // pack voxel index
  // note: don't convert between packed voxel indices and packed grid indices directly
  inline packedIndex packedVoxelIndex(int vi, int vj, int vk) const;
  inline packedIndex packedVoxelIndex(const Vec3i & voxelIndex) const;

  inline void unpackVoxelIndex(packedIndex packedVoxelIndex, int & vi, int & vj, int & vk) const;
  inline Vec3i unpackVoxelIndex(packedIndex packedVoxelIndex) const;

  inline BoundingBox getVoxelBoundingBox(packedIndex packedVoxelIndex) const;

  inline bool voxelInsideBox(int i, int j, int k) const; // whether voxel index is inside the voxel grid
  inline bool voxelInsideBox(const Vec3i & voxelIndex) const; // whether voxel index is inside the voxel grid
  inline bool voxelInsideBox(packedIndex packedVoxelIndex) const { return packedVoxelIndex < numVoxels; }
  inline packedIndex getNumVoxels() const { return numVoxels; }
  inline Vec3d getVoxelCenter(packedIndex packedVoxelIndex) const;

  inline double getVoxelDiagonal() const { return voxelDiagonal; }

  inline const Vec3d & getVoxelSides() const { return voxelSize; }

  // pack grid index
  // note: don't convert between packed voxel indices and packed grid indices directly
  inline packedIndex packedGridIndex(int gi, int gj, int gk) const;
  inline packedIndex packedGridIndex(const Vec3i & gridIndex) const;
  inline void unpackGridIndex(packedIndex packedGridIndex, int & gi, int & gj, int & gk) const;
  inline Vec3i unpackGridIndex(packedIndex packedGridIndex) const;

  // given grid point index, get the position of this index
  inline Vec3d getGridPosition(int gi, int gj, int gk) const;
  inline Vec3d getGridPosition(const Vec3i & gridIndex) const { return getGridPosition(gridIndex[0], gridIndex[1], gridIndex[2]); }
  inline Vec3d getGridPosition(const packedIndex packedGridIndex) const { Vec3i g = unpackGridIndex(packedGridIndex); return getGridPosition(g); }

  inline bool gridInsideBox(int gi, int gj, int gk) const;
  inline bool gridInsideBox(const Vec3i & gridIndex) const;
  inline bool gridInsideBox(packedIndex packedGridIndex) const { return packedGridIndex < numGridPoints; }
  inline packedIndex getNumGrids() const { return numGridPoints; }

protected:
  Vec3i resolution{0}; // #voxels on each axis
  Vec3i gridResolution{0};
  BoundingBox bb;
  Vec3d bbside{0.0};
  Vec3d voxelSize{0.0}, invVoxelSize{0.0};
  Vec3d halfVoxelSize{0.0};
  double voxelDiagonal = 0.0;
  packedIndex numVoxels = 0;
  packedIndex numGridPoints = 0;
};

// a voxel grid
template<typename T>
class VoxelFunction : public VoxelGrid
{
public:
  VoxelFunction() {}

  void initialize(const Vec3d & bmin, const Vec3d & bmax, const Vec3i & resolution, const T & initValue);

  // assuming this voxelIndex is inside the box
  typename std::vector<T>::reference value(const Vec3i & voxelIndex)  { return buffer[packedVoxelIndex(voxelIndex)]; }
  typename std::vector<T>::const_reference value(const Vec3i & voxelIndex) const { return buffer[packedVoxelIndex(voxelIndex)]; }

  typename std::vector<T>::reference value(packedIndex packedVoxelIndex)  { return buffer[packedVoxelIndex]; }
  typename std::vector<T>::const_reference value(packedIndex packedVoxelIndex) const { return buffer[packedVoxelIndex]; }

  void setAllValues(const T & value) { buffer.assign(buffer.size(), value); }

protected:
  std::vector<T> buffer;
};

//=======================================================================
//                    Below are implementations
//=======================================================================

///////////////////////////////////////////////////////////////
//                       VoxelGrid
///////////////////////////////////////////////////////////////

inline void VoxelGrid::initialize(const Vec3d & bmin, const Vec3d & bmax, const Vec3i & resolution)
{
  this->resolution = resolution;
  gridResolution = Vec3i(resolution[0]+1, resolution[1]+1, resolution[2]+1);
  bb = BoundingBox(bmin, bmax);
  bbside = bmax - bmin;
  for(int d = 0; d < 3; d++)
  {
    voxelSize[d] = bbside[d] / resolution[d];
    halfVoxelSize[d] = voxelSize[d] * 0.5;
    invVoxelSize[d] = 1.0 / voxelSize[d];
  }
  voxelDiagonal = len(voxelSize);

  numVoxels = (packedIndex)resolution[0] * resolution[1] * resolution[2];
  numGridPoints = (packedIndex)gridResolution[0] * gridResolution[1] * gridResolution[2];
}

inline Vec3i VoxelGrid::voxelIndex(const Vec3d & pos) const
{
  Vec3i index;
  index[0] = (int)((pos[0] - bb.bmin()[0]) * invVoxelSize[0]);
  index[1] = (int)((pos[1] - bb.bmin()[1]) * invVoxelSize[1]);
  index[2] = (int)((pos[2] - bb.bmin()[2]) * invVoxelSize[2]);
  return index;
}

inline VoxelGrid::packedIndex VoxelGrid::packedVoxelIndex(const Vec3d & position) const
{
  Vec3i index = voxelIndex(position);
  return packedVoxelIndex(index);
}

inline VoxelGrid::packedIndex VoxelGrid::packedVoxelIndex(const Vec3i & index) const
{
  return ((packedIndex)index[2] * resolution[1] + index[1]) * resolution[0] + index[0];
}

inline VoxelGrid::packedIndex VoxelGrid::packedVoxelIndex(int i, int j, int k) const
{
  return ((packedIndex)k * resolution[1] + j) * resolution[0] + i;
}

inline void VoxelGrid::unpackVoxelIndex(packedIndex packedIndex, int & i, int & j, int & k) const
{
  i = packedIndex % resolution[0];
  packedIndex = packedIndex / resolution[0];
  j = packedIndex % resolution[1];
  k = packedIndex / resolution[1];
}

inline Vec3i VoxelGrid::unpackVoxelIndex(packedIndex packedIndex) const
{
  Vec3i v;
  unpackVoxelIndex(packedIndex, v[0], v[1], v[2]);
  return v;
}

inline Vec3d VoxelGrid::getGridPosition(int i, int j, int k) const
{
  return Vec3d(bb.bmin()[0] + i * voxelSize[0], bb.bmin()[1] + j * voxelSize[1], bb.bmin()[2] + k * voxelSize[2]);
}

inline VoxelGrid::packedIndex VoxelGrid::packedGridIndex(int gi, int gj, int gk) const
{
  return ((packedIndex)gk * gridResolution[1] + gj) * gridResolution[0] + gi;
}

inline VoxelGrid::packedIndex VoxelGrid::packedGridIndex(const Vec3i & gridIndex) const
{
  return ((packedIndex)gridIndex[2] * gridResolution[1] + gridIndex[1]) * gridResolution[0] + gridIndex[0];
}

inline void VoxelGrid::unpackGridIndex(packedIndex packedGridIndex, int & gi, int & gj, int & gk) const
{
  gi = packedGridIndex % gridResolution[0];
  packedGridIndex = packedGridIndex / gridResolution[0];
  gj = packedGridIndex % gridResolution[1];
  gk = packedGridIndex / gridResolution[1];
}

inline Vec3i VoxelGrid::unpackGridIndex(packedIndex packedGridIndex) const
{
  Vec3i g;
  unpackGridIndex(packedGridIndex, g[0], g[1], g[2]);
  return g;
}

inline BoundingBox VoxelGrid::getVoxelBoundingBox(packedIndex packedIndex) const
{
  Vec3i v;
  unpackVoxelIndex(packedIndex, v[0], v[1], v[2]);
  return BoundingBox(getGridPosition(v[0], v[1], v[2]), getGridPosition(v[0]+1,v[1]+1,v[2]+1));
}

inline bool VoxelGrid::voxelInsideBox(const Vec3i & index) const
{
  return index[0] >= 0 && index[0] < resolution[0] && index[1] >= 0 && index[1] < resolution[1] && index[2] >= 0 && index[2] < resolution[2];
}

inline bool VoxelGrid::voxelInsideBox(int i, int j, int k) const
{
  return i >= 0 && i < resolution[0] && j >= 0 && j < resolution[1] && k >= 0 && k < resolution[2];
}

inline Vec3d VoxelGrid::getVoxelCenter(packedIndex packedIndex) const
{
  int i, j, k;
  unpackVoxelIndex(packedIndex, i, j, k);
  return getGridPosition(i, j, k) + halfVoxelSize;
}

inline bool VoxelGrid::gridInsideBox(int gi, int gj, int gk) const
{
  return gi >= 0 && gi < gridResolution[0] && gj >= 0 && gj < gridResolution[1] && gk >= 0 && gk < gridResolution[2];
}

inline bool VoxelGrid::gridInsideBox(const Vec3i & gridIndex) const
{
  return gridIndex[0] >= 0 && gridIndex[0] < gridResolution[0] &&
      gridIndex[1] >= 0 && gridIndex[1] < gridResolution[1] &&
      gridIndex[2] >= 0 && gridIndex[2] < gridResolution[2];
}

///////////////////////////////////////////////////////////////
//                     VoxelFunction
///////////////////////////////////////////////////////////////
template<typename T>
void VoxelFunction<T>::initialize(const Vec3d & bmin, const Vec3d & bmax, const Vec3i & resolution, const T & initValue)
{
  VoxelGrid::initialize(bmin, bmax, resolution);
  buffer.assign(numVoxels, initValue);
}

#endif
