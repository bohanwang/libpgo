/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#pragma once

namespace pgo
{
namespace SolidDeformationModel
{
enum FD_TEST_TYPE : int
{
  FD_TT_ELEMENT = 1,
  FD_TT_ASSEM_E = 2,
  FD_TT_ASSEM_EFP = 4,
};
int fdTestTetMesh(const char *tetMeshFilename, int numTestDOFs, int fdTestFlag);
int fdTestShellMesh(const char *surfaceMeshFilename, int numTestDOFs, int fdTestFlag);
//int fdTestBending(const char *surfaceMeshFilename, const char *basisFilename, const char *fiberDirectionFilename, int numDOFs, int testType);

}  // namespace SolidDeformationModel
}  // namespace pgo