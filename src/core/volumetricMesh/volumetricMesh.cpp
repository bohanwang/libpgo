/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "volumetricMesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Yijing Li                                *
 * http://www.jernejbarbic.com/vega                                      *
 *                                                                       *
 * Research: Jernej Barbic, Hongyi Xu, Yijing Li,                        *
 *           Danyong Zhao, Bohan Wang,                                   *
 *           Fun Shing Sin, Daniel Schroeder,                            *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC,                *
 *          Sloan Foundation, Okawa Foundation,                          *
 *          USC Annenberg Foundation                                     *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#include "volumetricMeshParser.h"
#include "volumetricMesh.h"
#include "volumetricMeshENuMaterial.h"
#include "volumetricMeshOrthotropicMaterial.h"
#include "volumetricMeshMooneyRivlinMaterial.h"

#include "range.h"
#include "stringHelper.h"
#include "pgoLogging.h"

#include <tbb/parallel_for.h>

#include <cfloat>
#include <cstring>
#include <cassert>
#include <iostream>
#include <map>
#include <numeric>

namespace pgo
{
namespace VolumetricMeshes
{
// parses the mesh, and returns the string corresponding to the element type
VolumetricMesh::VolumetricMesh(const char *filename, fileFormatType fileFormat, int numElementVertices_, elementType *elementType_, int verbose):
  numElementVertices(numElementVertices_)
{
  if (verbose) {
    printf("Opening file %s.\n", filename);
    fflush(nullptr);
  }

  if (fileFormat == VolumetricMesh::BY_EXT) {
    fileFormat = VolumetricMesh::getFileFormatTypeByExt(filename);
    // follow the bahavior of VolumetricMesh::getElementType(), which assumes ASCII if unknown file format found
    if (fileFormat == VolumetricMesh::NUM_FILE_FORMATS)
      fileFormat = VolumetricMesh::ASCII;
  }

  switch (fileFormat) {
  case ASCII:
    loadFromAscii(filename, elementType_);
    break;

  case BINARY:
    loadFromBinary(filename, elementType_);
    break;

  default:
    printf("Error in VolumetricMesh::VolumetricMesh: file format is unknown.\n");
    break;
  }
}

// parses the mesh, and returns the string corresponding to the element type
VolumetricMesh::VolumetricMesh(void *binaryInputStream, int numElementVertices_, elementType *elementType_, int memoryLoad):
  numElementVertices(numElementVertices_)
{
  if (memoryLoad)
    loadFromMemory((unsigned char *)binaryInputStream, elementType_);
  else
    loadFromBinary((FILE *)binaryInputStream, elementType_);
}

VolumetricMesh::~VolumetricMesh()
{
  for (int i = 0; i < numMaterials; i++)
    delete (materials[i]);
}

void VolumetricMesh::assignMaterialsToElements(int verbose)
{
  elementMaterial.assign(numElements, numMaterials);

  propagateRegionsToElements();

  // seek for unassigned elements
  std::set<int> unassignedElements;
  for (int el = 0; el < numElements; el++) {
    if (elementMaterial[el] == numMaterials)
      unassignedElements.insert(el);
  }

  if (unassignedElements.size() > 0) {
    // assign set and region to the unassigned elements

    // create a material if none exists
    if (numMaterials == 0) {
      numMaterials++;
      materials.assign(1, nullptr);
      materials[0] = new ENuMaterial("defaultMaterial", density_default, E_default, nu_default);
    }

    numSets++;
    sets.resize(numSets, Set("unassignedSet"));
    for (std::set<int>::iterator iter = unassignedElements.begin(); iter != unassignedElements.end(); iter++)
      sets[numSets - 1].insert(*iter);  // elements in sets are 0-indexed

    // create a new region for the unassigned elements
    numRegions++;
    regions.resize(numRegions, Region(numMaterials - 1, numSets - 1));

    for (std::set<int>::iterator iter = unassignedElements.begin(); iter != unassignedElements.end(); iter++)
      elementMaterial[*iter] = numMaterials - 1;

    if (verbose)
      printf("Warning: %d elements were not found in any of the regions. Using default material parameters for these elements.\n", (int)unassignedElements.size());
  }
}

void VolumetricMesh::loadFromAscii(const char *filename, elementType *elementType_, int verbose)
{
  // parse the .veg file
  VolumetricMeshParser volumetricMeshParser;

  if (volumetricMeshParser.open(filename) != 0) {
    printf("Error: could not open file %s.\n", filename);
    throw 1;
  }

  // === First pass: parse vertices and elements, and count the number of materials, sets and regions  ===

  int countNumVertices = 0;
  int countNumElements = 0;

  numElements = -1;
  numMaterials = 0;
  numSets = 1;  // set 0 is "allElements"
  numRegions = 0;
  if (elementType_)
    *elementType_ = INVALID;
  int parseState = 0;
  char lineBuffer[1024];

  int oneIndexedVertices = 1;
  int oneIndexedElements = 1;
  while (volumetricMeshParser.getNextLine(lineBuffer, 0, 0) != nullptr) {
    // lineBuffer now contains the next line
    // printf("%s\n", lineBuffer);

    // find *VERTICES
    if ((parseState == 0) && (strncmp(lineBuffer, "*VERTICES", 9) == 0)) {
      parseState = 1;

      if (volumetricMeshParser.getNextLine(lineBuffer, 0, 0) != nullptr) {
        // format is numVertices, 3, 0, 0
        sscanf(lineBuffer, "%d", &numVertices);  // ignore 3, 0, 0
        vertices.resize(numVertices);
      }
      else {
        printf("Error: file %s is not in the .veg format. Offending line:\n%s\n", filename, lineBuffer);
        throw 2;
      }

      continue;
    }

    // find *ELEMENTS
    if ((parseState == 1) && (strncmp(lineBuffer, "*ELEMENTS", 9) == 0)) {
      parseState = 2;

      // parse element type
      if (volumetricMeshParser.getNextLine(lineBuffer) != nullptr) {
        volumetricMeshParser.removeWhitespace(lineBuffer);

        if (strncmp(lineBuffer, "TET", 3) == 0) {
          if (elementType_)
            *elementType_ = TET;
        }
        else if (strncmp(lineBuffer, "CUBIC", 5) == 0) {
          if (elementType_)
            *elementType_ = CUBIC;
        }
        else {
          printf("Error: unknown mesh type %s in file %s\n", lineBuffer, filename);
          throw 3;
        }
      }
      else {
        printf("Error: file %s is not in the .veg format. Offending line:\n%s\n", filename, lineBuffer);
        throw 4;
      }

      // parse the number of elements
      if (volumetricMeshParser.getNextLine(lineBuffer, 0, 0) != nullptr) {
        // format is: numElements, numElementVertices, 0
        sscanf(lineBuffer, "%d", &numElements);  // only use numElements; ignore numElementVertices, 0
        elements.resize(numElements * numElementVertices);
      }
      else {
        printf("Error: file %s is not in the .veg format. Offending line:\n%s\n", filename, lineBuffer);
        throw 5;
      }

      continue;
    }

    if ((parseState == 2) && (lineBuffer[0] == '*')) {
      parseState = 3;  // end of elements
    }

    if (parseState == 1) {
      // read the vertex position
      if (countNumVertices >= numVertices) {
        printf("Error: mismatch in the number of vertices in %s.\n", filename);
        throw 6;
      }

      // ignore space, comma or tab
      char *ch = lineBuffer;
      while ((*ch == ' ') || (*ch == ',') || (*ch == '\t'))
        ch++;

      int index;
      sscanf(ch, "%d", &index);
      // seek next separator
      while ((*ch != ' ') && (*ch != ',') && (*ch != '\t') && (*ch != 0))
        ch++;

      if (index == 0)
        oneIndexedVertices = 0;  // input mesh has 0-indexed vertices

      double pos[3];
      for (int i = 0; i < 3; i++) {
        // ignore space, comma or tab
        while ((*ch == ' ') || (*ch == ',') || (*ch == '\t'))
          ch++;

        if (*ch == 0) {
          printf("Error parsing line %s in file %s.\n", lineBuffer, filename);
          throw 7;
        }

        sscanf(ch, "%lf", &pos[i]);

        // seek next separator
        while ((*ch != ' ') && (*ch != ',') && (*ch != '\t') && (*ch != 0))
          ch++;
      }

      vertices[countNumVertices] = Vec3d(pos);
      countNumVertices++;
    }

    if (parseState == 2) {
      // read the element vertices
      if (countNumElements >= numElements) {
        printf("Error: mismatch in the number of elements in %s.\n", filename);
        throw 8;
      }

      // ignore space, comma or tab
      char *ch = lineBuffer;
      while ((*ch == ' ') || (*ch == ',') || (*ch == '\t'))
        ch++;

      int index;
      sscanf(ch, "%d", &index);

      if (index == 0)
        oneIndexedElements = 0;  // input mesh has 0-indexed elements

      // seek next separator
      while ((*ch != ' ') && (*ch != ',') && (*ch != '\t') && (*ch != 0))
        ch++;

      for (int i = 0; i < numElementVertices; i++) {
        // ignore space, comma or tab
        while ((*ch == ' ') || (*ch == ',') || (*ch == '\t'))
          ch++;

        if (*ch == 0) {
          printf("Error parsing line %s in file %s.\n", lineBuffer, filename);
          throw 9;
        }

        int index = -1;
        sscanf(ch, "%d", &index);
        // if vertices were 1-numbered in the .veg file, convert to 0-numbered
        elements[countNumElements * numElementVertices + i] = index - oneIndexedVertices;

        // seek next separator
        while ((*ch != ' ') && (*ch != ',') && (*ch != '\t') && (*ch != 0))
          ch++;
      }

      countNumElements++;
    }

    if (strncmp(lineBuffer, "*MATERIAL", 9) == 0) {
      numMaterials++;
      continue;
    }

    if (strncmp(lineBuffer, "*SET", 4) == 0) {
      numSets++;
      continue;
    }

    if (strncmp(lineBuffer, "*REGION", 7) == 0) {
      numRegions++;
      continue;
    }
  }

  if (numElements < 0) {
    printf("Error: incorrect number of elements.  File %s may not be in the .veg format.\n", filename);
    throw 10;
  }

  // === Second pass: parse materials, sets and regions ===

  volumetricMeshParser.rewindToStart();

  if (verbose) {
    if (numMaterials == 0)
      printf("Warning: no materials encountered in %s.\n", filename);

    if (numRegions == 0)
      printf("Warning: no regions encountered in %s.\n", filename);
  }

  materials.assign(numMaterials, nullptr);
  sets.clear();
  regions.clear();
  sets.resize(numSets);
  regions.resize(numRegions);

  // create the "allElements" set, containing all the elements
  sets[0] = generateAllElementsSet(numElements);

  int countNumMaterials = 0;
  int countNumSets = 1;  // set 0 is "allElements"
  int countNumRegions = 0;

  std::map<std::string, int> materialMap;  // establishes a relationship between material's string name and its index in the "materials" array
  std::map<std::string, int> setMap;       // establishes a relationship between set's string name and its index in the "sets" array
  setMap.insert(std::pair<std::string, int>(sets[0].getName(), 0));

  parseState = 0;

  while (volumetricMeshParser.getNextLine(lineBuffer, 0, 0) != nullptr) {
    // printf("%s\n", lineBuffer);

    // exit parsing of comma-separated set elements upon the new * command
    if ((parseState == 11) && (lineBuffer[0] == '*'))
      parseState = 0;

    // parse material

    if ((parseState == 0) && (strncmp(lineBuffer, "*MATERIAL", 9) == 0)) {
      volumetricMeshParser.removeWhitespace(lineBuffer);

      // read material name
      char materialNameC[4096];
      strcpy(materialNameC, &lineBuffer[9]);

      // read the material type
      char materialSpecification[4096];
      if (volumetricMeshParser.getNextLine(lineBuffer) != nullptr) {
        volumetricMeshParser.removeWhitespace(lineBuffer);
        sscanf(lineBuffer, "%s", materialSpecification);
      }
      else {
        printf("Error: incorrect material in file %s. Offending line:\n%s\n", filename, lineBuffer);
        throw 11;
      }

      // seek for first comma in the material type specification
      char *ch = materialSpecification;
      while ((*ch != ',') && (*ch != 0))
        ch++;

      if (*ch == 0) {
        printf("Error parsing file %s. Offending line: %s.\n", filename, lineBuffer);
        throw 12;
      }

      // ch now points to the first comma
      // set the materialType (the string up to the first comma)
      char materialType[4096];
      unsigned int materialTypeLength = ch - materialSpecification;
      memcpy(materialType, materialSpecification, sizeof(unsigned char) * materialTypeLength);
      *(materialType + materialTypeLength) = 0;
      // materialType is now set

      ch++;
      // now, ch points to the first character after the comma

      if (strcmp(materialType, "ENU") == 0) {
        // material specified by E, nu, density
        double density, E, nu;
        sscanf(ch, "%lf,%lf,%lf", &density, &E, &nu);

        if ((E > 0) && (nu > -1.0) && (nu < 0.5) && (density > 0)) {
          // create new material
          std::string name(materialNameC);
          materials[countNumMaterials] = new ENuMaterial(name, density, E, nu);
          materialMap.insert(std::pair<std::string, int>(name, countNumMaterials));
        }
        else {
          printf("Error: incorrect material specification in file %s. Offending line: %s\n", filename, lineBuffer);
          throw 13;
        }
      }
      else if (strncmp(materialType, "ORTHOTROPIC", 11) == 0) {
        // orthotropic material
        double density = 0.0, E1 = 0.0, E2 = 0.0, E3 = 0.0, nu12 = 0.0, nu23 = 0.0, nu31 = 0.0, G12 = 0.0, G23 = 0.0, G31 = 0.0;
        double nu = 0.0, G = 1.0;
        bool useNuAndG = false;
        double R[9];               // rotation matrix, stored in row-major format
        memset(R, 0, sizeof(R));
        R[0] = R[4] = R[8] = 1.0;  // default to identity

        // find subtype
        char *subType = materialType + 11;
        bool enoughParameters = false;

        if ((sscanf(ch, "%lf,%lf,%lf,%lf", &density, &E1, &E2, &E3) == 4) && ((E1 > 0) && (E2 > 0) && (E3 > 0) && (density > 0))) {
          // move ch to the next parameters (skip over the four ones that were just read)
          for (int i = 0; i < 4; i++) {
            while ((*ch != ',') && (*ch != 0))
              ch++;
            if (*ch == 0)
              break;
            ch++;
          }
          // ch now points to the first character of the fifth parameter, or \0 if there is no fifth parameter
          // subtype points to the fitst character of the subtype suffix, or \0 if there is no suffix

          if ((*subType == 0) || (strcmp(subType, "_N3G3R9") == 0)) {
            // there is no suffix (*ORTHOTROPIC), or suffix is *ORTHOTROPIC_N3G3R9
            // parse nu12,nu23,nu31,G12,G23,G31,R0,R1,R2,R3,R4,R5,R6,R7,R8
            if (sscanf(ch, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                  &nu12, &nu23, &nu31, &G12, &G23, &G31,
                  &R[0], &R[1], &R[2], &R[3], &R[4], &R[5], &R[6], &R[7], &R[8]) == 15)
              enoughParameters = true;
          }
          else if (strcmp(subType, "_N3G3") == 0) {
            // *ORTHOTROPIC_N3G3
            // parse nu12,nu23,nu31,G12,G23,G31
            // R is default (identity)
            if (sscanf(ch, "%lf,%lf,%lf,%lf,%lf,%lf", &nu12, &nu23, &nu31, &G12, &G23, &G31) == 6)
              enoughParameters = true;
          }
          else if (strcmp(subType, "_N1G1R9") == 0) {
            // *ORTHOTROPIC_N1G1R9
            // parse nu,G,R0,R1,R2,R3,R4,R5,R6,R7,R8
            if (sscanf(ch, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                  &nu, &G, &R[0], &R[1], &R[2], &R[3], &R[4], &R[5], &R[6], &R[7], &R[8]) == 11) {
              useNuAndG = true;
              enoughParameters = true;
            }
          }
          else if (strcmp(subType, "_N1G1") == 0) {
            // *ORTHOTROPIC_N1G1
            // parse nu,G
            // R is default (identity)
            if (sscanf(ch, "%lf,%lf", &nu, &G) == 2) {
              useNuAndG = true;
              enoughParameters = true;
            }
          }
          else if (strcmp(subType, "_N1R9") == 0) {
            // *ORTHOTROPIC_N1R9
            // parse nu,R0,R1,R2,R3,R4,R5,R6,R7,R8
            // G is default (1.0)
            if (sscanf(ch, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                  &nu, &R[0], &R[1], &R[2], &R[3], &R[4], &R[5], &R[6], &R[7], &R[8]) == 10) {
              useNuAndG = true;
              enoughParameters = true;
            }
          }
          else if (strcmp(subType, "_N1") == 0) {
            // *ORTHOTROPIC_N1
            // parse nu
            // G is default (1.0)
            // R is default (identity)
            if (sscanf(ch, "%lf", &nu) == 1) {
              useNuAndG = true;
              enoughParameters = true;
            }
          }
          else {
            printf("Error: incorrect orthortropic material type \"%s\" in file %s. Offending line: %s\n", subType, filename, lineBuffer);
            throw 14;
          }

          if (useNuAndG && (nu > -1.0) && (nu < 0.5))  // if nu is not in this range, then G_ij is still 0, finally produce an error
          {
            // formulas from:
            // Yijing Li, Jernej Barbic: Stable Corotational Materials, Symposium on Computer Animation 2014
            nu12 = nu * sqrt(E1 / E2);
            nu23 = nu * sqrt(E2 / E3);
            nu31 = nu * sqrt(E3 / E1);
            G12 = G * sqrt(E1 * E2) / (2.0 * (1.0 + nu));
            G23 = G * sqrt(E2 * E3) / (2.0 * (1.0 + nu));
            G31 = G * sqrt(E3 * E1) / (2.0 * (1.0 + nu));
          }
        }  // end if (sscanf(ch,"%lf,%lf,%lf,%lf", &density, &E1, &E2, &E3) == 4 && (E1 > 0 && E2 > 0 && E3 > 0) )
        // if the if clause failed, enoughParameters is now false

        if (enoughParameters && (G12 > 0) && (G23 > 0) && (G31 > 0)) {
          // create new material
          std::string name(materialNameC);
          materials[countNumMaterials] = new OrthotropicMaterial(std::string(materialNameC), density, E1, E2, E3, nu12, nu23, nu31, G12, G23, G31, R);
          materialMap.insert(std::pair<std::string, int>(name, countNumMaterials));
        }
        else {
          printf("Error: incorrect material specification in file %s. Offending line: %s\n", filename, lineBuffer);
          throw 14;
        }
      }
      else if (strncmp(materialType, "MOONEYRIVLIN", 12) == 0) {
        // mu01, m10, v1, density
        double density, mu01, mu10, v1;
        sscanf(ch, "%lf,%lf,%lf,%lf", &density, &mu01, &mu10, &v1);

        if (density > 0) {
          // create new material
          std::string name(materialNameC);
          materials[countNumMaterials] = new MooneyRivlinMaterial(name, density, mu01, mu10, v1);
          materialMap.insert(std::pair<std::string, int>(name, countNumMaterials));
        }
        else {
          printf("Error: incorrect material specification in file %s. Offending line:\n%s\n", filename, lineBuffer);
          throw 15;
        }
      }
      else {
        printf("Error: incorrect material specification in file %s. Offending line:\n%s\n", filename, lineBuffer);
        throw 16;
      }

      countNumMaterials++;
    }

    // parse region
    if ((parseState == 0) && (strncmp(lineBuffer, "*REGION,", 7) == 0)) {
      volumetricMeshParser.removeWhitespace(lineBuffer);

      char setNameC[4096];
      char materialNameC[4096];

      if (volumetricMeshParser.getNextLine(lineBuffer) != nullptr) {
        volumetricMeshParser.removeWhitespace(lineBuffer);

        // format is set, material
        // seek for first comma
        char *ch = lineBuffer;
        while ((*ch != ',') && (*ch != 0))
          ch++;

        if (*ch == 0) {
          printf("Error parsing file %s. Offending line: %s.\n", filename, lineBuffer);
          throw 17;
        }

        // now, ch points to the comma
        *ch = 0;
        strcpy(setNameC, lineBuffer);
        *ch = ',';  // restore the lineBuffer
        ch++;
        strcpy(materialNameC, ch);
      }
      else {
        printf("Error: file %s is not in the .veg format. Offending line:\n%s\n", filename, lineBuffer);
        throw 18;
      }

      // seek for setNameC
      int setNum = -1;
      std::map<std::string, int>::iterator it = setMap.find(std::string(setNameC));
      if (it != setMap.end()) {
        setNum = it->second;
      }
      else {
        printf("Error: set \"%s\" not found among the sets.\n", setNameC);
        printf("All existing sets: \n");
        for (const auto &p : setMap)
          printf("%i: \"%s\"\n", p.second, p.first.c_str());
        printf("\n");
        throw 19;
      }

      // seek for materialNameC
      int materialNum = -1;
      it = materialMap.find(std::string(materialNameC));
      if (it != materialMap.end()) {
        materialNum = it->second;
      }
      else {
        printf("Error: material %s not found among the materials.\n", materialNameC);
        throw 20;
      }

      // create a new region
      regions[countNumRegions] = Region(materialNum, setNum);
      countNumRegions++;
    }

    // parse set elements
    if (parseState == 11) {
      // parse the next line of the comma-separated elements in the set
      // we know that lineBuffer[0] != '*' (i.e., not the end of the list), as that case was already previously handled

      volumetricMeshParser.removeWhitespace(lineBuffer);
      // printf("%s\n", lineBuffer);

      // parse the comma-separated line
      char *pch;
      pch = strtok(lineBuffer, ",");
      while ((pch != nullptr) && (isdigit(*pch))) {
        int newElement = atoi(pch);
        int ind = newElement - oneIndexedElements;
        if (ind >= numElements || ind < 0) {
          printf("Error: set element index: %d out of bounds.\n", newElement);
          throw 21;
        }
        sets[countNumSets - 1].insert(ind);  // sets are 0-indexed, but .veg files may be 1-indexed (oneIndexedElements == 1)
        pch = strtok(nullptr, ",");
      }
    }

    // parse set
    if ((parseState == 0) && (strncmp(lineBuffer, "*SET", 4) == 0)) {
      volumetricMeshParser.removeWhitespace(lineBuffer);

      std::string name(BasicAlgorithms::stripLight(&lineBuffer[4]));
      sets[countNumSets] = Set(name);
      setMap.insert(std::pair<std::string, int>(name, countNumSets));
      countNumSets++;
      parseState = 11;
    }
  }

  // === assign materials to elements and seek for unassigned elements ===
  assignMaterialsToElements(verbose);

  volumetricMeshParser.close();
}

VolumetricMesh::VolumetricMesh(int numVertices_, const double *vertices_,
  int numElements_, int numElementVertices_, const int *elements_,
  double E, double nu, double density):
  numElementVertices(numElementVertices_)
{
  numElements = numElements_;
  numVertices = numVertices_;

  numMaterials = 1;
  numSets = 1;
  numRegions = 1;

  vertices.resize(numVertices);
  elements.resize(numElementVertices * numElements);
  elementMaterial.assign(numElements, 0);
  materials.assign(numMaterials, nullptr);
  sets.resize(numSets);
  regions.assign(numRegions, Region(0, 0));

  memcpy(vertices.data(), vertices_, sizeof(double) * numVertices * 3);
  memcpy(elements.data(), elements_, sizeof(int) * numElementVertices * numElements);

  Material *material = new ENuMaterial("defaultMaterial", density, E, nu);
  materials[0] = material;

  sets[0] = generateAllElementsSet(numElements);
  // we don't need to call propagateRegionsToElements here because elementMaterial has been set
}

VolumetricMesh::VolumetricMesh(int numVertices_, const double *vertices_,
  int numElements_, int numElementVertices_, const int *elements_,
  int numMaterials_, const Material *const *materials_,
  int numSets_, const Set *sets_,
  int numRegions_, const Region *regions_):
  numElementVertices(numElementVertices_)
{
  numElements = numElements_;
  numVertices = numVertices_;

  numMaterials = numMaterials_;
  numSets = numSets_;
  numRegions = numRegions_;

  vertices.resize(numVertices);
  elements.resize(numElementVertices * numElements);
  elementMaterial.assign(numElements, 0);
  materials.assign(numMaterials, nullptr);
  sets.resize(numSets);
  regions.resize(numRegions);

  memcpy(vertices.data(), vertices_, sizeof(double) * numVertices * 3);
  memcpy(elements.data(), elements_, sizeof(int) * numElementVertices * numElements);

  for (int i = 0; i < numMaterials; i++)
    materials[i] = materials_[i]->clone();

  for (int i = 0; i < numSets; i++)
    sets[i] = sets_[i];

  for (int i = 0; i < numRegions; i++)
    regions[i] = regions_[i];

  // set elementMaterial:
  propagateRegionsToElements();
}

void VolumetricMesh::loadFromMemory(unsigned char *binaryInputStream, elementType *elementType_)
{
  int memoryLoad = 1;
  loadFromBinaryGeneric((void *)binaryInputStream, elementType_, memoryLoad);
}

void VolumetricMesh::loadFromBinary(FILE *binaryInputStream, elementType *elementType_)
{
  int memoryLoad = 0;
  loadFromBinaryGeneric((void *)binaryInputStream, elementType_, memoryLoad);
}

void VolumetricMesh::loadFromBinaryGeneric(void *binaryInputStream_, elementType *elementType_, int memoryLoad)
{
  unsigned int (*genericRead)(void *, unsigned int, unsigned int, void *);

  void *binaryInputStream;
  if (memoryLoad) {
    genericRead = &VolumetricMesh::readFromMemory;
    binaryInputStream = &(binaryInputStream_);
  }
  else {
    genericRead = &VolumetricMesh::readFromFile;
    binaryInputStream = binaryInputStream_;
  }

  double version;
  if ((int)genericRead(&version, sizeof(double), 1, binaryInputStream) != 1) {
    printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read version.\n");
    throw 0;
  };

  int eleType;
  if ((int)genericRead(&eleType, sizeof(int), 1, binaryInputStream) != 1) {
    printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read element type.\n");
    throw 0;
  };

  switch (eleType) {
  case TET:
    if (elementType_)
      *elementType_ = TET;
    break;
  case CUBIC:
    if (elementType_)
      *elementType_ = CUBIC;
    break;
  default:
    printf("Error in VolumetricMesh::loadFromBinaryGeneric: unknown mesh type %d in file stream\n", eleType);
    throw 2;
    break;
  }

  // input the number of vertices
  if ((int)genericRead(&numVertices, sizeof(int), 1, binaryInputStream) != 1) {
    printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read numVertices.\n");
    throw 0;
  }

  if (numVertices < 0) {
    printf("Error in VolumetricMesh::loadFromBinaryGeneric: incorrect number of vertices.\n");
    throw 3;
  }

  vertices.resize(numVertices);

  // input all the vertices
  if ((int)genericRead(vertices.data(), sizeof(double), 3 * numVertices, binaryInputStream) != 3 * numVertices) {
    printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read vertex coordinates.\n");
    throw 0;
  }

  // input the number of elements
  if ((int)genericRead(&numElements, sizeof(int), 1, binaryInputStream) != 1) {
    printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read numElements.\n");
    throw 0;
  }

  if (numElements < 0) {
    printf("Error in VolumetricMesh::loadFromBinaryGeneric: incorrect number of elements.\n");
    throw 4;
  }
  // input the number of vertices of every element
  int numEleVer;
  if ((int)genericRead(&numEleVer, sizeof(int), 1, binaryInputStream) != 1) {
    printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read number of vertices in every element.\n");
    throw 0;
  }

  if (numElementVertices > 0 && numEleVer != numElementVertices) {
    printf("Error in VolumetricMesh::loadFromBinaryGeneric: mismatch in the number of vertices of every element from file stream.\n");
    throw 5;
  }

  // input all elements
  elements.resize(numElements * numElementVertices);
  if ((int)genericRead(elements.data(), sizeof(int), numElementVertices * numElements, binaryInputStream) != numElementVertices * numElements) {
    printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the elements.\n");
    throw 0;
  }

  // input number of materials
  if ((int)genericRead(&numMaterials, sizeof(int), 1, binaryInputStream) != 1) {
    printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the numMaterials.\n");
    throw 0;
  }

  if (numMaterials < 0) {
    printf("Error in VolumetricMesh::loadFromBinaryGeneric: incorrect number of materials.\n");
    throw 6;
  }
  // input all the materials
  materials.assign(numMaterials, nullptr);
  for (int materialIndex = 0; materialIndex < numMaterials; materialIndex++) {
    // input material name
    char materialName[4096];
    int length;
    if ((int)genericRead(&length, sizeof(int), 1, binaryInputStream) != 1) {
      printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the length of material name.\n");
      throw 0;
    }
    if ((int)genericRead(materialName, sizeof(char), length, binaryInputStream) != length) {
      printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the material name.\n");
      throw 0;
    }
    materialName[length] = '\0';

    // input material type
    int matType;
    if ((int)genericRead(&matType, sizeof(int), 1, binaryInputStream) != 1) {
      printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the length of material type.\n");
      throw 0;
    }

    switch (matType) {
    case Material::ENU: {
      double materialProperty[Material::ENU_NUM_PROPERTIES];
      if ((int)genericRead(materialProperty, sizeof(double), Material::ENU_NUM_PROPERTIES, binaryInputStream) != Material::ENU_NUM_PROPERTIES) {
        printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the material properties (ENU).\n");
        throw 0;
      }

      if ((materialProperty[Material::ENU_E] > 0) && (materialProperty[Material::ENU_NU] > -1.0) && (materialProperty[Material::ENU_NU] < 0.5) && (materialProperty[Material::ENU_DENSITY] > 0)) {
        // create new material
        materials[materialIndex] = new ENuMaterial(materialName, materialProperty[Material::ENU_DENSITY], materialProperty[Material::ENU_E], materialProperty[Material::ENU_NU]);
      }
      else {
        printf("Error in VolumetricMesh::loadFromBinaryGeneric: incorrect material specification in file stream.\n");
        throw 7;
      }
    } break;

    case Material::ORTHOTROPIC: {
      double materialProperty[Material::ORTHOTROPIC_NUM_PROPERTIES];
      if ((int)genericRead(materialProperty, sizeof(double), Material::ORTHOTROPIC_NUM_PROPERTIES, binaryInputStream) != Material::ORTHOTROPIC_NUM_PROPERTIES) {
        printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the material properties (ORTHOTROPIC).\n");
        throw 0;
      }
      double R[9];
      if ((int)genericRead(R, sizeof(double), 9, binaryInputStream) != 9) {
        printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the material properties (ORTHOTROPIC).\n");
        throw 0;
      }
      if ((materialProperty[Material::ORTHOTROPIC_E1] > 0) && (materialProperty[Material::ORTHOTROPIC_E2] > 0) && (materialProperty[Material::ORTHOTROPIC_E3] > 0) && (materialProperty[Material::ORTHOTROPIC_G12] > 0) && (materialProperty[Material::ORTHOTROPIC_G23] > 0) && (materialProperty[Material::ORTHOTROPIC_G31] > 0) && (materialProperty[Material::ORTHOTROPIC_DENSITY] > 0)) {
        // create new material
        materials[materialIndex] = new OrthotropicMaterial(materialName, materialProperty[Material::ORTHOTROPIC_DENSITY], materialProperty[Material::ORTHOTROPIC_E1], materialProperty[Material::ORTHOTROPIC_E2], materialProperty[Material::ORTHOTROPIC_E3], materialProperty[Material::ORTHOTROPIC_NU12], materialProperty[Material::ORTHOTROPIC_NU23], materialProperty[Material::ORTHOTROPIC_NU31], materialProperty[Material::ORTHOTROPIC_G12], materialProperty[Material::ORTHOTROPIC_G23], materialProperty[Material::ORTHOTROPIC_G31], R);
      }
      else {
        printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the material properties (ORTHOTROPIC).\n");
        throw 14;
      }
    } break;

    case Material::MOONEYRIVLIN: {
      double materialProperty[Material::MOONEYRIVLIN_NUM_PROPERTIES];
      if ((int)genericRead(materialProperty, sizeof(double), Material::MOONEYRIVLIN_NUM_PROPERTIES, binaryInputStream) != Material::MOONEYRIVLIN_NUM_PROPERTIES) {
        printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the material properties (MOONEYRIVLIN).\n");
        throw 0;
      }
      if (materialProperty[Material::MOONEYRIVLIN_DENSITY] > 0) {
        // create new material
        materials[materialIndex] = new MooneyRivlinMaterial(materialName, materialProperty[Material::MOONEYRIVLIN_DENSITY], materialProperty[Material::MOONEYRIVLIN_MU01], materialProperty[Material::MOONEYRIVLIN_MU10], materialProperty[Material::MOONEYRIVLIN_V1]);
      }
      else {
        printf("Error in VolumetricMesh::loadFromBinaryGeneric: incorrect material specification in file stream. \n");
        throw 8;
      }
    } break;

    default:
      printf("Error in VolumetricMesh::loadFromBinaryGeneric: material type %d is unknown.\n", matType);
      throw 9;
      break;
    }
  }  // for materialIndex

  // input the number of sets
  if ((int)genericRead(&numSets, sizeof(int), 1, binaryInputStream) != 1) {
    printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the numSets.\n");
    throw 0;
  }

  sets.clear();
  sets.resize(numSets);

  // create the "allElements" set, containing all the elements
  sets[0] = generateAllElementsSet(numElements);

  std::vector<int> intTempVec;
  // input all the other sets (must start from set 1)
  for (int setIndex = 1; setIndex < numSets; setIndex++) {
    // input the name of the set
    char setName[4096];
    int length;
    if ((int)genericRead(&length, sizeof(int), 1, binaryInputStream) != 1) {
      printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the length of the set name.\n");
      throw 0;
    }

    if ((int)genericRead(setName, sizeof(char), length, binaryInputStream) != length) {
      printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the set name.\n");
      throw 0;
    }
    setName[length] = '\0';
    sets[setIndex] = Set(setName);

    // input number of elements in the current set
    int cardinality;
    if ((int)genericRead(&cardinality, sizeof(int), 1, binaryInputStream) != 1) {
      printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the number of elements in current set.\n");
      throw 0;
    }

    intTempVec.resize(cardinality);
    // input all the elements in the current set
    if ((int)genericRead(intTempVec.data(), sizeof(int), cardinality, binaryInputStream) != cardinality) {
      printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the elements in current set.\n");
      throw 0;
    }

    for (int setElementIndex = 0; setElementIndex < cardinality; setElementIndex++)
      sets[setIndex].insert(intTempVec[setElementIndex]);
  }

  // input the number of regions
  if ((int)genericRead(&numRegions, sizeof(int), 1, binaryInputStream) != 1) {
    printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the numRegions.\n");
    throw 0;
  }

  // input all the regions, for each region (material, set)
  regions.clear();
  regions.resize(numRegions);
  for (int regionIndex = 0; regionIndex < numRegions; regionIndex++) {
    int materialIndex;
    if ((int)genericRead(&materialIndex, sizeof(int), 1, binaryInputStream) != 1) {
      printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the materialIndex.\n");
      throw 0;
    }
    int setIndex;
    if ((int)genericRead(&setIndex, sizeof(int), 1, binaryInputStream) != 1) {
      printf("Error in VolumetricMesh::loadFromBinaryGeneric: cannot read the setIndex.\n");
      throw 0;
    }
    regions[regionIndex] = Region(materialIndex, setIndex);
  }  // for regionIndex

  // === assign materials to elements and handle the unassigned elements ===
  int verbose = 0;
  assignMaterialsToElements(verbose);
}

VolumetricMesh::VolumetricMesh(const VolumetricMesh &volumetricMesh)
{
  numVertices = volumetricMesh.numVertices;
  vertices = volumetricMesh.vertices;
  numElementVertices = volumetricMesh.numElementVertices;
  numElements = volumetricMesh.numElements;
  elements = volumetricMesh.elements;
  numMaterials = volumetricMesh.numMaterials;
  numSets = volumetricMesh.numSets;
  numRegions = volumetricMesh.numRegions;

  materials.assign(numMaterials, nullptr);
  for (int i = 0; i < numMaterials; i++)
    materials[i] = (volumetricMesh.materials)[i]->clone();

  sets = volumetricMesh.sets;
  regions = volumetricMesh.regions;

  elementMaterial = volumetricMesh.elementMaterial;
}

void VolumetricMesh::loadFromBinary(const char *filename, elementType *elementType_)
{
  FILE *fin = fopen(filename, "rb");

  if (fin == nullptr) {
    printf("Error in VolumetricMesh::loadFromBinary: could not open file %s.\n", filename);
    throw 1;
  }
  loadFromBinary(fin, elementType_);

  fclose(fin);
}

int VolumetricMesh::saveToAscii(const char *filename, elementType elementType_) const  // saves the mesh to a .veg file
{
  FILE *fout = fopen(filename, "w");
  if (!fout) {
    printf("Error: could not write to %s.\n", filename);
    return 1;
  }

  fprintf(fout, "# Vega mesh file.\n");
  fprintf(fout, "# %d vertices, %d elements\n", numVertices, numElements);
  fprintf(fout, "\n");

  // write vertices
  fprintf(fout, "*VERTICES\n");
  fprintf(fout, "%d 3 0 0\n", numVertices);

  for (int i = 0; i < numVertices; i++) {
    const Vec3d &v = getVertex(i);
    fprintf(fout, "%d %.15G %.15G %.15G\n", i + 1, v[0], v[1], v[2]);
  }
  fprintf(fout, "\n");

  // write elements
  fprintf(fout, "*ELEMENTS\n");

  char elementName[4096] = "INVALID";
  if (elementType_ == TET)
    strcpy(elementName, "TET");
  else if (elementType_ == CUBIC)
    strcpy(elementName, "CUBIC");
  fprintf(fout, "%s\n", elementName);

  fprintf(fout, "%d %d 0\n", numElements, numElementVertices);

  for (int el = 0; el < numElements; el++) {
    fprintf(fout, "%d ", el + 1);
    for (int j = 0; j < numElementVertices; j++) {
      fprintf(fout, "%d", getVertexIndex(el, j) + 1);
      if (j != numElementVertices - 1)
        fprintf(fout, " ");
    }
    fprintf(fout, "\n");
  }
  fprintf(fout, "\n");

  // write materials
  for (int materialIndex = 0; materialIndex < numMaterials; materialIndex++) {
    std::string name = materials[materialIndex]->getName();
    fprintf(fout, "*MATERIAL %s\n", name.c_str());

    if (materials[materialIndex]->getType() == Material::ENU) {
      ENuMaterial *material = downcastENuMaterial(materials[materialIndex]);
      double density = material->getDensity();
      double E = material->getE();
      double nu = material->getNu();
      fprintf(fout, "ENU, %.15G, %.15G, %.15G\n", density, E, nu);
    }

    if (materials[materialIndex]->getType() == Material::ORTHOTROPIC) {
      OrthotropicMaterial *material = downcastOrthotropicMaterial(materials[materialIndex]);
      double density = material->getDensity();

      double E1 = material->getE1();
      double E2 = material->getE2();
      double E3 = material->getE3();
      double nu12 = material->getNu12();
      double nu23 = material->getNu23();
      double nu31 = material->getNu31();
      double G12 = material->getG12();
      double G23 = material->getG23();
      double G31 = material->getG31();
      double R[9];
      material->getR(R);

      fprintf(fout, "ORTHOTROPIC, %.15G, %.15G, %.15G, %.15G, %.15G, %.15G, %.15G, %.15G, %.15G, %.15G, %.15G, %.15G, %.15G, %.15G, %.15G, %.15G, %.15G, %.15G, %.15G\n", density, E1, E2, E3, nu12, nu23, nu31, G12, G23, G31, R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8]);
    }

    if (materials[materialIndex]->getType() == Material::MOONEYRIVLIN) {
      MooneyRivlinMaterial *material = downcastMooneyRivlinMaterial(materials[materialIndex]);
      double density = material->getDensity();
      double mu01 = material->getmu01();
      double mu10 = material->getmu10();
      double v1 = material->getv1();
      fprintf(fout, "MOONEYRIVLIN, %.15G, %.15G, %.20G %.15G\n", density, mu01, mu10, v1);
    }
    fprintf(fout, "\n");
  }

  // write sets (skip the allElements set)
  for (int setIndex = 1; setIndex < numSets; setIndex++) {
    const std::string &name = sets[setIndex].getName();
    fprintf(fout, "*SET %s\n", name.c_str());
    const std::set<int> &setElements = sets[setIndex].getElements();
    int count = 0;
    for (std::set<int>::iterator iter = setElements.begin(); iter != setElements.end(); iter++) {
      fprintf(fout, "%d, ", *iter + 1);  // .veg files are 1-indexed
      count++;
      if (count == 8) {
        fprintf(fout, "\n");
        count = 0;
      }
    }
    if (count != 0)
      fprintf(fout, "\n");
    fprintf(fout, "\n");
  }

  // write regions
  for (int regionIndex = 0; regionIndex < numRegions; regionIndex++) {
    int materialIndex = regions[regionIndex].getMaterialIndex();
    int setIndex = regions[regionIndex].getSetIndex();

    fprintf(fout, "*REGION\n");
    fprintf(fout, "%s, %s\n", sets[setIndex].getName().c_str(), materials[materialIndex]->getName().c_str());
    fprintf(fout, "\n");
  }

  fclose(fout);
  return 0;
}

int VolumetricMesh::save(const char *filename) const  //  for backward compatibility
{
  fileFormatType fileType = getFileFormatTypeByExt(filename);
  if (fileType == ASCII)
    return saveToAscii(filename);
  if (fileType == BINARY)
    return saveToBinary(filename);
  return saveToAscii(filename);
}

int VolumetricMesh::saveToBinary(const char *filename, unsigned int *bytesWritten, elementType elementType_) const
{
  FILE *fout = fopen(filename, "wb");
  if (!fout) {
    printf("Error: could not write to %s.\n", filename);
    return 1;
  }
  int code = saveToBinary(fout, bytesWritten, elementType_);
  fclose(fout);
  return code;
}

unsigned int VolumetricMesh::readFromFile(void *buf, unsigned int elementSize, unsigned int numElements, void *fin)
{
  return fread(buf, elementSize, numElements, (FILE *)fin);
}

unsigned int VolumetricMesh::readFromMemory(void *buf, unsigned int elementSize, unsigned int numElements, void *memoryLocation_)
{
  unsigned char *memoryLocation = (unsigned char *)(*(void **)(memoryLocation_));
  unsigned int numBytes = elementSize * numElements;
  memcpy(buf, memoryLocation, numBytes);
  (*(void **)(memoryLocation_)) = (void *)((unsigned char *)(*(void **)(memoryLocation_)) + numBytes);
  return numElements;
}

int VolumetricMesh::saveToBinary(FILE *binaryOutputStream, unsigned int *bytesWritten, elementType elementType_, bool countBytesOnly) const
{
  unsigned int totalBytesWritten = 0;
  unsigned int itemsWritten;

  // output the binary file version (1x double)
  const double version = 1.0;
  itemsWritten = 1;
  if (!countBytesOnly)
    itemsWritten = fwrite(&version, sizeof(double), 1, binaryOutputStream);
  if (itemsWritten < 1)
    return 1;
  totalBytesWritten += itemsWritten * sizeof(double);

  // output the element type (1x int)
  int eleType = elementType_;
  itemsWritten = 1;
  if (!countBytesOnly)
    itemsWritten = fwrite(&eleType, sizeof(int), 1, binaryOutputStream);
  if (itemsWritten < 1)
    return 1;
  totalBytesWritten += itemsWritten * sizeof(int);

  // output the number of vertices (1x int)
  itemsWritten = 1;
  if (!countBytesOnly)
    itemsWritten = fwrite(&numVertices, sizeof(int), 1, binaryOutputStream);
  if (itemsWritten < 1)
    return 1;
  totalBytesWritten += itemsWritten * sizeof(int);

  // output the vertices (3 doubles x numVertices)
  for (int vertexIndex = 0; vertexIndex < numVertices; vertexIndex++) {
    const Vec3d &v = getVertex(vertexIndex);
    itemsWritten = 3;
    if (!countBytesOnly)
      itemsWritten = fwrite(v.data(), sizeof(double), 3, binaryOutputStream);
    if (itemsWritten < 3)
      return 1;
    totalBytesWritten += itemsWritten * sizeof(double);
  }

  // output the number of elements (1x int)
  itemsWritten = 1;
  if (!countBytesOnly)
    itemsWritten = fwrite(&numElements, sizeof(int), 1, binaryOutputStream);
  if (itemsWritten < 1)
    return 1;
  totalBytesWritten += itemsWritten * sizeof(int);

  // output the number of vertices of every element (1x int)
  itemsWritten = 1;
  if (!countBytesOnly)
    itemsWritten = fwrite(&numElementVertices, sizeof(int), 1, binaryOutputStream);
  if (itemsWritten < 1)
    return 1;
  totalBytesWritten += itemsWritten * sizeof(int);

  // output the vertex indices of every element
  for (int elementIndex = 0; elementIndex < numElements; elementIndex++) {
    for (int vertexIndex = 0; vertexIndex < numElementVertices; vertexIndex++) {
      int temp = getVertexIndex(elementIndex, vertexIndex);
      itemsWritten = 1;
      if (!countBytesOnly)
        itemsWritten = fwrite(&temp, sizeof(int), 1, binaryOutputStream);
      if (itemsWritten < 1)
        return 1;
      totalBytesWritten += itemsWritten * sizeof(int);
    }
  }

  // output the number of materials
  itemsWritten = 1;
  if (!countBytesOnly)
    itemsWritten = fwrite(&numMaterials, sizeof(int), 1, binaryOutputStream);
  if (itemsWritten < 1)
    return 1;
  totalBytesWritten += itemsWritten * sizeof(int);

  // output materials
  for (int materialIndex = 0; materialIndex < numMaterials; materialIndex++) {
    std::string name = materials[materialIndex]->getName();
    unsigned int length = name.size();

    // output material name
    itemsWritten = 1;
    if (!countBytesOnly)
      itemsWritten = fwrite(&length, sizeof(int), 1, binaryOutputStream);  // write the length of the material name
    if (itemsWritten < 1)
      return 1;
    totalBytesWritten += itemsWritten * sizeof(int);

    itemsWritten = length;
    if (!countBytesOnly)
      itemsWritten = fwrite(name.c_str(), sizeof(char), length, binaryOutputStream);
    if (itemsWritten < length)
      return 1;
    totalBytesWritten += itemsWritten * sizeof(char);

    // output material type (1x int)
    int matType = materials[materialIndex]->getType();

    itemsWritten = 1;
    if (!countBytesOnly)
      itemsWritten = fwrite(&matType, sizeof(int), 1, binaryOutputStream);
    if (itemsWritten < 1)
      return 1;
    totalBytesWritten += itemsWritten * sizeof(int);
    switch (matType) {
    case Material::ENU: {
      ENuMaterial *material = downcastENuMaterial(materials[materialIndex]);
      double materialProperty[Material::ENU_NUM_PROPERTIES];
      materialProperty[Material::ENU_DENSITY] = material->getDensity();
      materialProperty[Material::ENU_E] = material->getE();
      materialProperty[Material::ENU_NU] = material->getNu();

      itemsWritten = Material::ENU_NUM_PROPERTIES;
      if (!countBytesOnly)
        itemsWritten = fwrite(materialProperty, sizeof(double), Material::ENU_NUM_PROPERTIES, binaryOutputStream);
      if (itemsWritten < Material::ENU_NUM_PROPERTIES)
        return 1;
      totalBytesWritten += itemsWritten * sizeof(double);
    } break;

    case Material::MOONEYRIVLIN: {
      MooneyRivlinMaterial *material = downcastMooneyRivlinMaterial(materials[materialIndex]);
      double materialProperty[Material::MOONEYRIVLIN_NUM_PROPERTIES];
      materialProperty[Material::MOONEYRIVLIN_DENSITY] = material->getDensity();
      materialProperty[Material::MOONEYRIVLIN_MU01] = material->getmu01();
      materialProperty[Material::MOONEYRIVLIN_MU10] = material->getmu10();
      materialProperty[Material::MOONEYRIVLIN_V1] = material->getv1();

      itemsWritten = Material::MOONEYRIVLIN_NUM_PROPERTIES;
      if (!countBytesOnly)
        itemsWritten = fwrite(materialProperty, sizeof(double), Material::MOONEYRIVLIN_NUM_PROPERTIES, binaryOutputStream);
      if (itemsWritten < Material::MOONEYRIVLIN_NUM_PROPERTIES)
        return 1;
      totalBytesWritten += itemsWritten * sizeof(double);
    } break;

    case Material::ORTHOTROPIC: {
      OrthotropicMaterial *material = downcastOrthotropicMaterial(materials[materialIndex]);
      double materialProperty[Material::ORTHOTROPIC_NUM_PROPERTIES];
      materialProperty[Material::ORTHOTROPIC_DENSITY] = material->getDensity();
      materialProperty[Material::ORTHOTROPIC_E1] = material->getE1();
      materialProperty[Material::ORTHOTROPIC_E2] = material->getE2();
      materialProperty[Material::ORTHOTROPIC_E3] = material->getE3();
      materialProperty[Material::ORTHOTROPIC_NU12] = material->getNu12();
      materialProperty[Material::ORTHOTROPIC_NU23] = material->getNu23();
      materialProperty[Material::ORTHOTROPIC_NU31] = material->getNu31();
      materialProperty[Material::ORTHOTROPIC_G12] = material->getG12();
      materialProperty[Material::ORTHOTROPIC_G23] = material->getG23();
      materialProperty[Material::ORTHOTROPIC_G31] = material->getG31();

      itemsWritten = Material::ORTHOTROPIC_NUM_PROPERTIES;
      if (!countBytesOnly)
        itemsWritten = fwrite(materialProperty, sizeof(double), Material::ORTHOTROPIC_NUM_PROPERTIES, binaryOutputStream);
      if (itemsWritten < Material::ORTHOTROPIC_NUM_PROPERTIES)
        return 1;
      totalBytesWritten += itemsWritten * sizeof(double);

      double R[9];
      material->getR(R);

      itemsWritten = 9;
      if (!countBytesOnly)
        itemsWritten = fwrite(R, sizeof(double), 9, binaryOutputStream);
      if (itemsWritten < 9)
        return 1;
      totalBytesWritten += itemsWritten * sizeof(double);
    } break;

    default: {
      printf("Error: material type %d is unknown.\n", matType);
      return 1;
    } break;
    }

  }  // for materialIndex

  // output the number of sets
  itemsWritten = 1;
  if (!countBytesOnly)
    itemsWritten = fwrite(&numSets, sizeof(int), 1, binaryOutputStream);
  if (itemsWritten < 1)
    return 1;
  totalBytesWritten += itemsWritten * sizeof(int);

  // output sets (skip the allElements set whose index is 0)
  for (int setIndex = 1; setIndex < numSets; setIndex++) {
    const auto &name = sets[setIndex].getName();

    // output the name of the set
    unsigned int length = name.size();
    itemsWritten = 1;
    if (!countBytesOnly)
      itemsWritten = fwrite(&length, sizeof(int), 1, binaryOutputStream);
    if (itemsWritten < 1)
      return 1;
    totalBytesWritten += itemsWritten * sizeof(int);

    itemsWritten = length;
    if (!countBytesOnly)
      itemsWritten = fwrite(name.c_str(), sizeof(char), length, binaryOutputStream);
    if (itemsWritten < length)
      return 1;
    totalBytesWritten += itemsWritten * sizeof(char);

    std::set<int> setElements;
    sets[setIndex].getElements(setElements);

    // output the number of elements in the current set
    unsigned int cardinality = setElements.size();
    itemsWritten = 1;
    if (!countBytesOnly)
      itemsWritten = fwrite(&cardinality, sizeof(int), 1, binaryOutputStream);
    if (itemsWritten < 1)
      return 1;
    totalBytesWritten += itemsWritten * sizeof(int);

    int *elementGroup = (int *)malloc(sizeof(int) * cardinality);
    int count = 0;
    for (std::set<int>::iterator iter = setElements.begin(); iter != setElements.end(); iter++) {
      elementGroup[count] = *iter;
      count++;
    }
    // output the current set
    itemsWritten = cardinality;
    if (!countBytesOnly)
      itemsWritten = fwrite(elementGroup, sizeof(int), cardinality, binaryOutputStream);
    if (itemsWritten < cardinality)
      return 1;
    totalBytesWritten += itemsWritten * sizeof(int);
    free(elementGroup);
  }  // for setIndex

  // output the number of regions
  itemsWritten = 1;
  if (!countBytesOnly)
    itemsWritten = fwrite(&numRegions, sizeof(int), 1, binaryOutputStream);
  if (itemsWritten < 1)
    return 1;
  totalBytesWritten += itemsWritten * sizeof(int);

  // output regions (for each region, output material and set)
  for (int regionIndex = 0; regionIndex < numRegions; regionIndex++) {
    int materialIndex = regions[regionIndex].getMaterialIndex();
    itemsWritten = 1;
    if (!countBytesOnly)
      itemsWritten = fwrite(&materialIndex, sizeof(int), 1, binaryOutputStream);
    if (itemsWritten < 1)
      return 1;
    totalBytesWritten += itemsWritten * sizeof(int);

    int setIndex = regions[regionIndex].getSetIndex();
    itemsWritten = 1;
    if (!countBytesOnly)
      itemsWritten = fwrite(&setIndex, sizeof(int), 1, binaryOutputStream);
    if (itemsWritten < 1)
      return 1;
    totalBytesWritten += itemsWritten * sizeof(int);
  }

  if (bytesWritten != nullptr)
    *bytesWritten = totalBytesWritten;

  return 0;
}

VolumetricMesh::elementType VolumetricMesh::getElementTypeASCII(const char *filename)
{
  // printf("Parsing %s... (for element type determination)\n",filename);fflush(nullptr);
  elementType elementType_;

  // parse the .veg file
  VolumetricMeshParser volumetricMeshParser;
  elementType_ = INVALID;

  if (volumetricMeshParser.open(filename) != 0) {
    printf("Error: could not open file %s.\n", filename);
    return elementType_;
  }

  char lineBuffer[1024];
  while (volumetricMeshParser.getNextLine(lineBuffer, 0, 0) != nullptr) {
    // printf("%s\n", lineBuffer);

    // seek for *ELEMENTS
    if (strncmp(lineBuffer, "*ELEMENTS", 9) == 0) {
      // parse element type
      if (volumetricMeshParser.getNextLine(lineBuffer) != nullptr) {
        volumetricMeshParser.removeWhitespace(lineBuffer);

        if (strncmp(lineBuffer, "TET", 3) == 0)
          elementType_ = TET;
        else if (strncmp(lineBuffer, "CUBIC", 5) == 0)
          elementType_ = CUBIC;
        else {
          printf("Error: unknown mesh type %s in file %s\n", lineBuffer, filename);
          return elementType_;
        }
      }
      else {
        printf("Error (getElementType): file %s is not in the .veg format. Offending line:\n%s\n", filename, lineBuffer);
        return elementType_;
      }
    }
  }

  volumetricMeshParser.close();

  if (elementType_ == INVALID)
    printf("Error: could not determine the mesh type in file %s. File may not be in .veg format.\n", filename);

  return elementType_;
}

VolumetricMesh::elementType VolumetricMesh::getElementTypeBinary(const char *filename)
{
  FILE *fin = fopen(filename, "rb");
  if (fin == nullptr) {
    printf("Error in VolumetricMesh::getElementTypeBinary: could not open file %s.\n", filename);
    exit(0);
  }

  elementType elementType_ = getElementType(fin);
  fclose(fin);

  return (elementType_);
}

VolumetricMesh::fileFormatType VolumetricMesh::getFileFormatTypeByExt(const char *filename)
{
  if (BasicAlgorithms::iendWith(filename, ".vegb"))
    return BINARY;
  else if (BasicAlgorithms::iendWith(filename, ".veg"))
    return ASCII;
  return NUM_FILE_FORMATS;
}

VolumetricMesh::elementType VolumetricMesh::getElementType(const char *filename, fileFormatType fileFormat)
{
  if (fileFormat == BY_EXT) {
    fileFormat = getFileFormatTypeByExt(filename);
    if (fileFormat == NUM_FILE_FORMATS) {
      printf("Unknown file extension when loading %s, try ASCII format...\n", filename);
      fileFormat = ASCII;
    }
  }

  switch (fileFormat) {
  case ASCII:
    return VolumetricMesh::getElementTypeASCII(filename);

  case BINARY:
    return VolumetricMesh::getElementTypeBinary(filename);

  default:
    printf("Error: the file format %d is unknown. \n", fileFormat);
    return VolumetricMesh::INVALID;
  }
}

VolumetricMesh::elementType VolumetricMesh::getElementType(void *fin_, int memoryLoad)
{
  unsigned int (*genericRead)(void *, unsigned int, unsigned int, void *);

  void *fin;
  if (memoryLoad) {
    genericRead = &VolumetricMesh::readFromMemory;
    fin = &(fin_);
  }
  else {
    genericRead = &VolumetricMesh::readFromFile;
    fin = fin_;
  }

  double version;
  if ((int)genericRead(&version, sizeof(double), 1, fin) != 1) {
    printf("Error in VolumetricMesh::getElementType: cannot read version.\n");
    throw 0;
  };

  int eleType;
  if ((int)genericRead(&eleType, sizeof(int), 1, fin) != 1) {
    printf("Error in VolumetricMesh::getElementType: cannot read the element type from the stream\n");
    throw 1;
  }

  // rewind the size of an integer and a double
  if (memoryLoad) {
    fin = (void *)((unsigned char *)fin - ((int)sizeof(int)) - ((int)sizeof(double)));
  }
  else {
    int offset = -((int)sizeof(int)) - ((int)sizeof(double));
    fseek((FILE *)fin, offset, SEEK_CUR);
  }

  switch (eleType) {
  case TET:
    return TET;
    break;

  case CUBIC:
    return CUBIC;
    break;

  default:
    printf("Error: could not determine the mesh type from the stream.\n");
    break;
  }

  return INVALID;
}

double VolumetricMesh::getVolume() const
{
  double vol = 0.0;
  for (int el = 0; el < numElements; el++)
    vol += getElementVolume(el);
  return vol;
}

void VolumetricMesh::getVertexVolumes(double *vertexVolumes) const
{
  memset(vertexVolumes, 0, sizeof(double) * numVertices);
  double factor = 1.0 / numElementVertices;
  for (int el = 0; el < numElements; el++) {
    double volume = getElementVolume(el);
    for (int j = 0; j < numElementVertices; j++)
      vertexVolumes[getVertexIndex(el, j)] += factor * volume;
  }
}

Vec3d VolumetricMesh::getElementCenter(int el) const
{
  Vec3d pos(0, 0, 0);
  for (int i = 0; i < numElementVertices; i++)
    pos += getVertex(el, i);

  pos *= 1.0 / numElementVertices;

  return pos;
}

void VolumetricMesh::getVerticesInElements(const std::vector<int> &elements_, std::vector<int> &vertices_) const
{
  vertices_.clear();
  for (unsigned int i = 0; i < elements_.size(); i++)
    for (int j = 0; j < numElementVertices; j++)
      vertices_.push_back(getVertexIndex(elements_[i], j));

  // sort and deduplicate vertices_
  std::sort(vertices_.begin(), vertices_.end());
  auto newEnd = std::unique(vertices_.begin(), vertices_.end());
  vertices_.resize(std::distance(vertices_.begin(), newEnd));
}

void VolumetricMesh::getElementsTouchingVertices(const std::vector<int> &vertices_, std::vector<int> &elements_) const
{
  elements_.clear();
  for (int i = 0; i < numElements; i++) {
    for (int j = 0; j < numElementVertices; j++) {
      // assumes input vector vertices_ is sorted
      if (std::binary_search(vertices_.begin(), vertices_.end(), getVertexIndex(i, j))) {
        elements_.push_back(i);
        break;
      }
    }
  }
}

void VolumetricMesh::getElementsWithOnlyVertices(const std::vector<int> &vertices_, std::vector<int> &elements_) const
{
  elements_.clear();
  for (int i = 0; i < numElements; i++) {
    bool withOnlyVertices = true;
    for (int j = 0; j < numElementVertices; j++) {
      // assumes input vector vertices_ is sorted
      if (std::binary_search(vertices_.begin(), vertices_.end(), getVertexIndex(i, j)) == false) {
        withOnlyVertices = false;
        break;
      }
    }
    if (withOnlyVertices)
      elements_.push_back(i);
  }
}

void VolumetricMesh::getVertexNeighborhood(const std::vector<int> &vertices_, std::vector<int> &neighborhood) const
{
  std::vector<int> elements_;
  getElementsTouchingVertices(vertices_, elements_);
  getVerticesInElements(elements_, neighborhood);
}

double VolumetricMesh::getMass() const
{
  double mass = 0.0;
  for (int i = 0; i < getNumRegions(); i++) {
    const Region &region = getRegion(i);
    double density = getMaterial(region.getMaterialIndex())->getDensity();
    std::set<int> setElements;  // elements in the region
    getSet(region.getSetIndex()).getElements(setElements);

    // over all elements in the region
    for (std::set<int>::iterator iter = setElements.begin(); iter != setElements.end(); iter++) {
      int element = *iter;
      double elementVolume = getElementVolume(element);
      double elementMass = elementVolume * density;
      mass += elementMass;
    }
  }

  return mass;
}

void VolumetricMesh::getInertiaParameters(double &mass, Vec3d &centerOfMass, Mat3d &inertiaTensor) const
{
  mass = 0.0;
  centerOfMass.setZero();
  inertiaTensor.setZero();

  // compute mass, center of mass, inertia tensor
  for (int i = 0; i < numElements; i++) {
    double density = materials[elementMaterial[i]]->getDensity();
    double elementVolume = getElementVolume(i);
    double elementMass = elementVolume * density;

    mass += elementMass;
    Vec3d elementCenter = getElementCenter(i);
    centerOfMass += elementMass * elementCenter;

    Mat3d elementITUnitDensity;
    getElementInertiaTensor(i, elementITUnitDensity);

    double a = elementCenter[0];
    double b = elementCenter[1];
    double c = elementCenter[2];

    Mat3d elementITCorrection = asMat3d(
      b * b + c * c, -a * b, -a * c,
      -a * b, a * a + c * c, -b * c,
      -a * c, -b * c, a * a + b * b);

    Mat3d elementIT = density * elementITUnitDensity + elementMass * elementITCorrection;

    inertiaTensor += elementIT;
  }

  // printf("final mass: %G\n",mass);
  centerOfMass /= mass;

  // correct inertia tensor so it's around the center of mass
  double a = centerOfMass[0];
  double b = centerOfMass[1];
  double c = centerOfMass[2];

  Mat3d correction = asMat3d(
    b * b + c * c, -a * b, -a * c,
    -a * b, a * a + c * c, -b * c,
    -a * c, -b * c, a * a + b * b);

  inertiaTensor -= mass * correction;
}

void VolumetricMesh::getMeshGeometricParameters(Vec3d &centroid, double *radius) const
{
  // compute centroid
  centroid = Vec3d(0, 0, 0);
  for (int i = 0; i < numVertices; i++)
    centroid += getVertex(i);

  centroid /= numVertices;

  // compute radius
  *radius = 0;
  for (int i = 0; i < numVertices; i++) {
    Vec3d vertex = getVertex(i);
    double dist = (vertex - centroid).norm();
    if (dist > *radius)
      *radius = dist;
  }
}

Mesh::BoundingBox VolumetricMesh::getBoundingBox() const
{
  return Mesh::BoundingBox(vertices);
}

int VolumetricMesh::getClosestVertex(Vec3d pos) const
{
  // linear scan
  double closestDist = DBL_MAX;
  int closestVertex = -1;

  for (int i = 0; i < numVertices; i++) {
    const Vec3d &vertexPosition = vertices[i];
    double dist = (pos - vertexPosition).norm();
    if (dist < closestDist) {
      closestDist = dist;
      closestVertex = i;
    }
  }

  return closestVertex;
}

int VolumetricMesh::getClosestElement(const Vec3d &pos) const
{
  // linear scan
  double closestDist = DBL_MAX;
  int closestElement = 0;
  for (int element = 0; element < numElements; element++) {
    Vec3d center = getElementCenter(element);
    double dist = (pos - center).norm();
    if (dist < closestDist) {
      closestDist = dist;
      closestElement = element;
    }
  }

  return closestElement;
}

int VolumetricMesh::getContainingElement(Vec3d pos) const
{
  // linear scan
  for (int element = 0; element < numElements; element++) {
    if (containsVertex(element, pos))
      return element;
  }

  return -1;
}

void VolumetricMesh::setSingleMaterial(double E, double nu, double density)
{
  // erase previous materials
  for (int i = 0; i < numMaterials; i++)
    delete (materials[i]);

  // add a single material
  numMaterials = 1;
  numSets = 1;
  numRegions = 1;

  materials.assign(numMaterials, nullptr);
  sets.assign(numSets, generateAllElementsSet(numElements));
  regions.assign(numRegions, Region(0, 0));

  Material *material = new ENuMaterial("defaultMaterial", density, E, nu);
  materials[0] = material;

  elementMaterial.assign(numElements, 0);
}

void VolumetricMesh::propagateRegionsToElements()
{
  for (int regionIndex = 0; regionIndex < numRegions; regionIndex++) {
    const Region &region = regions[regionIndex];
    int materialIndex = region.getMaterialIndex();

    const std::set<int> &setElements = sets[region.getSetIndex()].getElements();
    for (const auto &elt : setElements)
      elementMaterial[elt] = materialIndex;
  }
}

int VolumetricMesh::generateInterpolationWeights(int numTargetLocations,
  const double *targetLocations, int *elements, int **vertices_, double **weights,
  double zeroThreshold, int verbose) const
{
  // allocate interpolation arrays
  *vertices_ = (int *)malloc(sizeof(int) * numElementVertices * numTargetLocations);
  *weights = (double *)malloc(sizeof(double) * numElementVertices * numTargetLocations);

  double *barycentricWeights = (double *)malloc(sizeof(double) * numElementVertices);

  for (int i = 0; i < numTargetLocations; i++)  // over all interpolation locations
  {
    if ((verbose) && (i % 100 == 0)) {
      printf("%d ", i);
      fflush(nullptr);
    }

    Vec3d pos = Vec3d(targetLocations[3 * i + 0],
      targetLocations[3 * i + 1],
      targetLocations[3 * i + 2]);

    int element = elements[i];
    if (element < 0) {
      printf("Error: invalid element index %d.\n", element);
      return 1;
    }

    computeBarycentricWeights(element, pos, barycentricWeights);

    if (zeroThreshold > 0) {
      // check whether vertex is close enough to the mesh
      double minDistance = DBL_MAX;
      int numElementVertices = getNumElementVertices();
      int assignedZero = 0;
      for (int ii = 0; ii < numElementVertices; ii++) {
        const Vec3d &vpos = getVertex(element, ii);
        if ((vpos - pos).norm() < minDistance) {
          minDistance = (vpos - pos).norm();
        }
      }

      if (minDistance > zeroThreshold) {
        // assign zero weights
        for (int ii = 0; ii < numElementVertices; ii++)
          barycentricWeights[ii] = 0.0;
        assignedZero++;
        continue;
      }
    }

    for (int ii = 0; ii < numElementVertices; ii++) {
      (*vertices_)[numElementVertices * i + ii] = getVertexIndex(element, ii);
      (*weights)[numElementVertices * i + ii] = barycentricWeights[ii];
    }
  }

  free(barycentricWeights);

  return 0;
}

int VolumetricMesh::generateContainingElements(int numTargetLocations, const double *targetLocations, int **elements, int useClosestElementIfOutside, int verbose) const
{
  int numExternalVertices = 0;

  (*elements) = (int *)malloc(sizeof(int) * numTargetLocations);

  // determine containing (or closest) elements
  for (int i = 0; i < numTargetLocations; i++)  // over all interpolation locations
  {
    Vec3d pos = Vec3d(targetLocations[3 * i + 0], targetLocations[3 * i + 1], targetLocations[3 * i + 2]);

    // find element containing pos
    int element = getContainingElement(pos);

    // use closest element if outside
    if (useClosestElementIfOutside && (element < 0)) {
      element = getClosestElement(pos);
      numExternalVertices++;
    }

    (*elements)[i] = element;
  }

  return numExternalVertices;
}

int VolumetricMesh::generateInterpolationWeights(int numTargetLocations, const double *targetLocations, int **vertices_, double **weights, double zeroThreshold, int **containingElements, int verbose) const
{
  int *elements = nullptr;

  int useClosestElementIfOutside = 1;
  int numExternalVertices = generateContainingElements(numTargetLocations, targetLocations, &elements, useClosestElementIfOutside);

  // allocate interpolation arrays
  *vertices_ = (int *)malloc(sizeof(int) * numElementVertices * numTargetLocations);
  *weights = (double *)malloc(sizeof(double) * numElementVertices * numTargetLocations);
  int code = generateInterpolationWeights(numTargetLocations, targetLocations, elements, vertices_, weights, zeroThreshold, verbose);

  if (containingElements == nullptr)
    free(elements);
  else
    *containingElements = elements;

  return (code == 0) ? numExternalVertices : -1;
}

int VolumetricMesh::getNumInterpolationElementVertices(const char *filename)
{
  FILE *fin = fopen(filename, "r");
  if (!fin) {
    printf("Error: unable to open file %s.\n", filename);
    return -1;
  }

  char s[1024];
  if (fgets(s, 1024, fin) == nullptr) {
    printf("Error: incorrect first line of file %s.\n", filename);
    return -2;
  }
  fclose(fin);

  VolumetricMeshParser::beautifyLine(s, 1);

  int slen = strlen(s);
  int count = 0;
  for (int i = 0; i < slen; i++)
    if (s[i] == ' ')
      count++;

  if (count % 2 == 1) {
    printf("Error: odd number of whitespaces in the first line of file %s.\n", filename);
    return -3;
  }

  return count / 2;
}

int VolumetricMesh::loadInterpolationWeights(const char *filename, int numTargetLocations, int numElementVertices_, int **vertices_, double **weights)
{
  FILE *fin = fopen(filename, "r");
  if (!fin) {
    printf("Error: unable to open file %s.\n", filename);
    return 2;
  }

  // allocate interpolation arrays
  *vertices_ = (int *)malloc(sizeof(int) * numElementVertices_ * numTargetLocations);
  *weights = (double *)malloc(sizeof(double) * numElementVertices_ * numTargetLocations);

  int numReadTargetLocations = -1;
  int currentVertex;

  // read the elements one by one and accumulate entries
  while (numReadTargetLocations < numTargetLocations - 1) {
    numReadTargetLocations++;

    if (feof(fin)) {
      printf("Error: interpolation file is too short. Num vertices in interp file: %d . Should be: %d .\n", numReadTargetLocations, numTargetLocations);
      free(*vertices_);
      free(*weights);
      *vertices_ = nullptr;
      *weights = nullptr;
      fclose(fin);
      return 1;
    }

    if (fscanf(fin, "%d", &currentVertex) < 1)
      printf("Warning: bad file syntax. Unable to read interpolation info.\n");

    if (currentVertex != numReadTargetLocations) {
      printf("Error: consecutive vertex index at position %d mismatch.\n", currentVertex);
      free(*vertices_);
      free(*weights);
      *vertices_ = nullptr;
      *weights = nullptr;
      fclose(fin);
      return 1;
    }

    for (int j = 0; j < numElementVertices_; j++) {
      if (fscanf(fin, "%d %lf", &((*vertices_)[currentVertex * numElementVertices_ + j]), &((*weights)[currentVertex * numElementVertices_ + j])) < 2)
        printf("Warning: bad file syntax. Unable to read interpolation info.\n");
    }

    if (fscanf(fin, "\n") < 0) {
      // printf("Warning: bad file syntax. Missing end of line in the interpolation file.\n");
      // do nothing
    }
  }

  fclose(fin);
  return 0;
}

int VolumetricMesh::loadInterpolationWeightsBinary(const char *filename, int *numTargetLocations, int *numElementVertices_, int **vertices_, double **weights)
{
  FILE *fin = fopen(filename, "rb");
  if (!fin) {
    printf("Error: unable to open file %s.\n", filename);
    return 2;
  }

  int code = loadInterpolationWeightsBinary(fin, numTargetLocations, numElementVertices_, vertices_, weights);
  fclose(fin);

  if (code != 0)
    printf("Error reading from file %s.\n", filename);

  return code;
}

int VolumetricMesh::loadInterpolationWeightsBinary(FILE *fin, int *numTargetLocations, int *numElementVertices_, int **vertices_, double **weights)
{
  int buffer[2];
  int readItems = (int)fread(buffer, sizeof(int), 2, fin);
  if (readItems < 2)
    return 1;

  *numTargetLocations = buffer[0];
  *numElementVertices_ = buffer[1];

  // allocate interpolation arrays
  *vertices_ = (int *)malloc(sizeof(int) * *numElementVertices_ * *numTargetLocations);
  *weights = (double *)malloc(sizeof(double) * *numElementVertices_ * *numTargetLocations);

  readItems = (int)fread(*vertices_, sizeof(int), *numElementVertices_ * *numTargetLocations, fin);
  if (readItems < *numElementVertices_ * *numTargetLocations)
    return 1;

  readItems = (int)fread(*weights, sizeof(double), *numElementVertices_ * *numTargetLocations, fin);
  if (readItems < *numElementVertices_ * *numTargetLocations)
    return 1;

  return 0;
}

int VolumetricMesh::saveInterpolationWeights(const char *filename, int numTargetLocations, int numElementVertices_, const int *vertices_, const double *weights)
{
  FILE *fout = fopen(filename, "w");
  if (!fout) {
    printf("Error: unable to open file %s.\n", filename);
    return 1;
  }

  for (int currentVertex = 0; currentVertex < numTargetLocations; currentVertex++) {
    fprintf(fout, "%d", currentVertex);

    for (int j = 0; j < numElementVertices_; j++)
      fprintf(fout, " %d %lf", vertices_[currentVertex * numElementVertices_ + j],
        weights[currentVertex * numElementVertices_ + j]);

    fprintf(fout, "\n");
  }

  fclose(fout);
  return 0;
}

int VolumetricMesh::saveInterpolationWeightsBinary(const char *filename, int numTargetLocations, int numElementVertices_, const int *vertices_, const double *weights)
{
  FILE *fout = fopen(filename, "wb");
  if (!fout) {
    printf("Error: unable to open file %s.\n", filename);
    return 1;
  }

  int code = saveInterpolationWeightsBinary(fout, numTargetLocations, numElementVertices_, vertices_, weights);
  fclose(fout);
  if (code != 0)
    printf("Error reading from file %s.\n", filename);
  return code;
}

int VolumetricMesh::saveInterpolationWeightsBinary(FILE *fout, int numTargetLocations, int numElementVertices_, const int *vertices_, const double *weights)
{
  int buffer[2];
  buffer[0] = numTargetLocations;
  buffer[1] = numElementVertices_;

  int writtenItems = (int)fwrite(buffer, sizeof(int), 2, fout);
  if (writtenItems < 2)
    return 1;

  writtenItems = (int)fwrite(vertices_, sizeof(int), numElementVertices_ * numTargetLocations, fout);
  if (writtenItems < numElementVertices_ * numTargetLocations)
    return 1;

  writtenItems = (int)fwrite(weights, sizeof(double), numElementVertices_ * numTargetLocations, fout);
  if (writtenItems < numElementVertices_ * numTargetLocations)
    return 1;

  return 0;
}

void VolumetricMesh::interpolate(const double *u, double *uTarget, int numTargetLocations, int numElementVertices_, const int *vertices_, const double *weights)
{
  tbb::parallel_for(0, numTargetLocations, [&](int i) {
    Vec3d defo(0, 0, 0);
    for (int j = 0; j < numElementVertices_; j++) {
      int volumetricMeshVertexIndex = vertices_[numElementVertices_ * i + j];
      defo += weights[numElementVertices_ * i + j] * asVec3d(u + 3 * volumetricMeshVertexIndex);
    }
    (Eigen::Map<Vec3d>(uTarget + 3 * i)) = defo;
  });
}

int VolumetricMesh::interpolateGradient(const double *U, int numFields, Vec3d pos, double *grad) const
{
  // find the element containing "pos"
  int externalVertex = 0;
  int element = getContainingElement(pos);
  if (element < 0) {
    element = getClosestElement(pos);
    externalVertex = 1;
  }

  interpolateGradient(element, U, numFields, pos, grad);

  return externalVertex;
}

void VolumetricMesh::exportMeshGeometry(int *numVertices_, double **vertices_, int *numElements_, int *numElementVertices_, int **elements_) const
{
  if (numVertices_ != nullptr)
    *numVertices_ = numVertices;
  if (numElements_ != nullptr)
    *numElements_ = numElements;
  if (numElementVertices_ != nullptr)
    *numElementVertices_ = numElementVertices;

  if (vertices_ != nullptr) {
    *vertices_ = (double *)malloc(sizeof(double) * 3 * numVertices);
    for (int i = 0; i < numVertices; i++) {
      const Vec3d &v = getVertex(i);
      (*vertices_)[3 * i + 0] = v[0];
      (*vertices_)[3 * i + 1] = v[1];
      (*vertices_)[3 * i + 2] = v[2];
    }
  }

  if (elements_ != nullptr) {
    *elements_ = (int *)malloc(sizeof(int) * numElementVertices * numElements);
    memcpy(*elements_, elements.data(), sizeof(int) * numElementVertices * numElements);
  }
}

void VolumetricMesh::exportMeshGeometry(std::vector<Vec3d> &vertexBuffer, std::vector<int> &elementsBuffer) const
{
  vertexBuffer = vertices;
  elementsBuffer = elements;
}

void VolumetricMesh::exportMeshGeometry(std::vector<Vec3d> &vertexBuffer) const
{
  vertexBuffer = vertices;
}

void VolumetricMesh::computeGravity(double *gravityForce, double g, bool addForce) const
{
  if (!addForce)
    memset(gravityForce, 0, sizeof(double) * 3 * numVertices);

  double invNumElementVertices = 1.0 / getNumElementVertices();

  for (int el = 0; el < numElements; el++) {
    double volume = getElementVolume(el);
    double density = getElementDensity(el);
    double mass = density * volume;
    for (int j = 0; j < getNumElementVertices(); j++)
      gravityForce[3 * getVertexIndex(el, j) + 1] -= invNumElementVertices * mass * g;  // gravity assumed to act in negative y-direction
  }
}

void VolumetricMesh::applyDeformation(const double *u)
{
  for (int i = 0; i < numVertices; i++) {
    Vec3d &v = getVertex(i);
    v[0] += u[3 * i + 0];
    v[1] += u[3 * i + 1];
    v[2] += u[3 * i + 2];
  }
}

// transforms every vertex as X |--> pos + R * X
void VolumetricMesh::applyLinearTransformation(double *pos, double *R)
{
  for (int i = 0; i < numVertices; i++) {
    Vec3d &v = getVertex(i);

    double newPos[3];
    for (int j = 0; j < 3; j++) {
      newPos[j] = pos[j];
      for (int k = 0; k < 3; k++)
        newPos[j] += R[3 * j + k] * v[k];
    }

    v[0] = newPos[0];
    v[1] = newPos[1];
    v[2] = newPos[2];
  }
}

void VolumetricMesh::setMaterial(int i, const Material *material)
{
  delete (materials[i]);
  materials[i] = material->clone();
}

VolumetricMesh::VolumetricMesh(const VolumetricMesh &volumetricMesh, int numElements_, int *elements_, std::map<int, int> *vertexMap_)
{
  // determine vertices in the submesh
  numElementVertices = volumetricMesh.getNumElementVertices();
  std::set<int> vertexSet;
  for (int i = 0; i < numElements_; i++)
    for (int j = 0; j < numElementVertices; j++)
      vertexSet.insert(volumetricMesh.getVertexIndex(elements_[i], j));

  // copy vertices into place and also into vertexMap
  numVertices = vertexSet.size();
  vertices.resize(numVertices);
  std::set<int>::iterator iter;
  int vertexNo = 0;
  std::map<int, int> vertexMap;
  for (iter = vertexSet.begin(); iter != vertexSet.end(); iter++) {
    vertices[vertexNo] = volumetricMesh.getVertex(*iter);
    vertexMap.insert(std::make_pair(*iter, vertexNo));
    vertexNo++;
  }

  if (vertexMap_ != nullptr)
    *vertexMap_ = vertexMap;

  // copy elements
  numElements = numElements_;
  elements.resize(numElements * numElementVertices);
  elementMaterial.resize(numElements);
  std::map<int, int> elementMap;
  for (int i = 0; i < numElements; i++) {
    for (int j = 0; j < numElementVertices; j++) {
      std::map<int, int>::iterator iter2 = vertexMap.find(volumetricMesh.elements[elements_[i] * numElementVertices + j]);
      if (iter2 == vertexMap.end()) {
        printf("Internal error 1.\n");
        throw 1;
      }
      elements[i * numElementVertices + j] = iter2->second;
    }

    elementMaterial[i] = (volumetricMesh.elementMaterial)[elements_[i]];
    elementMap.insert(std::make_pair(elements_[i], i));
  }

  // copy materials
  numMaterials = volumetricMesh.getNumMaterials();
  numSets = volumetricMesh.getNumSets();
  numRegions = volumetricMesh.getNumRegions();

  materials.assign(numMaterials, nullptr);
  for (int i = 0; i < numMaterials; i++)
    materials[i] = volumetricMesh.getMaterial(i)->clone();

  // copy element sets; restrict element sets to the new mesh, also rename vertices to reflect new vertex indices
  std::vector<Set> newSets;
  std::map<int, int> oldToNewSetIndex;
  for (int oldSetIndex = 0; oldSetIndex < volumetricMesh.getNumSets(); oldSetIndex++) {
    const Set &oldSet = volumetricMesh.getSet(oldSetIndex);
    std::set<int> oldElements;
    oldSet.getElements(oldElements);

    for (std::set<int>::iterator iter = oldElements.begin(); iter != oldElements.end(); iter++) {
      if (*iter < 0) {
        printf("Internal error 2.\n");
        exit(1);
      }
    }

    // construct the element list
    std::vector<int> newElements;
    for (std::set<int>::iterator iter = oldElements.begin(); iter != oldElements.end(); iter++) {
      std::map<int, int>::iterator iter2 = elementMap.find(*iter);
      if (iter2 != elementMap.end())
        newElements.push_back(iter2->second);
    }

    // if there is at least one element in the new set, create a set for it
    if (newElements.size() > 0) {
      Set newSet(oldSet.getName());
      for (unsigned int j = 0; j < newElements.size(); j++) {
        if (newElements[j] < 0) {
          printf("Internal error 3.\n");
          exit(1);
        }
        newSet.insert(newElements[j]);
      }
      newSets.emplace_back(std::move(newSet));
      oldToNewSetIndex.insert(std::make_pair(oldSetIndex, (int)newSets.size() - 1));
    }
  }

  numSets = newSets.size();
  sets = std::move(newSets);

  // printf("numSets: %d\n", numSets);

  // copy regions; remove empty ones
  std::vector<Region> vregions;
  for (int i = 0; i < numRegions; i++) {
    const Region &sregion = volumetricMesh.getRegion(i);
    std::map<int, int>::iterator iter = oldToNewSetIndex.find(sregion.getSetIndex());
    if (iter != oldToNewSetIndex.end()) {
      vregions.emplace_back(sregion.getMaterialIndex(), iter->second);
    }
  }

  numRegions = vregions.size();
  regions = move(vregions);

  // sanity check
  // seek each element in all the regions
  for (int el = 0; el < numElements; el++) {
    int found = 0;
    for (int region = 0; region < numRegions; region++) {
      int elementSet = regions[region].getSetIndex();

      // seek for element in elementSet
      if (sets[elementSet].isMember(el)) {
        if (found != 0)
          printf("Warning: element %d (1-indexed) is in more than one region.\n", el + 1);
        else
          found = 1;
      }
    }
    if (found == 0)
      printf("Warning: element %d (1-indexed) is not in any of the regions.\n", el + 1);
  }

  // sanity check: make sure all elements are between bounds
  for (int i = 0; i < numSets; i++) {
    std::set<int> elts;
    sets[i].getElements(elts);
    for (std::set<int>::iterator iter = elts.begin(); iter != elts.end(); iter++) {
      if (*iter < 0)
        printf("Warning: encountered negative element index in element set %d.\n", i);
      if (*iter >= numElements)
        printf("Warning: encountered too large element index in element set %d.\n", i);
    }
  }
}

// if vertexMap is non-null, it also returns a renaming datastructure: vertexMap[big mesh vertex] is the vertex index in the subset mesh
void VolumetricMesh::setToSubsetMesh(const std::set<int> &subsetElements, int toRemoveIsolatedVertices, std::map<int, int> *vertexMap)
{
  // oldEleID -> newEleID
  std::vector<int> lookupTable(numElements, -1);

  int head = 0;
  for (int tail = 0; tail < numElements; tail++) {
    if (subsetElements.find(tail) != subsetElements.end()) {
      if (head != tail) {
        memcpy(&elements[head * numElementVertices], &elements[tail * numElementVertices], sizeof(int) * numElementVertices);
        elementMaterial[head] = elementMaterial[tail];
      }
      lookupTable[tail] = head;  // update to new index
      head++;
    }
  }
  PGO_ALOG(head == (int)subsetElements.size());

  numElements = subsetElements.size();
  elements.resize(numElements * numElementVertices);
  elementMaterial.resize(numElements);

  for (int setIndex = 0; setIndex < numSets; setIndex++) {
    std::set<int> oldSet = move(sets[setIndex].getElements());
    sets[setIndex].clear();
    for (int eleID : oldSet) {
      if (subsetElements.find(eleID) == subsetElements.end())  // not found!!
        continue;
      int newIndex = lookupTable[eleID];
      sets[setIndex].insert(newIndex);
    }
  }

  if (toRemoveIsolatedVertices) {
    removeIsolatedVertices(vertexMap);
  }

  // note: empty sets and regions are not removed
}

void VolumetricMesh::removeIsolatedVertices(std::map<int, int> *old2NewVertexIDMap)
{
  std::vector<int> retainedVertices;
  for (int el = 0; el < numElements; el++)
    for (int j = 0; j < numElementVertices; j++)
      retainedVertices.push_back(getVertexIndex(el, j));
  // sort and deduplicate retainedVertices
  std::sort(retainedVertices.begin(), retainedVertices.end());
  auto newEnd = std::unique(retainedVertices.begin(), retainedVertices.end());
  retainedVertices.resize(std::distance(retainedVertices.begin(), newEnd));

  // oldvtxID -> newVtxID
  std::vector<int> old2newVtxID(numVertices, -1);
  if (old2NewVertexIDMap != nullptr)
    old2NewVertexIDMap->clear();

  int head = 0;
  for (int tail = 0; tail < numVertices; tail++) {
    if (std::binary_search(retainedVertices.begin(), retainedVertices.end(), tail)) {
      old2newVtxID[tail] = head;
      if (old2NewVertexIDMap != nullptr)
        old2NewVertexIDMap->emplace(tail, head);
      if (head != tail)
        vertices[head] = vertices[tail];
      head++;
    }
  }
  PGO_ALOG(head == int(retainedVertices.size()));

  // rename vertices inside the elements
  for (int el = 0, size = numElements * numElementVertices; el < size; el++)
    elements[el] = old2newVtxID[elements[el]];

  numVertices = retainedVertices.size();
  vertices.resize(numVertices);
}

int VolumetricMesh::exportToEle(const char *baseFilename, int includeRegions) const
{
  char s[1024];
  sprintf(s, "%s.ele", baseFilename);

  FILE *fout = fopen(s, "w");
  if (!fout) {
    printf("Error: could not write to %s.\n", s);
    return 1;
  }

  int *elementRegion = nullptr;
  if (includeRegions) {
    elementRegion = (int *)malloc(sizeof(int) * getNumElements());
    for (int el = 0; el < getNumElements(); el++) {
      int found = 0;
      for (int region = 0; region < numRegions; region++) {
        int elementSet = (regions[region]).getSetIndex();

        // seek for element in elementSet
        if (sets[elementSet].isMember(el)) {
          if (found != 0)
            printf("Warning: element %d (1-indexed) is in more than one region.\n", el + 1);
          else
            found = region + 1;
        }
      }
      if (found == 0)
        printf("Warning: element %d (1-indexed) is not in any of the regions.\n", el + 1);
      elementRegion[el] = found;
    }
  }

  if (includeRegions)
    fprintf(fout, "%d %d %d\n", numElements, numElementVertices, 1);
  else
    fprintf(fout, "%d %d %d\n", numElements, numElementVertices, 0);

  for (int el = 0; el < numElements; el++) {
    fprintf(fout, "%d ", el + 1);
    for (int j = 0; j < numElementVertices; j++) {
      fprintf(fout, "%d", getVertexIndex(el, j) + 1);
      if (j != numElementVertices - 1)
        fprintf(fout, " ");
    }
    if (includeRegions) {
      fprintf(fout, " %d", elementRegion[el]);
    }
    fprintf(fout, "\n");
  }

  fprintf(fout, "# generated by the volumetricMesh class\n");

  fclose(fout);

  if (includeRegions)
    free(elementRegion);

  sprintf(s, "%s.node", baseFilename);

  fout = fopen(s, "w");
  if (!fout) {
    printf("Error: could not write to %s.\n", s);
    return 1;
  }

  fprintf(fout, "%d %d %d %d\n", numVertices, 3, 0, 0);
  for (int v = 0; v < numVertices; v++) {
    fprintf(fout, "%d ", v + 1);
    const Vec3d &vtx = getVertex(v);
    fprintf(fout, "%.15f %.15f %.15f\n", vtx[0], vtx[1], vtx[2]);
  }

  fprintf(fout, "# generated by the volumetricMesh class\n");

  fclose(fout);

  return 0;
}

void VolumetricMesh::renumberVertices(const std::vector<int> &permutation)
{
  // renumber vertices
  std::vector<Vec3d> newVertices(numVertices);
  for (int i = 0; i < numVertices; i++)
    newVertices[permutation[i]] = vertices[i];
  vertices = std::move(newVertices);

  // renumber tets
  for (int &vtxID : elements)
    vtxID = permutation[vtxID];
}

void VolumetricMesh::addMaterial(const Material *material, const Set &newSet, bool removeEmptySets, bool removeEmptyMaterials)
{
  // add new material to materials
  numMaterials++;
  materials.resize(numMaterials);
  materials[numMaterials - 1] = material->clone();

  // remove indices in the sets that belong to the newSet
  const std::set<int> &newElements = newSet.getElements();
  std::vector<bool> eleCovered(numElements, false);
  for (int i = 1; i < numSets; i++)  // skip the first set, which is always allElements
  {
    std::set<int> &s = sets[i].getElements();

    for (std::set<int>::iterator it = s.begin(); it != s.end();) {
      int ele = *it;
      eleCovered[ele] = true;
      if (newElements.find(ele) != newElements.end()) {
        std::set<int>::iterator it2 = it;
        it++;
        s.erase(it2);
      }
      else
        it++;
    }
  }

  // we have to be careful here: if previously the entire mesh is in default set: allElements,
  // adding a new set will invalidate the elements that are previously in allElements
  // a newSet is needed to cover those elements
  std::set<int> restSet;
  for (int i = 0; i < numElements; i++)
    if (eleCovered[i] == false)  // this element is covered only by allElements
      restSet.insert(i);
  if (restSet.size() > 0)
    numSets++;

  // add the new Set
  numSets++;
  sets.resize(numSets);
  sets[numSets - 1] = Set(newSet);
  if (restSet.size() > 0)
    sets[numSets - 2] = Set("restElements", restSet);

  // create a new Region
  numRegions++;
  regions.resize(numRegions);
  regions[numRegions - 1] = Region(numMaterials - 1, numSets - 1);

  if (restSet.size() > 0) {
    // first find the material that used for allElements
    for (int i = 0; i < numRegions - 1; i++)
      if (regions[i].getSetIndex() == 0)  // this is the index of the allElements set
        regions[i].setSetIndex(numSets - 2);
  }

  // modify elementMaterial
  for (std::set<int>::const_iterator it = newElements.begin(); it != newElements.end(); it++) {
    int el = *it;
    PGO_ALOG(el >= 0 && el < numElements);
    elementMaterial[el] = numMaterials - 1;
  }

  if (removeEmptySets) {
    bool hasEmptySet = false;
    std::vector<int> setIndexChange(numSets, 0);  // old set index -> new set index
    int newIndex = 0;                             // store the next available set index
    for (int i = 0; i < numSets; i++) {
      if (sets[i].getNumElements() == 0) {
        setIndexChange[i] = -1;  // this set will be deleted, so its new set index is -1
        hasEmptySet = true;
      }
      else {
        setIndexChange[i] = newIndex;  // this set is remained, its new index is the next available index
        if (newIndex != i)             // this means there's already at least one set deleted
        {
          PGO_ALOG(newIndex < i);
          sets[newIndex] = sets[i];  // assign the pointer to the current set to the location at newIndex
        }
        newIndex++;
      }
    }

    if (hasEmptySet) {
      PGO_ALOG(newIndex < numSets);
      numSets = newIndex;
      sets.resize(numSets);

      int newRegionIdx = 0;
      for (int i = 0; i < numRegions; i++) {
        int oldSetIdx = regions[i].getSetIndex();
        PGO_ALOG((size_t)oldSetIdx < setIndexChange.size());
        if (setIndexChange[oldSetIdx] != -1)  // this set has been deleted
        {
          regions[i].setSetIndex(setIndexChange[oldSetIdx]);
          if (newRegionIdx != i)
            regions[newRegionIdx] = regions[i];
          newRegionIdx++;
        }
      }
      numRegions = newRegionIdx;
      regions.resize(numRegions);
    }  // end if (hasEmptySet)
    else {
      PGO_ALOG(newIndex == numSets);
    }
  }  // end if (removeEmptySets)

  // remove material
  if (removeEmptyMaterials) {
    // count #Elements for each material
    std::vector<int> elementsWithMaterial(numMaterials, 0);

    for (int i = 0; i < numElements; i++) {
      int matIdx = elementMaterial[i];
      PGO_ALOG(matIdx >= 0 && matIdx < numMaterials);
      elementsWithMaterial[matIdx]++;
    }

    int newMatIdx = 0;
    std::vector<int> matIndexChange(numMaterials, 0);  // old material index -> new material index
    bool hasEmptyMat = false;
    for (int i = 0; i < numMaterials; i++) {
      if (elementsWithMaterial[i] == 0) {
        matIndexChange[i] = -1;
        delete materials[i];
        materials[i] = nullptr;
        hasEmptyMat = true;
      }
      else {
        matIndexChange[i] = newMatIdx;
        if (newMatIdx != i)
          materials[newMatIdx] = materials[i];
        newMatIdx++;
      }
    }

    if (hasEmptyMat) {
      numMaterials = newMatIdx;
      materials.resize(numMaterials);

      // we also need to modify and delete invalid regions
      bool hasInvalidRegion = false;
      int newRegionIdx = 0;
      for (int i = 0; i < numRegions; i++) {
        int oldMatIndex = regions[i].getMaterialIndex();
        if (matIndexChange[oldMatIndex] < 0) {
          hasInvalidRegion = true;
        }
        else {
          regions[i].setMaterialIndex(matIndexChange[oldMatIndex]);
          if (newRegionIdx != i)
            regions[newRegionIdx] = regions[i];
          newRegionIdx++;
        }
      }

      if (hasInvalidRegion) {
        numRegions = newRegionIdx;
        regions.resize(numRegions);
      }

      // reassign the correct material index to each element
      propagateRegionsToElements();
    }
  }  // end if (removeEmptyMaterials)
}

VolumetricMesh::Set VolumetricMesh::generateAllElementsSet(int numElements)
{
  Set set(allElementsSetName);
  for (int i = 0; i < numElements; i++)
    set.insert(i);
  return set;
}

}  // namespace VolumetricMeshes
}  // namespace pgo
