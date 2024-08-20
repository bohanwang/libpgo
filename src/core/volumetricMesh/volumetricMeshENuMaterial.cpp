/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "volumetricMesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
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

#include "volumetricMeshENuMaterial.h"

namespace pgo
{
namespace VolumetricMeshes
{
VolumetricMesh::Material *VolumetricMesh::ENuMaterial::clone() const
{
  return new VolumetricMesh::ENuMaterial(*this);
}

// performs a check via getType and returns NULL if material is not ENU
VolumetricMesh::ENuMaterial *downcastENuMaterial(VolumetricMesh::Material *material)
{
  if (material->getType() != VolumetricMesh::Material::ENU)
    return nullptr;

  return dynamic_cast<VolumetricMesh::ENuMaterial *>(material);
}

const VolumetricMesh::ENuMaterial *downcastENuMaterial(const VolumetricMesh::Material *material)
{
  if (material->getType() != VolumetricMesh::Material::ENU)
    return nullptr;

  return dynamic_cast<const VolumetricMesh::ENuMaterial *>(material);
}

}  // namespace VolumetricMeshes
}  // namespace pgo
