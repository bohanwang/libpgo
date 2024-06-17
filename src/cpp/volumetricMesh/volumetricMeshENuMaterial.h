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
#pragma once

#include "volumetricMesh.h"

namespace pgo
{
namespace VolumetricMeshes
{

// stores an isotropic material specified by E (Young's modulus), nu (Poisson's ratio), and density
// such a material specification is very common: (corotational) linear FEM, StVK, etc.
class VolumetricMesh::ENuMaterial : public VolumetricMesh::Material
{
public:
  ENuMaterial(std::string name, double density = density_default, double E = E_default, double nu = nu_default);
  ENuMaterial(const ENuMaterial &eNuMaterial);
  virtual ~ENuMaterial() {}
  virtual Material *clone() const override;
  virtual Material::materialType getType() const override { return Material::ENU; }

  inline double getE() const;       // Young's modulus
  inline double getNu() const;      // Poisson's ratio
  inline double getLambda() const;  // Lame's lambda coefficient
  inline double getMu() const;      // Lame's mu coefficient
  inline void setE(double E);
  inline void setNu(double nu);

protected:
  double E_, nu_;
};

inline VolumetricMesh::ENuMaterial::ENuMaterial(std::string name, double density, double E, double nu):
  VolumetricMesh::Material(name, density), E_(E), nu_(nu)
{
}
inline VolumetricMesh::ENuMaterial::ENuMaterial(const ENuMaterial &eNuMaterial):
  VolumetricMesh::Material(eNuMaterial.getName(), eNuMaterial.getDensity()), E_(eNuMaterial.getE()), nu_(eNuMaterial.getNu())
{
}
inline double VolumetricMesh::ENuMaterial::getE() const
{
  return E_;
}
inline double VolumetricMesh::ENuMaterial::getNu() const
{
  return nu_;
}
inline double VolumetricMesh::ENuMaterial::getLambda() const
{
  return (nu_ * E_) / ((1 + nu_) * (1 - 2 * nu_));
}
inline double VolumetricMesh::ENuMaterial::getMu() const
{
  return E_ / (2 * (1 + nu_));
}
inline void VolumetricMesh::ENuMaterial::setE(double E)
{
  E_ = E;
}
inline void VolumetricMesh::ENuMaterial::setNu(double nu)
{
  nu_ = nu;
}

// obtain pointer to ENuMaterial (necessary inside classes that assume ENu material)
// performs a check via getType and returns nullptr if material is not ENU
VolumetricMesh::ENuMaterial *downcastENuMaterial(VolumetricMesh::Material *material);
const VolumetricMesh::ENuMaterial *downcastENuMaterial(const VolumetricMesh::Material *material);
}  // namespace VolumetricMeshes
}  // namespace pgo