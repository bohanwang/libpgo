/*
copyright to Bohan Wang
*/

#pragma once

namespace pgo
{
namespace Contact
{
namespace CIPC
{

struct PTPair
{
  int p, t0, t1, t2;  // vertex indices
  double weight;      // area(point) * area(triangle)
};

struct EEPair
{
  int ea0, ea1, eb0, eb1;  // vertex indices
  double weight;           // length(edgeA) * length(edgeB)
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
