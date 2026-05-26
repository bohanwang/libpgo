/*
copyright to Bohan Wang
*/

#pragma once

namespace pgo
{
namespace Contact
{
namespace IPC
{

namespace barrier
{

double b(double s, double shat);
double dbds(double s, double shat);
double d2bds2(double s, double shat);

}  // namespace barrier

}  // namespace IPC
}  // namespace Contact
}  // namespace pgo
