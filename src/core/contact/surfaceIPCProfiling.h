#pragma once

#include <string_view>

namespace pgo::Contact::SurfaceIPCProfileSections
{

inline constexpr std::string_view kPairBuildStatic = "contact.surface.pair_build.static";
inline constexpr std::string_view kPairBuildSwept = "contact.surface.pair_build.swept";
inline constexpr std::string_view kMaxStepPT = "contact.surface.max_step_pt";
inline constexpr std::string_view kMaxStepEE = "contact.surface.max_step_ee";
inline constexpr std::string_view kEnergy = "contact.surface.energy";
inline constexpr std::string_view kGradient = "contact.surface.gradient";
inline constexpr std::string_view kHessian = "contact.surface.hessian";
inline constexpr std::string_view kCombined = "contact.surface.combined";
inline constexpr std::string_view kWrapperSync = "contact.wrapper.sync";
inline constexpr std::string_view kFloorPostPass = "contact.wrapper.floor_post_pass";

}  // namespace pgo::Contact::SurfaceIPCProfileSections
