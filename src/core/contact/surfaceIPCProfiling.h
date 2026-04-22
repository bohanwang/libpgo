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
inline constexpr std::string_view kAdapterFunc = "contact.adapter.func";
inline constexpr std::string_view kAdapterGradient = "contact.adapter.gradient";
inline constexpr std::string_view kAdapterHessianDirect = "contact.adapter.hessian_direct";
inline constexpr std::string_view kAdapterMaxStep = "contact.adapter.max_step";
inline constexpr std::string_view kAdapterMapToSurface = "contact.adapter.map_to_surface";
inline constexpr std::string_view kAdapterPullbackGradient = "contact.adapter.pullback_gradient";
inline constexpr std::string_view kAdapterPullbackHessian = "contact.adapter.pullback_hessian";

}  // namespace pgo::Contact::SurfaceIPCProfileSections
