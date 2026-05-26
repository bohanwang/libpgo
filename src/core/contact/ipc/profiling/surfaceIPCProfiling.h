#pragma once

#include <string_view>

namespace pgo::Contact::SurfaceIPCProfileSections
{

inline constexpr std::string_view kPairBuildStatic = "contact.surface.pair_build.static";
inline constexpr std::string_view kPairBuildSwept = "contact.surface.pair_build.swept";
inline constexpr std::string_view kPairBuildSelfAABB = "contact.surface.pair_build.self_aabb";
inline constexpr std::string_view kPairBuildSelfPTHashQuery = "contact.surface.pair_build.self_pt_hash_query";
inline constexpr std::string_view kPairBuildSelfPTHashInsert = "contact.surface.pair_build.self_pt_hash_insert";
inline constexpr std::string_view kPairBuildSelfPTQuery = "contact.surface.pair_build.self_pt_query";
inline constexpr std::string_view kPairBuildSelfPTHashCandidates = "contact.surface.pair_build.self_pt.hash_candidates";
inline constexpr std::string_view kPairBuildSelfPTDistanceTests = "contact.surface.pair_build.self_pt.distance_tests";
inline constexpr std::string_view kPairBuildSelfPTAcceptedPairs = "contact.surface.pair_build.self_pt.accepted_pairs";
inline constexpr std::string_view kPairBuildSelfEEHashQuery = "contact.surface.pair_build.self_ee_hash_query";
inline constexpr std::string_view kPairBuildSelfEEHashInsert = "contact.surface.pair_build.self_ee_hash_insert";
inline constexpr std::string_view kPairBuildSelfEEQuery = "contact.surface.pair_build.self_ee_query";
inline constexpr std::string_view kPairBuildSelfEEHashCandidates = "contact.surface.pair_build.self_ee.hash_candidates";
inline constexpr std::string_view kPairBuildSelfEEDistanceTests = "contact.surface.pair_build.self_ee.distance_tests";
inline constexpr std::string_view kPairBuildSelfEEAcceptedPairs = "contact.surface.pair_build.self_ee.accepted_pairs";
inline constexpr std::string_view kPairBuildExternal = "contact.surface.pair_build.external";
inline constexpr std::string_view kPairBuildExternalAABB = "contact.surface.pair_build.external_aabb";
inline constexpr std::string_view kPairBuildExternalPT = "contact.surface.pair_build.external_pt";
inline constexpr std::string_view kPairBuildExternalTP = "contact.surface.pair_build.external_tp";
inline constexpr std::string_view kPairBuildExternalEE = "contact.surface.pair_build.external_ee";
inline constexpr std::string_view kPairBuildExternalOverlappingObstacles = "contact.surface.pair_build.external.overlapping_obstacles";
inline constexpr std::string_view kPairBuildExternalPTHashCandidates = "contact.surface.pair_build.external_pt.hash_candidates";
inline constexpr std::string_view kPairBuildExternalPTDistanceTests = "contact.surface.pair_build.external_pt.distance_tests";
inline constexpr std::string_view kPairBuildExternalPTAcceptedPairs = "contact.surface.pair_build.external_pt.accepted_pairs";
inline constexpr std::string_view kPairBuildExternalTPHashCandidates = "contact.surface.pair_build.external_tp.hash_candidates";
inline constexpr std::string_view kPairBuildExternalTPDistanceTests = "contact.surface.pair_build.external_tp.distance_tests";
inline constexpr std::string_view kPairBuildExternalTPAcceptedPairs = "contact.surface.pair_build.external_tp.accepted_pairs";
inline constexpr std::string_view kPairBuildExternalEEHashCandidates = "contact.surface.pair_build.external_ee.hash_candidates";
inline constexpr std::string_view kPairBuildExternalEEDistanceTests = "contact.surface.pair_build.external_ee.distance_tests";
inline constexpr std::string_view kPairBuildExternalEEAcceptedPairs = "contact.surface.pair_build.external_ee.accepted_pairs";
inline constexpr std::string_view kMaxStepPT = "contact.surface.max_step_pt";
inline constexpr std::string_view kMaxStepEE = "contact.surface.max_step_ee";
inline constexpr std::string_view kMaxStepSelfPTHashCandidates = "contact.surface.max_step.self_pt.hash_candidates";
inline constexpr std::string_view kMaxStepSelfPTCCDTests = "contact.surface.max_step.self_pt.ccd_tests";
inline constexpr std::string_view kMaxStepSelfEEHashCandidates = "contact.surface.max_step.self_ee.hash_candidates";
inline constexpr std::string_view kMaxStepSelfEECCDTests = "contact.surface.max_step.self_ee.ccd_tests";
inline constexpr std::string_view kMaxStepExternalOverlappingObstacles = "contact.surface.max_step.external.overlapping_obstacles";
inline constexpr std::string_view kMaxStepExternalPTHashCandidates = "contact.surface.max_step.external_pt.hash_candidates";
inline constexpr std::string_view kMaxStepExternalPTCCDTests = "contact.surface.max_step.external_pt.ccd_tests";
inline constexpr std::string_view kMaxStepExternalTPHashCandidates = "contact.surface.max_step.external_tp.hash_candidates";
inline constexpr std::string_view kMaxStepExternalTPCCDTests = "contact.surface.max_step.external_tp.ccd_tests";
inline constexpr std::string_view kMaxStepExternalEEHashCandidates = "contact.surface.max_step.external_ee.hash_candidates";
inline constexpr std::string_view kMaxStepExternalEECCDTests = "contact.surface.max_step.external_ee.ccd_tests";
inline constexpr std::string_view kEnergy = "contact.surface.energy";
inline constexpr std::string_view kGradient = "contact.surface.gradient";
inline constexpr std::string_view kHessian = "contact.surface.hessian";
inline constexpr std::string_view kCombined = "contact.surface.combined";
inline constexpr std::string_view kBuildActiveSet = "contact.surface.build_active_set";
inline constexpr std::string_view kActiveSetEnergy = "contact.surface.active_set_energy";
inline constexpr std::string_view kActiveSetGradient = "contact.surface.active_set_gradient";
inline constexpr std::string_view kActiveSetHessian = "contact.surface.active_set_hessian";
inline constexpr std::string_view kActiveSetCombined = "contact.surface.active_set_combined";
inline constexpr std::string_view kActiveSetSelfCombined = "contact.surface.active_set_combined.self";
inline constexpr std::string_view kActiveSetSelfPTCombined = "contact.surface.active_set_combined.self_pt";
inline constexpr std::string_view kActiveSetSelfEECombined = "contact.surface.active_set_combined.self_ee";
inline constexpr std::string_view kActiveSetSelfPTPairCount = "contact.surface.active_set_combined.self_pt.pairs";
inline constexpr std::string_view kActiveSetSelfEEPairCount = "contact.surface.active_set_combined.self_ee.pairs";
inline constexpr std::string_view kActiveSetExternalCombined = "contact.surface.active_set_combined.external";
inline constexpr std::string_view kActiveSetExternalPTCombined = "contact.surface.active_set_combined.external_pt";
inline constexpr std::string_view kActiveSetExternalTPCombined = "contact.surface.active_set_combined.external_tp";
inline constexpr std::string_view kActiveSetExternalEECombined = "contact.surface.active_set_combined.external_ee";
inline constexpr std::string_view kActiveSetExternalPTPairCount = "contact.surface.active_set_combined.external_pt.pairs";
inline constexpr std::string_view kActiveSetExternalTPPairCount = "contact.surface.active_set_combined.external_tp.pairs";
inline constexpr std::string_view kActiveSetExternalEEPairCount = "contact.surface.active_set_combined.external_ee.pairs";
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
