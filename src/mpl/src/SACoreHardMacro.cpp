// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#include "SACoreHardMacro.h"

#include <cmath>
#include <vector>

#include "MplObserver.h"
#include "SimulatedAnnealingCore.h"
#include "clusterEngine.h"
#include "object.h"
#include "odb/db.h"
#include "util.h"
#include "utl/Logger.h"

namespace mpl {
using utl::MPL;

//////////////////////////////////////////////////////////////////
// Class SACoreHardMacro
// constructors
SACoreHardMacro::SACoreHardMacro(PhysicalHierarchy* tree,
                                 const Rect& outline,
                                 const std::vector<HardMacro>& macros,
                                 const SACoreWeights& core_weights,
                                 // probability of each action
                                 float pos_swap_prob,
                                 float neg_swap_prob,
                                 float double_swap_prob,
                                 float exchange_prob,
                                 float flip_prob,
                                 // Fast SA hyperparameter
                                 float init_prob,
                                 int max_num_step,
                                 int num_perturb_per_step,
                                 unsigned seed,
                                 MplObserver* graphics,
                                 utl::Logger* logger,
                                 odb::dbBlock* block)
    : SimulatedAnnealingCore<HardMacro>(tree,
                                        outline,
                                        macros,
                                        core_weights,
                                        pos_swap_prob,
                                        neg_swap_prob,
                                        double_swap_prob,
                                        exchange_prob,
                                        init_prob,
                                        max_num_step,
                                        num_perturb_per_step,
                                        seed,
                                        graphics,
                                        logger,
                                        block)
{
  flip_prob_ = flip_prob;
}

void SACoreHardMacro::run()
{
  if (graphics_) {
    graphics_->startSA();
  }

  fastSA();

  if (graphics_) {
    graphics_->endSA(calNormCost());
  }
}

float SACoreHardMacro::calNormCost() const
{
  float cost = 0.0;  // Initialize cost
  if (norm_area_penalty_ > 0.0) {
    cost += core_weights_.area * getAreaPenalty();
  }
  if (norm_outline_penalty_ > 0.0) {
    cost += core_weights_.outline * outline_penalty_ / norm_outline_penalty_;
  }
  if (norm_wirelength_ > 0.0) {
    cost += core_weights_.wirelength * wirelength_ / norm_wirelength_;
  }
  if (norm_guidance_penalty_ > 0.0) {
    cost += core_weights_.guidance * guidance_penalty_ / norm_guidance_penalty_;
  }
  if (norm_fence_penalty_ > 0.0) {
    cost += core_weights_.fence * fence_penalty_ / norm_fence_penalty_;
  }
  return cost;
}

void SACoreHardMacro::calPenalty()
{
  calOutlinePenalty();
  calWirelength();
  calGuidancePenalty();
  calFencePenalty();
  if (graphics_) {
    graphics_->setAreaPenalty(
        {"Area", core_weights_.area, getAreaPenalty(), norm_area_penalty_});
    graphics_->penaltyCalculated(calNormCost());
  }
}

void SACoreHardMacro::flipAllMacros()
{
  for (auto& macro_id : pos_seq_) {
    macros_[macro_id].flip(false);
  }
}

void SACoreHardMacro::perturb()
{
  if (macros_.empty()) {
    return;
  }

  // Keep back up
  pre_pos_seq_ = pos_seq_;
  pre_neg_seq_ = neg_seq_;
  pre_width_ = width_;
  pre_height_ = height_;
  pre_outline_penalty_ = outline_penalty_;
  pre_wirelength_ = wirelength_;
  pre_guidance_penalty_ = guidance_penalty_;
  pre_fence_penalty_ = fence_penalty_;

  // generate random number (0 - 1) to determine actions
  const float op = distribution_(generator_);
  const float action_prob_1 = pos_swap_prob_;
  const float action_prob_2 = action_prob_1 + neg_swap_prob_;
  const float action_prob_3 = action_prob_2 + double_swap_prob_;
  const float action_prob_4 = action_prob_3 + exchange_prob_;
  if (op <= action_prob_1) {
    action_id_ = 1;
    singleSeqSwap(true);  // Swap two macros in pos_seq_
  } else if (op <= action_prob_2) {
    action_id_ = 2;
    singleSeqSwap(false);  // Swap two macros in neg_seq_;
  } else if (op <= action_prob_3) {
    action_id_ = 3;
    doubleSeqSwap();  // Swap two macros in pos_seq_ and
                      // other two macros in neg_seq_
  } else if (op <= action_prob_4) {
    action_id_ = 4;
    exchangeMacros();  // exchange two macros in the sequence pair
  } else {
    action_id_ = 5;
    pre_macros_ = macros_;
    flipAllMacros();
  }

  // update the macro locations based on Sequence Pair
  packFloorplan();
  // Update all the penalties
  calPenalty();
}

void SACoreHardMacro::restore()
{
  if (macros_.empty()) {
    return;
  }

  // To reduce the runtime, here we do not call PackFloorplan
  // again. So when we need to generate the final floorplan out,
  // we need to call PackFloorplan again at the end of SA process
  if (action_id_ == 5) {
    macros_ = pre_macros_;
    // macros_[macro_id_] = pre_macros_[macro_id_];
  } else if (action_id_ == 1) {
    pos_seq_ = pre_pos_seq_;
  } else if (action_id_ == 2) {
    neg_seq_ = pre_neg_seq_;
  } else {
    pos_seq_ = pre_pos_seq_;
    neg_seq_ = pre_neg_seq_;
  }

  width_ = pre_width_;
  height_ = pre_height_;
  outline_penalty_ = pre_outline_penalty_;
  wirelength_ = pre_wirelength_;
  guidance_penalty_ = pre_guidance_penalty_;
  fence_penalty_ = pre_fence_penalty_;
}

void SACoreHardMacro::initialize()
{
  initSequencePair();

  std::vector<float> area_penalty_list;
  std::vector<float> outline_penalty_list;
  std::vector<float> wirelength_list;
  std::vector<float> guidance_penalty_list;
  std::vector<float> fence_penalty_list;
  std::vector<float> width_list;
  std::vector<float> height_list;

  // We don't want to stop in the normalization factor setup
  MplObserver* save_graphics = graphics_;
  graphics_ = nullptr;

  for (int i = 0; i < num_perturb_per_step_; i++) {
    perturb();
    // store current penalties
    width_list.push_back(width_);
    height_list.push_back(height_);
    area_penalty_list.push_back(width_ * height_);
    outline_penalty_list.push_back(outline_penalty_);
    wirelength_list.push_back(wirelength_);
    guidance_penalty_list.push_back(guidance_penalty_);
    fence_penalty_list.push_back(fence_penalty_);
  }
  graphics_ = save_graphics;

  norm_area_penalty_ = calAverage(area_penalty_list);
  norm_outline_penalty_ = calAverage(outline_penalty_list);
  norm_wirelength_ = calAverage(wirelength_list);
  norm_guidance_penalty_ = calAverage(guidance_penalty_list);
  norm_fence_penalty_ = calAverage(fence_penalty_list);

  if (norm_area_penalty_ <= 1e-4) {
    norm_area_penalty_ = 1.0;
  }

  if (norm_outline_penalty_ <= 1e-4) {
    norm_outline_penalty_ = 1.0;
  }

  if (norm_wirelength_ <= 1e-4) {
    norm_wirelength_ = 1.0;
  }

  if (norm_guidance_penalty_ <= 1e-4) {
    norm_guidance_penalty_ = 1.0;
  }

  if (norm_fence_penalty_ <= 1e-4) {
    norm_fence_penalty_ = 1.0;
  }

  // Calculate initial temperature
  std::vector<float> cost_list;
  for (int i = 0; i < outline_penalty_list.size(); i++) {
    width_ = width_list[i];
    height_ = height_list[i];
    outline_penalty_ = outline_penalty_list[i];
    wirelength_ = wirelength_list[i];
    guidance_penalty_ = guidance_penalty_list[i];
    fence_penalty_ = fence_penalty_list[i];
    cost_list.push_back(calNormCost());
  }
  float delta_cost = 0.0;
  for (int i = 1; i < cost_list.size(); i++) {
    delta_cost += std::abs(cost_list[i] - cost_list[i - 1]);
  }
  if (cost_list.size() > 1 && delta_cost > 0.0) {
    init_temperature_
        = (-1.0) * (delta_cost / (cost_list.size() - 1)) / std::log(init_prob_);
  } else {
    init_temperature_ = 1.0;
  }
}

void SACoreHardMacro::setWeights(const SACoreWeights& weights)
{
  core_weights_ = weights;
}

void SACoreHardMacro::printResults() const
{
  reportCoreWeights();
  reportTotalCost();
  if (logger_->debugCheck(MPL, "hierarchical_macro_placement", 2)) {
    reportLocations();
  }
}

}  // namespace mpl
