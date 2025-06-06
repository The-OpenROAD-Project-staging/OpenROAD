// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

#pragma once

#include <map>
#include <vector>

#include "clusterEngine.h"
#include "object.h"
#include "odb/geom.h"
#include "util.h"

namespace mpl {

class Cluster;

class MplObserver
{
 public:
  MplObserver() = default;
  virtual ~MplObserver() = default;

  virtual void startCoarse() {}
  virtual void startFine() {}

  virtual void startSA() {}
  virtual void saStep(const std::vector<SoftMacro>& macros) {}
  virtual void saStep(const std::vector<HardMacro>& macros) {}
  virtual void endSA(float norm_cost) {}
  virtual void drawResult() {}

  virtual void finishedClustering(PhysicalHierarchy* tree) {}

  virtual void setMaxLevel(int max_level) {}
  virtual void setMacroBlockages(const std::vector<mpl::Rect>& macro_blockages)
  {
  }
  virtual void setPlacementBlockages(
      const std::vector<mpl::Rect>& placement_blockages)
  {
  }
  virtual void setBundledNets(const std::vector<BundledNet>& bundled_nets) {}
  virtual void setShowBundledNets(bool show_bundled_nets) {}
  virtual void setShowClustersIds(bool show_clusters_ids) {}
  virtual void setSkipSteps(bool skip_steps) {}
  virtual void doNotSkip() {}
  virtual void setOnlyFinalResult(bool skip_to_end) {}
  virtual void setTargetClusterId(int target_cluster_id) {}
  virtual void setCurrentCluster(Cluster* current_cluster) {}

  virtual void setOutline(const odb::Rect& outline) {}
  virtual void setGuides(const std::map<int, Rect>& guides) {}
  virtual void setFences(const std::map<int, Rect>& fences) {}
  virtual void setIOConstraintsMap(
      const ClusterToBoundaryRegionMap& io_cluster_to_constraint)
  {
  }
  virtual void setBlockedRegionsForPins(
      const std::vector<odb::Rect>& blocked_regions_for_pins)
  {
  }
  virtual void setAvailableRegionsForUnconstrainedPins(
      const BoundaryRegionList& regions)
  {
  }

  virtual void setAreaPenalty(const PenaltyData& penalty) {}
  virtual void setBoundaryPenalty(const PenaltyData& penalty) {}
  virtual void setFencePenalty(const PenaltyData& penalty) {}
  virtual void setGuidancePenalty(const PenaltyData& penalty) {}
  virtual void setMacroBlockagePenalty(const PenaltyData& penalty) {}
  virtual void setNotchPenalty(const PenaltyData& penalty) {}
  virtual void setOutlinePenalty(const PenaltyData& penalty) {}
  virtual void setWirelengthPenalty(const PenaltyData& penalty) {}
  virtual void penaltyCalculated(float norm_cost) {}

  virtual void eraseDrawing() {}
};

}  // namespace mpl
