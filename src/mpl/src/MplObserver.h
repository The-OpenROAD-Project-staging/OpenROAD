///////////////////////////////////////////////////////////////////////////
//
// BSD 3-Clause License
//
// Copyright (c) 2023, Google LLC
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <optional>
#include <vector>

#include "clusterEngine.h"
#include "object.h"
#include "odb/geom.h"
#include "util.h"
#include "utl/Logger.h"

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
