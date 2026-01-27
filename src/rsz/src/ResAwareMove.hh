// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025, The OpenROAD Authors

#pragma once

#include "BaseMove.hh"

namespace rsz {

class ResAwareMove : public BaseMove
{
 public:
  using BaseMove::BaseMove;

  bool doMove(const sta::Path* drvr_path,
              int drvr_index,
              sta::Slack drvr_slack,
              sta::PathExpanded* expanded,
              float setup_slack_margin) override;

  const char* name() override { return "ResAwareMove"; }
};

}  // namespace rsz
