// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025, The OpenROAD Authors

#include "ResAwareMove.hh"

#include "grt/GlobalRouter.h"
#include "rsz/Resizer.hh"
#include "sta/Path.hh"

namespace rsz {

bool ResAwareMove::doMove(const sta::Path* drvr_path, float setup_slack_margin)
{
  if (!resizer_->global_router_ || !resizer_->global_router_->haveRoutes()) {
    return false;
  }

  sta::Pin* drvr_pin = drvr_path->pin(this);
  if (!drvr_pin) {
    return false;
  }

  sta::Instance* drvr_inst = network_->instance(drvr_pin);
  if (resizer_->dontTouch(drvr_inst)) {
    return false;
  }

  if (hasMoves(drvr_inst)) {
    return false;
  }

  sta::Net* net = network_->net(drvr_pin);
  if (!net) {
    return false;
  }

  odb::dbNet* db_net = db_network_->staToDb(net);
  if (!db_net) {
    return false;
  }

  // Check if we should even try to reroute
  if (resizer_->global_router_->isNetResAware(db_net)) {
    return false;
  }

  // Check net resistance before accepting the move
  // If there is no improvement, skip it
  float resistance = resizer_->global_router_->getFRNetResistance(db_net);

  resizer_->global_router_->setResistanceAware(true);
  resizer_->global_router_->addDirtyNet(db_net);
  resizer_->global_router_->setNetIsResAware(db_net, true);

  grt::IncrementalGRoute incr_groute(resizer_->global_router_,
                                     resizer_->getDbBlock());
  incr_groute.updateRoutes();

  float new_resistance = resizer_->global_router_->getFRNetResistance(db_net);
  if (resistance <= new_resistance) {
    return false;
  }

  estimate_parasitics_->parasiticsInvalid(db_net);

  countMove(drvr_inst);
  return true;
}

}  // namespace rsz
