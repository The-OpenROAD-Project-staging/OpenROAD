// Copyright 2023 Google LLC
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file or at
// https://developers.google.com/open-source/licenses/bsd

#include <tcl.h>
#include <unistd.h>

#include <array>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <string>

#include "db_sta/MakeDbSta.hh"
#include "db_sta/dbNetwork.hh"
#include "db_sta/dbSta.hh"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "odb/lefin.h"
#include "sta/Corner.hh"
#include "sta/FuncExpr.hh"
#include "sta/Graph.hh"
#include "sta/Liberty.hh"
#include "sta/PathAnalysisPt.hh"
#include "sta/Search.hh"
#include "sta/Sta.hh"
#include "sta/Units.hh"
#include "utl/Logger.h"
#include "utl/deleter.h"

namespace odb {

std::once_flag init_sta_flag;

/*
  Extract the hierarchical information in human readable format.
  Shows the dbNet and dbModNet view of the database.
*/

void DbStrDebugHierarchy(dbBlock* block, std::stringstream& str_db)
{
  char tmp_str[10000];
  sprintf(tmp_str,
          "Debug: Data base tables for block at %s:\n",
          block->getName().c_str());
  str_db << tmp_str;

  sprintf(tmp_str, "Db nets (The Flat db view)\n");
  str_db << tmp_str;

  for (auto dbnet : block->getNets()) {
    sprintf(tmp_str,
            "dbnet %s (id %u)\n",
            dbnet->getName().c_str(),
            dbnet->getId());
    str_db << tmp_str;

    for (auto db_iterm : dbnet->getITerms()) {
      sprintf(tmp_str,
              "\t-> Db Iterm %u (%s)\n",
              db_iterm->getId(),
              db_iterm->getName().c_str());
      str_db << tmp_str;
    }
    for (auto db_bterm : dbnet->getBTerms()) {
      sprintf(tmp_str, "\t-> Db Bterm %u\n", db_bterm->getId());
      str_db << tmp_str;
    }
  }

  sprintf(tmp_str, "Block ports\n");
  // got through the ports and their owner
  sprintf(tmp_str, "\t\tBTerm Ports +++\n");
  str_db << tmp_str;
  dbSet<dbBTerm> block_bterms = block->getBTerms();
  for (auto bt : block_bterms) {
    sprintf(tmp_str,
            "\t\tBterm (%u) %s Net %s (%u)  Mod Net %s (%u) \n",
            bt->getId(),
            bt->getName().c_str(),
            bt->getNet() ? bt->getNet()->getName().c_str() : "",
            bt->getNet() ? bt->getNet()->getId() : 0,
            bt->getModNet() ? bt->getModNet()->getName() : "",
            bt->getModNet() ? bt->getModNet()->getId() : 0);
    str_db << tmp_str;
  }
  sprintf(tmp_str, "\t\tBTerm Ports ---\n");
  str_db << tmp_str;

  sprintf(tmp_str, "The hierarchical db view:\n");
  str_db << tmp_str;
  dbSet<dbModule> block_modules = block->getModules();
  sprintf(tmp_str, "Content size %u modules\n", block_modules.size());
  str_db << tmp_str;
  for (auto mi : block_modules) {
    dbModule* cur_obj = mi;
    if (cur_obj == block->getTopModule()) {
      sprintf(tmp_str, "Top Module\n");
      str_db << tmp_str;
    }
    sprintf(tmp_str,
            "\tModule %s %s\n",
            (cur_obj == block->getTopModule()) ? "(Top Module)" : "",
            ((dbModule*) cur_obj)->getName());
    str_db << tmp_str;
    // in case of top level, care as the bterms double up as pins
    if (cur_obj == block->getTopModule()) {
      for (auto bterm : block->getBTerms()) {
        sprintf(tmp_str,
                "Top B-term %s dbNet %s (%d) modNet %s (%d)\n",
                bterm->getName().c_str(),
                bterm->getNet() ? bterm->getNet()->getName().c_str() : "",
                bterm->getNet() ? bterm->getNet()->getId() : -1,
                bterm->getModNet() ? bterm->getModNet()->getName() : "",
                bterm->getModNet() ? bterm->getModNet()->getId() : -1);
        str_db << tmp_str;
      }
    }
    // got through the module ports and their owner
    sprintf(tmp_str, "\t\tModBTerm Ports +++\n");
    str_db << tmp_str;

    dbSet<dbModBTerm> module_ports = cur_obj->getModBTerms();
    for (auto module_port : module_ports) {
      sprintf(
          tmp_str,
          "\t\tPort %s Net %s (%d)\n",
          module_port->getName(),
          (module_port->getModNet()) ? (module_port->getModNet()->getName())
                                     : "No-modnet",
          (module_port->getModNet()) ? module_port->getModNet()->getId() : -1);
      str_db << tmp_str;
      sprintf(tmp_str,
              "\t\tPort parent %s\n\n",
              module_port->getParent()->getName());
      str_db << tmp_str;
    }
    sprintf(tmp_str, "\t\tModBTermPorts ---\n");
    str_db << tmp_str;

    sprintf(tmp_str, "\t\tModule instances +++\n");
    str_db << tmp_str;
    dbSet<dbModInst> module_instances = mi->getModInsts();
    for (auto module_inst : module_instances) {
      sprintf(tmp_str, "\t\tMod inst %s ", module_inst->getName());
      str_db << tmp_str;
      dbModule* master = module_inst->getMaster();
      sprintf(
          tmp_str, "\t\tMaster %s\n\n", module_inst->getMaster()->getName());
      str_db << tmp_str;
      dbBlock* owner = master->getOwner();
      if (owner != block) {
        sprintf(tmp_str, "\t\t\tMaster owner in wrong block\n");
        str_db << tmp_str;
      }
      sprintf(tmp_str, "\t\tConnections\n");
      str_db << tmp_str;
      for (dbModITerm* miterm_pin : module_inst->getModITerms()) {
        sprintf(tmp_str,
                "\t\t\tModIterm : %s (%u) Mod Net %s (%u) \n",
                miterm_pin->getName(),
                miterm_pin->getId(),
                miterm_pin->getModNet() ? (miterm_pin->getModNet()->getName())
                                        : "No-net",
                miterm_pin->getModNet() ? miterm_pin->getModNet()->getId() : 0);
        str_db << tmp_str;
      }
    }
    sprintf(tmp_str, "\t\tModule instances ---\n");
    str_db << tmp_str;
    sprintf(tmp_str, "\t\tDb instances +++\n");
    str_db << tmp_str;
    for (dbInst* db_inst : cur_obj->getInsts()) {
      sprintf(tmp_str, "\t\tdb inst %s\n", db_inst->getName().c_str());
      str_db << tmp_str;
      sprintf(tmp_str, "\t\tdb iterms:\n");
      str_db << tmp_str;
      dbSet<dbITerm> iterms = db_inst->getITerms();
      for (auto iterm : iterms) {
        dbMTerm* mterm = iterm->getMTerm();
        sprintf(
            tmp_str,
            "\t\t\t\t iterm: %s (%u) Net: %s Mod net : %s (%d)\n",
            mterm->getName().c_str(),
            iterm->getId(),
            iterm->getNet() ? iterm->getNet()->getName().c_str() : "unk-dbnet",
            iterm->getModNet() ? iterm->getModNet()->getName() : "unk-modnet",
            iterm->getModNet() ? iterm->getModNet()->getId() : -1);
        str_db << tmp_str;
      }
    }
    sprintf(tmp_str, "\t\tDb instances ---\n");
    str_db << tmp_str;
    sprintf(tmp_str, "\tModule nets (modnets) +++ \n");
    str_db << tmp_str;
    sprintf(tmp_str,
            "\t# mod nets %u in %s \n",
            cur_obj->getModNets().size(),
            cur_obj->getName());

    str_db << tmp_str;
    dbSet<dbModNet> mod_nets = cur_obj->getModNets();
    for (auto mod_net : mod_nets) {
      sprintf(
          tmp_str, "\t\tNet: %s (%u)\n", mod_net->getName(), mod_net->getId());
      str_db << tmp_str;
      sprintf(tmp_str,
              "\t\tConnections -> modIterms/modbterms/bterms/iterms:\n");
      str_db << tmp_str;
      sprintf(
          tmp_str, "\t\t -> %u moditerms\n", mod_net->getModITerms().size());
      str_db << tmp_str;
      for (dbModITerm* modi_term : mod_net->getModITerms()) {
        sprintf(tmp_str, "\t\t\t%s\n", modi_term->getName());
        str_db << tmp_str;
      }
      sprintf(
          tmp_str, "\t\t -> %u modbterms\n", mod_net->getModBTerms().size());
      str_db << tmp_str;
      for (dbModBTerm* modb_term : mod_net->getModBTerms()) {
        sprintf(tmp_str, "\t\t\t%s\n", modb_term->getName());
        str_db << tmp_str;
      }
      sprintf(tmp_str, "\t\t -> %u iterms\n", mod_net->getITerms().size());
      str_db << tmp_str;
      for (dbITerm* db_iterm : mod_net->getITerms()) {
        sprintf(tmp_str, "\t\t\t%s\n", db_iterm->getName().c_str());
        str_db << tmp_str;
      }
      sprintf(tmp_str, "\t\t -> %u bterms\n", mod_net->getBTerms().size());
      str_db << tmp_str;
      for (dbBTerm* db_bterm : mod_net->getBTerms()) {
        sprintf(tmp_str, "\t\t\t%s\n", db_bterm->getName().c_str());
        str_db << tmp_str;
      }
    }
  }
}

/*
   Test fixture for hierarchical connect test
   ------------------------------------------

   Two modules, 1 with 1 inverter, second with 4 inverters.

   inv1 is instantiated in the top level
   inv4 is instantiated three levels down:
   level0 -- instantiation in root
   level1
   level2 -- contains 4 inverters

   Objective of this test is to stress the hierarchical connection code.

   root {inputs: i0; outputs: o0, o1, o2,o3}
       --> inv1 (module containing 1 inverter) this drives inv1_level0_inst
       --> inv4_level0_inst  {wrapper input: i0, outputs o0,o1,o2,o3}
           --> inv4_level1_inst {wrapper input: i0, outputs o0,o1,o2,o3}
               --> inv4_level2_inst {contains 4 inverters. Has input i0, outputs
   o0,01,o2,o3}
*/

class TestHconn : public ::testing::Test
{
 protected:
  void SetUp() override
  {
    // this will be so much easier with read_def
    db_ = utl::UniquePtrWithDeleter<odb::dbDatabase>(odb::dbDatabase::create(),
                                                     &odb::dbDatabase::destroy);
    std::call_once(init_sta_flag, []() { sta::initSta(); });
    sta_ = std::unique_ptr<sta::dbSta>(sta::makeDbSta());
    sta_->initVars(Tcl_CreateInterp(), db_.get(), &logger_);
    auto path = std::filesystem::canonical("./Nangate45/Nangate45_typ.lib");
    library_ = sta_->readLiberty(path.string().c_str(),
                                 sta_->findCorner("default"),
                                 sta::MinMaxAll::all(),
                                 /*infer_latches=*/false);
    odb::lefin lef_parser(
        db_.get(), &logger_, /*ignore_non_routing_layers*/ false);
    const char* lib_name = "Nangate45.lef";
    lib_ = lef_parser.createTechAndLib(
        "tech", lib_name, "./Nangate45/Nangate45.lef");

    sta_->postReadLef(/*tech=*/nullptr, lib_);

    sta::Units* units = library_->units();
    power_unit_ = units->powerUnit();
    db_network_ = sta_->getDbNetwork();
    // turn on hierarchy
    db_network_->setHierarchy();
    db_->setLogger(&logger_);

    // create a chain consisting of 4 buffers
    odb::dbChip* chip = odb::dbChip::create(db_.get());
    block_ = odb::dbBlock::create(chip, "top");
    db_network_->setBlock(block_);
    block_->setDieArea(odb::Rect(0, 0, 1000, 1000));
    // register proper callbacks for timer like read_def
    sta_->postReadDef(block_);

    root_mod_ = dbModule::create(block_, "root_mod");

    // The bterms are created below during wiring
    // Note a bterm without a parent is a root bterm.

    // Make the inv1 module (contains 1 inverter).
    inv1_mod_master_ = dbModule::create(block_, "inv1_master_mod");
    inv1_mod_i0_port_ = dbModBTerm::create(inv1_mod_master_, "i0");
    inv1_mod_o0_port_ = dbModBTerm::create(inv1_mod_master_, "o0");

    inv1_mod_inst_
        = dbModInst::create(root_mod_, inv1_mod_master_, "inv1_mod_inst");

    inv1_mod_inst_i0_miterm_
        = dbModITerm::create(inv1_mod_inst_, "i0", inv1_mod_i0_port_);
    inv1_mod_inst_o0_miterm_
        = dbModITerm::create(inv1_mod_inst_, "o0", inv1_mod_o0_port_);
    // correlate the iterms and bterms
    inv1_mod_inst_i0_miterm_->setChildModBTerm(inv1_mod_i0_port_);
    inv1_mod_i0_port_->setParentModITerm(inv1_mod_inst_i0_miterm_);

    inv1_mod_inst_o0_miterm_->setChildModBTerm(inv1_mod_o0_port_);
    inv1_mod_o0_port_->setParentModITerm(inv1_mod_inst_o0_miterm_);

    inv1_1_ = dbInst::create(block_,
                             lib_->findMaster("INV_X1"),
                             "inv1_mod_inst/inst1",
                             false,
                             inv1_mod_master_);

    inv1_1_inst_ip0_ = inv1_1_->findITerm("A");
    inv1_1_inst_op0_ = inv1_1_->findITerm("ZN");

    inv4_mod_level0_master_ = dbModule::create(block_, "inv4_master_level0");
    inv4_mod_level0_master_i0_port_
        = dbModBTerm::create(inv4_mod_level0_master_, "i0");
    inv4_mod_level0_master_o0_port_
        = dbModBTerm::create(inv4_mod_level0_master_, "o0");
    inv4_mod_level0_master_o1_port_
        = dbModBTerm::create(inv4_mod_level0_master_, "o1");
    inv4_mod_level0_master_o2_port_
        = dbModBTerm::create(inv4_mod_level0_master_, "o2");
    inv4_mod_level0_master_o3_port_
        = dbModBTerm::create(inv4_mod_level0_master_, "o3");

    inv4_mod_level1_master_ = dbModule::create(block_, "inv4_master_level1");
    inv4_mod_level1_master_i0_port_
        = dbModBTerm::create(inv4_mod_level1_master_, "i0");
    inv4_mod_level1_master_o0_port_
        = dbModBTerm::create(inv4_mod_level1_master_, "o0");
    inv4_mod_level1_master_o1_port_
        = dbModBTerm::create(inv4_mod_level1_master_, "o1");
    inv4_mod_level1_master_o2_port_
        = dbModBTerm::create(inv4_mod_level1_master_, "o2");
    inv4_mod_level1_master_o3_port_
        = dbModBTerm::create(inv4_mod_level1_master_, "o3");

    inv4_mod_level2_master_ = dbModule::create(block_, "inv4_master_level2");
    inv4_mod_level2_master_i0_port_
        = dbModBTerm::create(inv4_mod_level2_master_, "i0");
    inv4_mod_level2_master_o0_port_
        = dbModBTerm::create(inv4_mod_level2_master_, "o0");
    inv4_mod_level2_master_o1_port_
        = dbModBTerm::create(inv4_mod_level2_master_, "o1");
    inv4_mod_level2_master_o2_port_
        = dbModBTerm::create(inv4_mod_level2_master_, "o2");
    inv4_mod_level2_master_o3_port_
        = dbModBTerm::create(inv4_mod_level2_master_, "o3");

    // During modinst creation we set the parent.
    inv4_mod_level0_inst_ = dbModInst::create(root_mod_,  // parent
                                              inv4_mod_level0_master_,
                                              "inv4_mod_level0_inst");
    inv4_mod_level0_inst_i0_miterm_
        = dbModITerm::create(inv4_mod_level0_inst_, "i0");
    inv4_mod_level0_inst_i0_miterm_->setChildModBTerm(
        inv4_mod_level0_master_i0_port_);
    inv4_mod_level0_master_i0_port_->setParentModITerm(
        inv4_mod_level0_inst_i0_miterm_);

    inv4_mod_level0_inst_o0_miterm_
        = dbModITerm::create(inv4_mod_level0_inst_, "o0");
    inv4_mod_level0_inst_o0_miterm_->setChildModBTerm(
        inv4_mod_level0_master_o0_port_);
    inv4_mod_level0_master_o0_port_->setParentModITerm(
        inv4_mod_level0_inst_o0_miterm_);

    inv4_mod_level0_inst_o1_miterm_
        = dbModITerm::create(inv4_mod_level0_inst_, "o1");
    inv4_mod_level0_inst_o1_miterm_->setChildModBTerm(
        inv4_mod_level0_master_o1_port_);
    inv4_mod_level0_master_o1_port_->setParentModITerm(
        inv4_mod_level0_inst_o1_miterm_);

    inv4_mod_level0_inst_o2_miterm_
        = dbModITerm::create(inv4_mod_level0_inst_, "o2");
    inv4_mod_level0_inst_o2_miterm_->setChildModBTerm(
        inv4_mod_level0_master_o2_port_);
    inv4_mod_level0_master_o2_port_->setParentModITerm(
        inv4_mod_level0_inst_o2_miterm_);

    inv4_mod_level0_inst_o3_miterm_
        = dbModITerm::create(inv4_mod_level0_inst_, "o3");
    inv4_mod_level0_inst_o3_miterm_->setChildModBTerm(
        inv4_mod_level0_master_o3_port_);
    inv4_mod_level0_master_o3_port_->setParentModITerm(
        inv4_mod_level0_inst_o3_miterm_);

    inv4_mod_level1_inst_
        = dbModInst::create(inv4_mod_level0_master_,  // parent
                            inv4_mod_level1_master_,
                            "inv4_mod_level1_inst");

    inv4_mod_level1_inst_i0_miterm_
        = dbModITerm::create(inv4_mod_level1_inst_, "i0");
    inv4_mod_level1_inst_i0_miterm_->setChildModBTerm(
        inv4_mod_level1_master_i0_port_);
    inv4_mod_level1_master_i0_port_->setParentModITerm(
        inv4_mod_level1_inst_i0_miterm_);

    inv4_mod_level1_inst_o0_miterm_
        = dbModITerm::create(inv4_mod_level1_inst_, "o0");
    inv4_mod_level1_inst_o0_miterm_->setChildModBTerm(
        inv4_mod_level1_master_o0_port_);
    inv4_mod_level1_master_o0_port_->setParentModITerm(
        inv4_mod_level1_inst_o0_miterm_);

    inv4_mod_level1_inst_o1_miterm_
        = dbModITerm::create(inv4_mod_level1_inst_, "o1");
    inv4_mod_level1_inst_o1_miterm_->setChildModBTerm(
        inv4_mod_level1_master_o1_port_);
    inv4_mod_level1_master_o1_port_->setParentModITerm(
        inv4_mod_level1_inst_o1_miterm_);

    inv4_mod_level1_inst_o2_miterm_
        = dbModITerm::create(inv4_mod_level1_inst_, "o2");
    inv4_mod_level1_inst_o2_miterm_->setChildModBTerm(
        inv4_mod_level1_master_o2_port_);
    inv4_mod_level1_master_o2_port_->setParentModITerm(
        inv4_mod_level1_inst_o2_miterm_);

    inv4_mod_level1_inst_o3_miterm_
        = dbModITerm::create(inv4_mod_level1_inst_, "o3");
    inv4_mod_level1_inst_o3_miterm_->setChildModBTerm(
        inv4_mod_level1_master_o3_port_);
    inv4_mod_level1_master_o3_port_->setParentModITerm(
        inv4_mod_level1_inst_o3_miterm_);

    inv4_mod_level2_inst_
        = dbModInst::create(inv4_mod_level1_master_,  // parent
                            inv4_mod_level2_master_,
                            "inv4_mod_level2_inst");

    inv4_mod_level2_inst_i0_miterm_
        = dbModITerm::create(inv4_mod_level2_inst_, "i0");
    inv4_mod_level2_inst_i0_miterm_->setChildModBTerm(
        inv4_mod_level2_master_i0_port_);
    inv4_mod_level2_master_i0_port_->setParentModITerm(
        inv4_mod_level2_inst_i0_miterm_);

    inv4_mod_level2_inst_o0_miterm_
        = dbModITerm::create(inv4_mod_level2_inst_, "o0");
    inv4_mod_level2_inst_o0_miterm_->setChildModBTerm(
        inv4_mod_level2_master_o0_port_);
    inv4_mod_level2_master_o0_port_->setParentModITerm(
        inv4_mod_level2_inst_o0_miterm_);

    inv4_mod_level2_inst_o1_miterm_
        = dbModITerm::create(inv4_mod_level2_inst_, "o1");
    inv4_mod_level2_inst_o1_miterm_->setChildModBTerm(
        inv4_mod_level2_master_o1_port_);
    inv4_mod_level2_master_o1_port_->setParentModITerm(
        inv4_mod_level2_inst_o1_miterm_);

    inv4_mod_level2_inst_o2_miterm_
        = dbModITerm::create(inv4_mod_level2_inst_, "o2");
    inv4_mod_level2_inst_o2_miterm_->setChildModBTerm(
        inv4_mod_level2_master_o2_port_);
    inv4_mod_level2_master_o2_port_->setParentModITerm(
        inv4_mod_level2_inst_o2_miterm_);

    inv4_mod_level2_inst_o3_miterm_
        = dbModITerm::create(inv4_mod_level2_inst_, "o3");
    inv4_mod_level2_inst_o3_miterm_->setChildModBTerm(
        inv4_mod_level2_master_o3_port_);
    inv4_mod_level2_master_o3_port_->setParentModITerm(
        inv4_mod_level2_inst_o3_miterm_);

    //
    // create the low level inverter instances, for now uniquely name them in
    // the scope of a block.

    inv4_1_ = dbInst::create(
        block_,
        lib_->findMaster("INV_X1"),
        "inv4_mod_level0_inst/inv4_mod_level1_inst/inv4_mod_level2_inst/inst1");

    // get the iterm off the instance. Just give the terminal name
    // offset used to find iterm.
    inv4_1_ip_ = inv4_1_->findITerm("A");
    inv4_1_op_ = inv4_1_->findITerm("ZN");

    inv4_2_ = dbInst::create(
        block_,
        lib_->findMaster("INV_X1"),
        "inv4_mod_level0_inst/inv4_mod_level1_inst/inv4_mod_level2_inst/inst2");
    inv4_2_ip_ = inv4_2_->findITerm("A");
    inv4_2_op_ = inv4_2_->findITerm("ZN");

    inv4_3_ = dbInst::create(
        block_,
        lib_->findMaster("INV_X1"),
        "inv4_mod_level0_inst/inv4_mod_level1_inst/inv4_mod_level2_inst/inst3");
    inv4_3_ip_ = inv4_3_->findITerm("A");
    inv4_3_op_ = inv4_3_->findITerm("ZN");

    inv4_4_ = dbInst::create(
        block_,
        lib_->findMaster("INV_X1"),
        "inv4_mod_level0_inst/inv4_mod_level1_inst/inv4_mod_level2_inst/inst4");
    inv4_4_ip_ = inv4_4_->findITerm("A");
    inv4_4_op_ = inv4_4_->findITerm("ZN");

    //
    // inv4_mod_level2_master is lowest level in hierarchy
    //
    inv4_mod_level2_master_->addInst(inv4_1_);
    inv4_mod_level2_master_->addInst(inv4_2_);
    inv4_mod_level2_master_->addInst(inv4_3_);
    inv4_mod_level2_master_->addInst(inv4_4_);

    inv4_mod_level2_inst_i0_mnet_ = dbModNet::create(
        inv4_mod_level2_master_, "inv4_mod_level2_inst_i0_mnet");
    inv4_mod_level2_inst_o0_mnet_ = dbModNet::create(
        inv4_mod_level2_master_, "inv4_mod_level2_inst_o0_mnet");
    inv4_mod_level2_inst_o1_mnet_ = dbModNet::create(
        inv4_mod_level2_master_, "inv4_mod_level2_inst_o1_mnet");
    inv4_mod_level2_inst_o2_mnet_ = dbModNet::create(
        inv4_mod_level2_master_, "inv4_mod_level2_inst_o2_mnet");
    inv4_mod_level2_inst_o3_mnet_ = dbModNet::create(
        inv4_mod_level2_master_, "inv4_mod_level2_inst_o3_mnet");
    inv4_1_ip_->connect(inv4_mod_level2_inst_i0_mnet_);
    inv4_2_ip_->connect(inv4_mod_level2_inst_i0_mnet_);
    inv4_3_ip_->connect(inv4_mod_level2_inst_i0_mnet_);
    inv4_4_ip_->connect(inv4_mod_level2_inst_i0_mnet_);

    inv4_1_op_->connect(inv4_mod_level2_inst_o0_mnet_);
    inv4_2_op_->connect(inv4_mod_level2_inst_o1_mnet_);
    inv4_3_op_->connect(inv4_mod_level2_inst_o2_mnet_);
    inv4_4_op_->connect(inv4_mod_level2_inst_o3_mnet_);

    inv4_mod_level2_master_i0_port_->connect(inv4_mod_level2_inst_i0_mnet_);
    inv4_mod_level2_master_o0_port_->connect(inv4_mod_level2_inst_o0_mnet_);
    inv4_mod_level2_master_o1_port_->connect(inv4_mod_level2_inst_o1_mnet_);
    inv4_mod_level2_master_o2_port_->connect(inv4_mod_level2_inst_o2_mnet_);
    inv4_mod_level2_master_o3_port_->connect(inv4_mod_level2_inst_o3_mnet_);

    // First make the flat view connectivity
    ip0_net_ = dbNet::create(block_, "ip0_flat_net", false);
    inv_op_net_ = dbNet::create(block_, "inv_op_flat_net", false);
    op0_net_ = dbNet::create(block_, "op0_flat_net", false);
    op1_net_ = dbNet::create(block_, "op1_flat_net", false);
    op2_net_ = dbNet::create(block_, "op2_flat_net", false);
    op3_net_ = dbNet::create(block_, "op3_flat_net", false);

    //
    // connections to the primary (root) ports
    // Note: a bterm without a parent is a root port.
    //
    ip0_bterm_ = dbBTerm::create(ip0_net_, "ip0");
    op0_bterm_ = dbBTerm::create(op0_net_, "op0");
    op1_bterm_ = dbBTerm::create(op1_net_, "op1");
    op2_bterm_ = dbBTerm::create(op2_net_, "op2");
    op3_bterm_ = dbBTerm::create(op3_net_, "op3");

    op0_bterm_->connect(inv4_mod_level2_inst_o0_mnet_);
    op1_bterm_->connect(inv4_mod_level2_inst_o1_mnet_);
    op2_bterm_->connect(inv4_mod_level2_inst_o2_mnet_);
    op3_bterm_->connect(inv4_mod_level2_inst_o3_mnet_);

    // flat connections:
    // inverter 1 in module inv1_1
    inv1_1_inst_ip0_->connect(ip0_net_);
    inv1_1_inst_op0_->connect(inv_op_net_);

    // now the 4 inverters in inv4_1
    inv4_1_ip_->connect(inv_op_net_);
    inv4_2_ip_->connect(inv_op_net_);
    inv4_3_ip_->connect(inv_op_net_);
    inv4_4_ip_->connect(inv_op_net_);

    // now core to external connections
    inv4_1_op_->connect(op0_net_);
    inv4_2_op_->connect(op1_net_);
    inv4_3_op_->connect(op2_net_);
    inv4_4_op_->connect(op3_net_);

    // std::stringstream str_str;
    //    DbStrDebugHierarchy(block_, str_str);
    //    printf("The Flat design created %s\n", str_str.str().c_str());

    // Now build the hierarchical "overlay"
    // What we are doing here is adding the modnets which hook up

    // wire in hierarchy for inv1

    // Contents of inv1 (modnets from ports to internals).
    inv1_mod_i0_modnet_
        = dbModNet::create(inv1_mod_master_, "inv1_mod_i0_mnet");
    inv1_mod_o0_modnet_
        = dbModNet::create(inv1_mod_master_, "inv1_mod_o0_mnet");
    // connection from port to net in module
    inv1_mod_i0_port_->connect(inv1_mod_i0_modnet_);
    // connection to inverter ip
    inv1_1_inst_ip0_->connect(inv1_mod_i0_modnet_);
    // from inverter op to port
    inv1_1_inst_op0_->connect(inv1_mod_o0_modnet_);
    inv1_mod_o0_port_->connect(inv1_mod_o0_modnet_);

    // Instantiation of inv1 connections to top level (root bterms to modnets,
    // modnets to moditerms on inv1).
    root_inv1_i0_mnet_ = dbModNet::create(root_mod_, "inv1_inst_i0_mnet");
    root_inv1_o0_mnet_ = dbModNet::create(root_mod_, "inv1_inst_o0_mnet");
    inv1_mod_inst_i0_miterm_->connect(root_inv1_i0_mnet_);
    inv1_mod_inst_o0_miterm_->connect(root_inv1_o0_mnet_);

    // top level connections for inv1.
    // root input to instance, via a modnet.
    ip0_bterm_->connect(root_inv1_i0_mnet_);
    // note the output from out inverter module is inv1_inst_o0_mnet
    // which now needs to be hooked up the inv4 hierarchy.

    // The inv4 hierarchy connections
    // inv1 -> inv4 connection
    inv4_mod_level0_inst_i0_miterm_->connect(root_inv1_o0_mnet_);

    // level 0 is top of hierarchy, in root
    // the level 0 instance -> root connections
    root_inv4_o0_mnet_ = dbModNet::create(root_mod_, "inv4_inst_o0_mnet");
    inv4_mod_level0_inst_o0_miterm_->connect(root_inv4_o0_mnet_);
    op0_bterm_->connect(root_inv4_o0_mnet_);

    root_inv4_o1_mnet_ = dbModNet::create(root_mod_, "inv4_inst_o1_mnet");
    inv4_mod_level0_inst_o1_miterm_->connect(root_inv4_o1_mnet_);
    op1_bterm_->connect(root_inv4_o1_mnet_);

    root_inv4_o2_mnet_ = dbModNet::create(root_mod_, "inv4_inst_o2_mnet");
    inv4_mod_level0_inst_o2_miterm_->connect(root_inv4_o2_mnet_);
    op2_bterm_->connect(root_inv4_o2_mnet_);

    root_inv4_o3_mnet_ = dbModNet::create(root_mod_, "inv4_inst_o3_mnet");
    inv4_mod_level0_inst_o3_miterm_->connect(root_inv4_o3_mnet_);
    op3_bterm_->connect(root_inv4_o3_mnet_);

    // level 1 is next level down
    // The level 1 instance connections within the scope of the
    // inv4_mod_level0_master
    inv4_mod_level1_inst_i0_mnet_ = dbModNet::create(
        inv4_mod_level0_master_, "inv4_mod_level1_inst_i0_mnet");

    inv4_mod_level1_inst_i0_miterm_->connect(inv4_mod_level1_inst_i0_mnet_);
    inv4_mod_level0_master_i0_port_->connect(inv4_mod_level1_inst_i0_mnet_);

    inv4_mod_level1_inst_o0_mnet_ = dbModNet::create(
        inv4_mod_level0_master_, "inv4_mod_level1_inst_o1_mnet");
    inv4_mod_level1_inst_o0_miterm_->connect(inv4_mod_level1_inst_o0_mnet_);
    inv4_mod_level0_master_o0_port_->connect(inv4_mod_level1_inst_o0_mnet_);

    inv4_mod_level1_inst_o1_mnet_ = dbModNet::create(
        inv4_mod_level0_master_, "inv4_mod_level1_inst_o1_mnet");
    inv4_mod_level1_inst_o1_miterm_->connect(inv4_mod_level1_inst_o1_mnet_);
    inv4_mod_level0_master_o1_port_->connect(inv4_mod_level1_inst_o1_mnet_);

    inv4_mod_level1_inst_o2_mnet_ = dbModNet::create(
        inv4_mod_level0_master_, "inv4_mod_level1_inst_o2_mnet");
    inv4_mod_level1_inst_o2_miterm_->connect(inv4_mod_level1_inst_o2_mnet_);
    inv4_mod_level0_master_o2_port_->connect(inv4_mod_level1_inst_o2_mnet_);

    inv4_mod_level1_inst_o3_mnet_ = dbModNet::create(
        inv4_mod_level0_master_, "inv4_mod_level1_inst_o3_mnet");
    inv4_mod_level1_inst_o3_miterm_->connect(inv4_mod_level1_inst_o3_mnet_);
    inv4_mod_level0_master_o3_port_->connect(inv4_mod_level1_inst_o3_mnet_);

    // The level 2 instance connections within the scope of the
    // inv4_mod_level1_master level 2 is the cell which contains the 4 inverters
    inv4_mod_level2_inst_i0_mnet_ = dbModNet::create(
        inv4_mod_level1_master_, "inv4_mod_level2_inst_i0_mnet");

    inv4_mod_level2_inst_i0_miterm_->connect(inv4_mod_level2_inst_i0_mnet_);
    inv4_mod_level1_master_i0_port_->connect(inv4_mod_level2_inst_i0_mnet_);

    inv4_mod_level2_inst_o0_mnet_ = dbModNet::create(
        inv4_mod_level1_master_, "inv4_mod_level2_inst_o0_mnet");

    inv4_mod_level2_inst_o0_miterm_->connect(inv4_mod_level2_inst_o0_mnet_);
    inv4_mod_level1_master_o0_port_->connect(inv4_mod_level2_inst_o0_mnet_);

    inv4_mod_level2_inst_o1_mnet_ = dbModNet::create(
        inv4_mod_level1_master_, "inv4_mod_level2_inst_o1_mnet");

    inv4_mod_level2_inst_o1_miterm_->connect(inv4_mod_level2_inst_o1_mnet_);
    inv4_mod_level1_master_o1_port_->connect(inv4_mod_level2_inst_o1_mnet_);

    inv4_mod_level2_inst_o2_mnet_ = dbModNet::create(
        inv4_mod_level1_master_, "inv4_mod_level2_inst_o2_mnet");

    inv4_mod_level2_inst_o2_miterm_->connect(inv4_mod_level2_inst_o2_mnet_);
    inv4_mod_level1_master_o2_port_->connect(inv4_mod_level2_inst_o2_mnet_);

    inv4_mod_level2_inst_o3_mnet_ = dbModNet::create(
        inv4_mod_level1_master_, "inv4_mod_level2_inst_o3_mnet");

    inv4_mod_level2_inst_o3_miterm_->connect(inv4_mod_level2_inst_o3_mnet_);
    inv4_mod_level1_master_o3_port_->connect(inv4_mod_level2_inst_o3_mnet_);

    // Uncomment this to see the full design
    //    std::stringstream full_design;
    //    DbStrDebugHierarchy(block_, full_design);
    //    printf("The  design created (flat and hierarchical) %s\n",
    //	   full_design.str().c_str());
  }

  utl::UniquePtrWithDeleter<odb::dbDatabase> db_;
  sta::Unit* power_unit_;
  std::unique_ptr<sta::dbSta> sta_;
  sta::LibertyLibrary* library_;

  utl::Logger logger_;
  sta::dbNetwork* db_network_;

  dbBlock* block_;
  odb::dbLib* lib_;

  dbModule* root_mod_;
  dbMTerm* root_mod_i0_mterm_;
  dbMTerm* root_mod_i1_mterm_;
  dbMTerm* root_mod_i2_mterm_;
  dbMTerm* root_mod_i3_mterm_;
  dbMTerm* root_mod_o0_mterm_;

  dbBTerm* root_mod_i0_bterm_;
  dbBTerm* root_mod_i1_bterm_;
  dbBTerm* root_mod_i2_bterm_;
  dbBTerm* root_mod_i3_bterm_;
  dbBTerm* root_mod_o0_bterm_;

  dbModule* inv1_mod_master_;
  dbModInst* inv1_mod_inst_;
  dbModITerm* inv1_mod_inst_i0_miterm_;
  dbModITerm* inv1_mod_inst_o0_miterm_;
  dbModBTerm* inv1_mod_i0_port_;
  dbModBTerm* inv1_mod_o0_port_;
  dbModNet* inv1_mod_i0_modnet_;
  dbModNet* inv1_mod_o0_modnet_;

  dbModule* inv4_mod_level0_master_;
  dbModule* inv4_mod_level1_master_;
  dbModule* inv4_mod_level2_master_;

  dbModBTerm* inv4_mod_level0_master_i0_port_;
  dbModBTerm* inv4_mod_level0_master_o0_port_;
  dbModBTerm* inv4_mod_level0_master_o1_port_;
  dbModBTerm* inv4_mod_level0_master_o2_port_;
  dbModBTerm* inv4_mod_level0_master_o3_port_;

  dbModBTerm* inv4_mod_level1_master_i0_port_;
  dbModBTerm* inv4_mod_level1_master_o0_port_;
  dbModBTerm* inv4_mod_level1_master_o1_port_;
  dbModBTerm* inv4_mod_level1_master_o2_port_;
  dbModBTerm* inv4_mod_level1_master_o3_port_;

  dbModBTerm* inv4_mod_level2_master_i0_port_;
  dbModBTerm* inv4_mod_level2_master_o0_port_;
  dbModBTerm* inv4_mod_level2_master_o1_port_;
  dbModBTerm* inv4_mod_level2_master_o2_port_;
  dbModBTerm* inv4_mod_level2_master_o3_port_;

  dbModInst* inv4_mod_level0_inst_;
  dbModInst* inv4_mod_level1_inst_;
  dbModInst* inv4_mod_level2_inst_;

  dbModITerm* inv4_mod_level0_inst_i0_miterm_;
  dbModITerm* inv4_mod_level0_inst_o0_miterm_;
  dbModITerm* inv4_mod_level0_inst_o1_miterm_;
  dbModITerm* inv4_mod_level0_inst_o2_miterm_;
  dbModITerm* inv4_mod_level0_inst_o3_miterm_;

  dbModITerm* inv4_mod_level1_inst_i0_miterm_;
  dbModITerm* inv4_mod_level1_inst_o0_miterm_;
  dbModITerm* inv4_mod_level1_inst_o1_miterm_;
  dbModITerm* inv4_mod_level1_inst_o2_miterm_;
  dbModITerm* inv4_mod_level1_inst_o3_miterm_;

  dbModITerm* inv4_mod_level2_inst_i0_miterm_;
  dbModITerm* inv4_mod_level2_inst_o0_miterm_;
  dbModITerm* inv4_mod_level2_inst_o1_miterm_;
  dbModITerm* inv4_mod_level2_inst_o2_miterm_;
  dbModITerm* inv4_mod_level2_inst_o3_miterm_;

  dbInst* inv1_1_;
  dbInst* inv4_1_;
  dbInst* inv4_2_;
  dbInst* inv4_3_;
  dbInst* inv4_4_;

  dbNet* ip0_net_;
  dbNet* inv_op_net_;
  dbNet* op0_net_;
  dbNet* op1_net_;
  dbNet* op2_net_;
  dbNet* op3_net_;

  dbBTerm* ip0_bterm_;
  dbBTerm* op0_bterm_;
  dbBTerm* op1_bterm_;
  dbBTerm* op2_bterm_;
  dbBTerm* op3_bterm_;

  dbModITerm* inv1_1_ip_;
  dbModITerm* inv1_1_op_;

  dbModNet* root_inv1_i0_mnet_;
  dbModNet* root_inv1_o0_mnet_;

  dbModNet* root_inv4_o0_mnet_;
  dbModNet* root_inv4_o1_mnet_;
  dbModNet* root_inv4_o2_mnet_;
  dbModNet* root_inv4_o3_mnet_;

  // first input inv4_mod_level_i0 comes from root_inv1_o0_modnet
  dbModNet* inv4_mod_level0_inst_o0_mnet_;
  dbModNet* inv4_mod_level0_inst_o1_mnet_;
  dbModNet* inv4_mod_level0_inst_o2_mnet_;
  dbModNet* inv4_mod_level0_inst_o3_mnet_;

  dbModNet* inv4_mod_level1_inst_i0_mnet_;
  dbModNet* inv4_mod_level1_inst_o0_mnet_;
  dbModNet* inv4_mod_level1_inst_o1_mnet_;
  dbModNet* inv4_mod_level1_inst_o2_mnet_;
  dbModNet* inv4_mod_level1_inst_o3_mnet_;

  dbModNet* inv4_mod_level2_inst_i0_mnet_;
  dbModNet* inv4_mod_level2_inst_o0_mnet_;
  dbModNet* inv4_mod_level2_inst_o1_mnet_;
  dbModNet* inv4_mod_level2_inst_o2_mnet_;
  dbModNet* inv4_mod_level2_inst_o3_mnet_;

  dbITerm* inv1_1_inst_ip0_;
  dbITerm* inv1_1_inst_op0_;

  dbITerm* inv4_1_ip_;
  dbITerm* inv4_2_ip_;
  dbITerm* inv4_3_ip_;
  dbITerm* inv4_4_ip_;
  dbITerm* inv4_1_op_;
  dbITerm* inv4_2_op_;
  dbITerm* inv4_3_op_;
  dbITerm* inv4_4_op_;
};

TEST_F(TestHconn, ConnectionMade)
{
  //  std::stringstream str_str_initial;
  //  DbStrDebugHierarchy(block_, str_str_initial);
  //  printf("The initial design: %s\n", str_str_initial.str().c_str());

  // ECO test: get initial state before we start modifying
  // the design. Then at end we undo everything and
  // validate initial state preserved

  size_t initial_db_net_count = block_->getNets().size();
  size_t initial_mod_net_count = block_->getModNets().size();
  odb::dbDatabase::beginEco(block_);

  //
  //
  // Hierarchical connection test case:
  //----------------------------------
  //
  // add a new inverter to inv1 and connect to the fanout of the
  // existing one.
  // remove the driver of the 4th inverter in inv4_4 leaf
  // and connect to it via hierarchy
  //
  //
  // Before
  //--inv1_1 -----------inv4_1_ip0
  //        |----------inv4_2_ip0
  //        |----------inv4_3_ip0
  //        |----------inv4_4_ip0
  //
  //
  // After:
  //
  //--inv1_1-----------inv4_1_ip0
  //        |----------inv4_2_ip0
  //        |----------inv4_3_ip0
  //
  //--inv1_2---------inv4_4_ip0
  //
  //
  //
  // Objects referenced in "before"
  //
  // inv1_1 -- dbInst in inv1_1
  // inv1_1_inst_op0 -- op iterm on inv1
  // inv1_mod_master -- dbModule parent for inv1_1
  // inv4_mod_level2_master -- dbModule parent for leaf level inverters in inv4
  // inv4_4_ip0 -- ip iterm on 4th inverter
  //
  //

  //
  // first create the new inverter in inv1_mod_master
  // and get its pins ready for connection.
  // Note we are declaring these new objects outside
  // of the test struct so it is obvious what is in the test harness
  // and what is in the test
  //
  dbInst* inv1_2 = dbInst::create(block_,
                                  lib_->findMaster("INV_X1"),
                                  "inv1_mod_inst/inst2",
                                  false,
                                  inv1_mod_master_);
  dbITerm* inv1_2_inst_ip0 = inv1_2->findITerm("A");
  dbITerm* inv1_2_inst_op0 = inv1_2->findITerm("ZN");

  // Plumb in the new input of the inverter
  // This is the ip0_net, which is hooked to the top level bterm.
  // And to the top level modnet (which connects the modbterm on
  // inv1_mod_master) to the inverter. Note the dual connection:
  // one flat, one hierarchical.

  inv1_2_inst_ip0->connect(ip0_net_);             // flat world
  inv1_2_inst_ip0->connect(inv1_mod_i0_modnet_);  // hierarchical world

  // now disconnect the 4th inverter in inv4_4 (this is in level 2, the 3rd
  // level, of the hierarchy).
  // This kills both the dbNet (flat) connection
  // and the dbModNet (hierarchical) connection.

  inv4_4_ip_->disconnect();

  //
  // Make the flat connection.
  // We keep two worlds: the modnet world
  // and the dbNet world. The modnet world
  // exists at the edges eg where a module ports
  //(dbModITerms) connect to core gates (eg dbITerms).
  //
  // The flat world is always there.
  //
  //
  std::string flat_net_name = inv1_2->getName() + inv4_4_ip_->getName('/');
  std::string hier_net_name = "test_hier_" + flat_net_name;

  dbNet* flat_net = dbNet::create(block_, flat_net_name.c_str(), false);

  inv1_2_inst_op0->connect(flat_net);

  inv4_4_ip_->connect(flat_net);

  //
  // The hierarchical test: connect the output of the inverter in
  // inv1_2_inst_op0 to the input of the inverter inv4_4_ip0;
  //

  db_network_->hierarchicalConnect(
      inv1_2_inst_op0, inv4_4_ip_, hier_net_name.c_str());

  // Uncomment this to see the final design
  //  std::stringstream str_str_final;
  //  DbStrDebugHierarchy(block_, str_str_final);
  //  printf("The final design: %s\n", str_str_final.str().c_str());

  // Example of how to turn on the call backs for all the bterms/iterms
  // used by the sta
  dbSet<dbModNet> mod_nets = block_->getModNets();
  for (auto mnet : mod_nets) {
    sta::NetSet visited_nets;
    // given one mod net go get all its low level objects to
    // issue call back on
    sta::Net* cur_net = db_network_->dbToSta(mnet);
    sta::NetConnectedPinIterator* npi
        = db_network_->connectedPinIterator(cur_net);
    while (npi->hasNext()) {
      const sta::Pin* cur_pin = npi->next();
      odb::dbModITerm* moditerm;
      odb::dbITerm* iterm;
      odb::dbBTerm* bterm;
      db_network_->staToDb(cur_pin, iterm, bterm, moditerm);
      if (iterm) {
        db_network_->connectPinAfter(const_cast<sta::Pin*>(cur_pin));
        sta_->connectPinAfter(const_cast<sta::Pin*>(cur_pin));
      }
      if (bterm) {
        db_network_->connectPinAfter(const_cast<sta::Pin*>(cur_pin));
        sta_->connectPinAfter(const_cast<sta::Pin*>(cur_pin));
      }
      if (moditerm) {
        ;
      }
    }
  }

  // Get the final design state statistics

  size_t final_db_net_count = block_->getNets().size();
  size_t final_mod_net_count = block_->getModNets().size();

  EXPECT_EQ(initial_db_net_count, 6);
  EXPECT_EQ(initial_mod_net_count, 23);
  EXPECT_EQ(final_mod_net_count, 26);
  EXPECT_EQ(final_db_net_count, 7);

  //
  // Journalling test.
  // Undo everything and check initial state preserved
  //
  odb::dbDatabase::endEco(block_);
  odb::dbDatabase::undoEco(block_);

  size_t restored_db_net_count = block_->getNets().size();
  size_t restored_mod_net_count = block_->getModNets().size();

  EXPECT_EQ(restored_mod_net_count, initial_mod_net_count);
  EXPECT_EQ(restored_db_net_count, initial_db_net_count);

  // Test deletion of a dbModInst
  size_t num_mod_insts_before_delete = block_->getModInsts().size();
  size_t num_mod_iterms_before_delete = block_->getModITerms().size();
  dbModInst::destroy(inv1_mod_inst_);
  size_t num_mod_insts_after_delete = block_->getModInsts().size();
  size_t num_mod_iterms_after_delete = block_->getModITerms().size();

  // Test deletion of a dbModule
  size_t num_mods_before_delete = block_->getModules().size();
  size_t num_mod_bterms_before_delete = block_->getModBTerms().size();
  dbModule::destroy(inv1_mod_master_);
  size_t num_mods_after_delete = block_->getModules().size();
  size_t num_mod_bterms_after_delete = block_->getModBTerms().size();

  // we have killed one module instance
  EXPECT_EQ((num_mods_before_delete - num_mods_after_delete), 1);
  // we have killed one module
  EXPECT_EQ((num_mod_insts_before_delete - num_mod_insts_after_delete), 1);
  // we are deleting an inverter, so expect to delete 2 ports
  EXPECT_EQ((num_mod_bterms_before_delete - num_mod_bterms_after_delete), 2);
  // and the number of moditerms reduced should be the same (2)
  EXPECT_EQ((num_mod_iterms_before_delete - num_mod_iterms_after_delete),
            (num_mod_bterms_before_delete - num_mod_bterms_after_delete));
}

}  // namespace odb
