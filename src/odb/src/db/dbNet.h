// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#pragma once

#include "dbCore.h"
#include "dbVector.h"
#include "odb/dbId.h"
#include "odb/dbTypes.h"
#include "odb/odb.h"

namespace odb {

class _dbITerm;
class _dbBTerm;
class _dbWire;
class _dbSWire;
class _dbCapNode;
class _dbRSeg;
class _dbCCSeg;
class _dbTechNonDefaultRule;
class _dbDatabase;
class _dbGroup;
class _dbGuide;
class _dbNetTrack;
class dbIStream;
class dbOStream;

struct _dbNetFlags
{
  dbSigType::Value _sig_type : 4;
  dbWireType::Value _wire_type : 4;
  uint _special : 1;
  uint _wild_connect : 1;
  uint _wire_ordered : 1;
  uint _unused2 : 1;       // free to reuse
  uint _disconnected : 1;  // this flag is only valid if wire_ordered == true
  uint _spef : 1;
  uint _select : 1;
  uint _mark : 1;
  uint _mark_1 : 1;
  uint _wire_altered : 1;
  uint _extracted : 1;
  uint _rc_graph : 1;
  uint _unused : 1;  // free to reuse
  uint _set_io : 1;
  uint _io : 1;
  uint _dont_touch : 1;
  uint _fixed_bump : 1;
  dbSourceType::Value _source : 4;
  uint _rc_disconnected : 1;
  uint _block_rule : 1;
  uint _has_jumpers : 1;
};

class _dbNet : public _dbObject
{
 public:
  enum Field  // dbJournal field name
  {
    FLAGS,
    NON_DEFAULT_RULE,
    TERM_EXTID,
    HEAD_CAPNODE,
    HEAD_RSEG,
    REVERSE_RSEG
  };

  // PERSISTANT-MEMBERS
  _dbNetFlags _flags;
  char* _name;
  union
  {
    float _gndc_calibration_factor;
    float _refCC;
  };
  union
  {
    float _cc_calibration_factor;
    float _dbCC;
    float _CcMatchRatio;
  };
  dbId<_dbNet> _next_entry;
  dbId<_dbITerm> _iterms;
  dbId<_dbBTerm> _bterms;
  dbId<_dbWire> _wire;
  dbId<_dbWire> _global_wire;
  dbId<_dbSWire> _swires;
  dbId<_dbCapNode> _cap_nodes;
  dbId<_dbRSeg> _r_segs;
  dbId<_dbTechNonDefaultRule> _non_default_rule;
  dbId<_dbGuide> guides_;
  dbId<_dbNetTrack> tracks_;
  dbVector<dbId<_dbGroup>> _groups;
  int _weight;
  int _xtalk;
  float _ccAdjustFactor;
  uint _ccAdjustOrder;
  // NON PERSISTANT-MEMBERS
  int _drivingIterm;

  _dbNet(_dbDatabase*);
  _dbNet(_dbDatabase*, const _dbNet& n);
  ~_dbNet();

  bool operator==(const _dbNet& rhs) const;
  bool operator!=(const _dbNet& rhs) const { return !operator==(rhs); }
  bool operator<(const _dbNet& rhs) const;
  void collectMemInfo(MemInfo& info);
};

dbOStream& operator<<(dbOStream& stream, const _dbNet& net);
dbIStream& operator>>(dbIStream& stream, _dbNet& net);

}  // namespace odb
