// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

// Generator Code Begin Header
#pragma once

#include "dbCore.h"
#include "odb/dbIterator.h"
#include "odb/odb.h"

namespace odb {
class _dbScanInst;

class dbScanListScanInstItr : public dbIterator
{
 public:
  dbScanListScanInstItr(dbTable<_dbScanInst>* scan_inst_tbl)
  {
    _scan_inst_tbl = scan_inst_tbl;
  }

  bool reversible() override;
  bool orderReversed() override;
  void reverse(dbObject* parent) override;
  uint sequential() override;
  uint size(dbObject* parent) override;
  uint begin(dbObject* parent) override;
  uint end(dbObject* parent) override;
  uint next(uint id, ...) override;
  dbObject* getObject(uint id, ...) override;

 private:
  dbTable<_dbScanInst>* _scan_inst_tbl;
};

}  // namespace odb
// Generator Code End Header