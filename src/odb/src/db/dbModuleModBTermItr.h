// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2020-2025, The OpenROAD Authors

// Generator Code Begin Header
#pragma once

#include "dbCore.h"
#include "odb/dbIterator.h"
#include "odb/odb.h"

namespace odb {
class _dbModBTerm;

class dbModuleModBTermItr : public dbIterator
{
 public:
  dbModuleModBTermItr(dbTable<_dbModBTerm>* modbterm_tbl)
  {
    _modbterm_tbl = modbterm_tbl;
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
  dbTable<_dbModBTerm>* _modbterm_tbl;
};

}  // namespace odb
   // Generator Code End Header
