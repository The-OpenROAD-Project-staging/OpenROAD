///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (c) 2020, The Regents of the University of California
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

// Generator Code Begin Header
#pragma once

#include "dbCore.h"
#include "dbVector.h"
#include "odb/dbSet.h"
#include "odb/odb.h"
// User Code Begin Includes
#include <string>
#include <unordered_map>

#include "dbHashTable.h"
#include "dbModulePortItr.h"
// User Code End Includes

namespace odb {
class dbIStream;
class dbOStream;
class _dbDatabase;
class _dbInst;
class _dbModInst;
class _dbModNet;
class _dbModBTerm;
// User Code Begin Classes
class dbInst;
// User Code End Classes

class _dbModule : public _dbObject
{
 public:
  _dbModule(_dbDatabase*);

  ~_dbModule();

  bool operator==(const _dbModule& rhs) const;
  bool operator!=(const _dbModule& rhs) const { return !operator==(rhs); }
  bool operator<(const _dbModule& rhs) const;
  void collectMemInfo(MemInfo& info);
  // User Code Begin Methods

  // This is only used when destroying an inst
  void removeInst(dbInst* inst);

  // User Code End Methods

  char* _name;
  dbId<_dbModule> _next_entry;
  dbId<_dbInst> _insts;
  dbId<_dbModInst> _mod_inst;
  dbId<_dbModInst> _modinsts;
  dbId<_dbModNet> _modnets;
  dbId<_dbModBTerm> _modbterms;

  // User Code Begin Fields
  // custom iterator for traversing ports
  // fast access
  std::unordered_map<std::string, dbId<_dbInst>> _dbinst_hash;
  std::unordered_map<std::string, dbId<_dbModInst>> _modinst_hash;
  std::unordered_map<std::string, dbId<_dbModBTerm>> _modbterm_hash;
  std::unordered_map<std::string, dbId<_dbModNet>> _modnet_hash;
  dbModulePortItr* _port_iter = nullptr;
  // User Code End Fields
};
dbIStream& operator>>(dbIStream& stream, _dbModule& obj);
dbOStream& operator<<(dbOStream& stream, const _dbModule& obj);
}  // namespace odb
   // Generator Code End Header
