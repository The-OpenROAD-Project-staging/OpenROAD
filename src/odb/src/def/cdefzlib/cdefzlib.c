/*******************************************************************************
 *******************************************************************************
 * Copyright 2014, Cadence Design Systems
 *
 * This  file  is  part  of  the  Cadence  LEF/DEF  Open   Source
 * Distribution,  Product Version 5.8.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *
 *        http://www.apache.org/licenses/LICENSE-2.0
 *
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 *    implied. See the License for the specific language governing
 *    permissions and limitations under the License.
 *
 * For updates, support, or to become part of the LEF/DEF Community,
 * check www.openeda.org for details.
 *******************************************************************************
 ******************************************************************************/

#include <sys/stat.h>
#include <sys/types.h>

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "defrReader.h"
#include "defzlib.h"
#include "zlib.h"

/*
 * Private functions:
 */
size_t defGZip_read(FILE* file, char* buf, size_t len)
{
  return gzread((gzFile) file, buf, (unsigned int) len);
}

/*
 * Public functions:
 */
defGZFile defGZipOpen(const char* gzipPath, const char* mode)
{
  defGZFile fptr;

  if (!gzipPath) {
    return NULL;
  }

  fptr = gzopen(gzipPath, mode);

  if (fptr) {
    /* successfully open the gzip file */
    /* set the read function to read from a compressed file */
    defrSetReadFunction(defGZip_read);
    return (defGZFile) fptr;
  } else {
    return NULL;
  }
}

int defGZipClose(defGZFile filePtr)
{
  return (gzclose((gzFile) filePtr));
}

int defrReadGZip(defGZFile file, const char* gzipFile, defiUserData uData)
{
  return defrRead((FILE*) file, gzipFile, uData, 1);
}
