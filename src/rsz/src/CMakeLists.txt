############################################################################
##
## BSD 3-Clause License
##
## Copyright (c) 2019, James Cherry, Parallax Software, Inc.
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
## * Redistributions of source code must retain the above copyright notice, this
##   list of conditions and the following disclaimer.
##
## * Redistributions in binary form must reproduce the above copyright notice,
##   this list of conditions and the following disclaimer in the documentation
##   and/or other materials provided with the distribution.
##
## * Neither the name of the copyright holder nor the names of its
##   contributors may be used to endorse or promote products derived from
##   this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
## ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
## LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
## CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
## SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
## INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
## CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
## ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
## POSSIBILITY OF SUCH DAMAGE.
##
############################################################################

include("openroad")

swig_lib(NAME      rsz
         NAMESPACE rsz
         I_FILE    Resizer.i
         SWIG_INCLUDES ${OPENSTA_HOME}
         SCRIPTS   Resizer.tcl
)

add_library(rsz_lib
    BufferedNet.cc
    PreChecks.cc      
    RecoverPower.cc    
    RepairDesign.cc
    RepairHold.cc
    RepairSetup.cc
    Rebuffer.cc
    SpefWriter.cc
    SteinerTree.cc
    EstimateWireParasitics.cc
    Resizer.cc
    OdbCallBack.cc
)

target_sources(rsz
  PRIVATE
    MakeResizer.cc
    SteinerRenderer.cc
)

target_include_directories(rsz
  PUBLIC
    ../include
  PRIVATE
    # A side-effect of including OpenSTA swig files
    ${OPENSTA_HOME}
    ${OPENSTA_HOME}/include/sta
)

target_include_directories(rsz_lib
  PUBLIC
    ../include
)

target_link_libraries(rsz_lib
  PUBLIC
    stt_lib
    odb
    OpenSTA
    dbSta_lib
    grt_lib
    utl_lib
)

target_link_libraries(rsz
  PUBLIC
    rsz_lib
    stt
    odb
    OpenSTA
    dbSta
    grt
    gui
    utl_lib
)

messages(
  TARGET rsz
  OUTPUT_DIR ..
)
