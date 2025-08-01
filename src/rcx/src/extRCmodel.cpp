// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include <filesystem>
#include <limits>
#include <map>
#include <vector>

#include "parse.h"
#include "rcx/extRCap.h"
#include "rcx/extprocess.h"
#include "rcx/grids.h"
#include "utl/Logger.h"

namespace rcx {

using odb::dbRSeg;
using utl::RCX;

int extRCModel::getMaxMetIndexOverUnder(int met, int layerCnt)
{
  int n = 0;
  for (uint u = met - 1; u > 0; u--) {
    for (uint o = met + 1; o < layerCnt; o++) {
      int metIndex = extRCModel::getMetIndexOverUnder(met, u, o, layerCnt);
      if (n < metIndex) {
        n = metIndex;
      }
    }
  }
  return n;
}

static double lineSegment(double X, double x1, double x2, double y1, double y2)
{
  double slope = (y2 - y1) / (x2 - x1);
  double retVal = y2 - slope * (x2 - X);

  return retVal;
}

void extDistRC::interpolate(uint d, extDistRC* rc1, extDistRC* rc2)
{
  sep_ = d;
  coupling_
      = lineSegment(d, rc1->sep_, rc2->sep_, rc1->coupling_, rc2->coupling_);
  fringe_ = lineSegment(d, rc1->sep_, rc2->sep_, rc1->fringe_, rc2->fringe_);
  res_ = lineSegment(d, rc1->sep_, rc2->sep_, rc1->res_, rc2->res_);
}

double extDistRC::interpolate_res(uint d, extDistRC* rc2)
{
  return lineSegment(d, coupling_, rc2->coupling_, res_, rc2->res_);
}

void extDistRC::set(uint d, double cc, double fr, double a, double r)
{
  sep_ = d;
  coupling_ = cc;
  fringe_ = fr;
  diag_ = a;
  res_ = r;
}

void extDistRC::readRC(Ath__parser* parser, double dbFactor)
{
  sep_ = lround(dbFactor * 1000 * parser->getDouble(0));
  coupling_ = parser->getDouble(1) / dbFactor;
  fringe_ = parser->getDouble(2) / dbFactor;
  res_ = parser->getDouble(3) / dbFactor;
}

void extDistRC::readRC_res2(Ath__parser* parser, double dbFactor)
{
  sep_ = lround(dbFactor * 1000 * parser->getDouble(1));
  coupling_ = lround(dbFactor * 1000 * parser->getDouble(0));
  fringe_ = parser->getDouble(2) / dbFactor;
  res_ = parser->getDouble(3) / dbFactor;
}

double extDistRC::getCoupling()
{
  return coupling_;
}

double extDistRC::getFringe()
{
  return fringe_;
}

double extDistRC::getFringeW()
{
  return fringeW_;
}

double extDistRC::getDiag()
{
  return diag_;
}

double extDistRC::getRes()
{
  return res_;
}

int extDistRC::getSep()
{
  return sep_;
}

void extDistRC::setCoupling(const double coupling)
{
  coupling_ = coupling;
}

void extDistRC::setRes(const double res)
{
  res_ = res;
}

void extDistRC::setFringe(double fringe)
{
  fringe_ = fringe;
}

void extDistRC::setFringeW(double fringew)
{
  fringeW_ = fringew;
}

void extDistRC::writeRC()
{
  logger_->info(RCX,
                208,
                "{} {} {} {}  {}",
                0.001 * sep_,
                coupling_,
                fringe_,
                res_,
                coupling_ + fringe_);
}

void extDistRC::writeRC(FILE* fp, bool bin)
{
  fprintf(fp, "%g %g %g %g\n", 0.001 * sep_, coupling_, fringe_, res_);
}

void extRCTable::makeCapTableOver()
{
  _over = true;

  for (uint jj = 1; jj < _maxCnt1; jj++) {
    _inTable[jj] = new Ath__array1D<extDistRC*>*[jj];

    for (uint kk = 0; kk < jj; kk++) {
      _inTable[jj][kk] = new Ath__array1D<extDistRC*>(32);
    }
  }
}

void extRCTable::makeCapTableUnder()
{
  _over = false;
  for (uint jj = 1; jj < _maxCnt1; jj++) {
    _inTable[jj] = new Ath__array1D<extDistRC*>*[_maxCnt1];

    for (uint ii = 0; ii <= jj; ii++) {
      _inTable[jj][ii] = nullptr;
    }
    for (uint kk = jj + 1; kk < _maxCnt1; kk++) {
      _inTable[jj][kk] = new Ath__array1D<extDistRC*>(32);
    }
  }
}

extDistRCTable::extDistRCTable(uint distCnt)
{
  uint n = 16 * (distCnt / 16 + 1);
  measureTable_ = new Ath__array1D<extDistRC*>(n);
  measureInR_ = false;

  computeTable_ = nullptr;

  for (int i = 0; i < 16; i++) {
    measureTableR_[i] = nullptr;
    computeTableR_[i] = nullptr;
  }
}

extDistRCTable::~extDistRCTable()
{
  delete measureTable_;
  delete computeTable_;

  for (int i = 0; i < 16; i++) {
    if (measureTableR_[i] != measureTable_) {
      delete measureTableR_[i];
    }
    if (computeTableR_[i] != computeTable_) {
      delete computeTableR_[i];
    }
  }
}

uint extDistRCTable::mapExtrapolate(uint loDist,
                                    extDistRC* rc2,
                                    uint distUnit,
                                    AthPool<extDistRC>* rcPool)
{
  uint cnt = 0;
  uint d1 = loDist;
  uint d2 = rc2->sep_;

  for (uint d = d1; d <= d2; d += distUnit) {
    extDistRC* rc = rcPool->alloc();

    rc->sep_ = d;
    rc->coupling_ = rc2->coupling_;
    rc->fringe_ = rc2->fringe_;
    rc->res_ = rc2->res_;

    uint n = d / distUnit;

    computeTable_->set(n, rc);

    cnt++;
  }
  return cnt;
}

uint extDistRCTable::mapInterpolate(extDistRC* rc1,
                                    extDistRC* rc2,
                                    uint distUnit,
                                    int maxDist,
                                    AthPool<extDistRC>* rcPool)
{
  uint cnt = 0;
  uint d1 = rc1->sep_;
  uint d2 = rc2->sep_;

  if ((int) d2 > maxDist) {
    d2 = maxDist;
  }

  for (uint d = d1; d <= d2; d += distUnit) {
    extDistRC* rc = rcPool->alloc();

    rc->sep_ = d;
    rc->interpolate(rc->sep_, rc1, rc2);

    uint n = d / distUnit;

    computeTable_->set(n, rc);

    cnt++;
  }
  return cnt;
}

uint extDistRCTable::interpolate(uint distUnit,
                                 int maxDist,
                                 AthPool<extDistRC>* rcPool)
{
  uint cnt = measureTable_->getCnt();
  uint Cnt = cnt;
  if (cnt == 0) {
    return 0;
  }

  if (maxDist < 0) {
    extDistRC* lastRC = measureTable_->get(cnt - 1);
    maxDist = lastRC->sep_;
    if (maxDist == 100000) {
      maxDist = measureTable_->get(cnt - 2)->sep_;
      if (maxDist == 99000) {
        maxDist = measureTable_->get(cnt - 3)->sep_;
        Cnt = cnt - 2;
      } else {
        Cnt = cnt - 1;
      }
      extDistRC* rc31 = rcPool->alloc();
      rc31->set(0, 0.0, 0.0, 0.0, 0.0);
      measureTable_->set(31, rc31);
    }
  }

  makeComputeTable(maxDist, distUnit);

  mapExtrapolate(0, measureTable_->get(0), distUnit, rcPool);

  for (uint ii = 0; ii < Cnt - 1; ii++) {
    extDistRC* rc1 = measureTable_->get(ii);
    extDistRC* rc2 = measureTable_->get(ii + 1);

    mapInterpolate(rc1, rc2, distUnit, maxDist, rcPool);
  }
  if (Cnt != cnt) {
    extDistRC* rc1 = measureTable_->get(Cnt);
    extDistRC* rc = rcPool->alloc();
    rc->set(rc1->sep_, rc1->coupling_, rc1->fringe_, 0.0, rc1->res_);
    computeTable_->set(computeTable_->getSize() - 1, rc);
  }

  return computeTable_->getCnt();
}

uint extDistRCTable::writeRules(FILE* fp,
                                Ath__array1D<extDistRC*>* table,
                                double w,
                                bool bin)
{
  bool modify_last_line = false;  // Dimitri 09152023

  uint cnt = table->getCnt();
  if (cnt > 0) {
    extDistRC* rc1 = table->get(cnt - 1);
    if (rc1 != nullptr && modify_last_line) {
      rc1->set(rc1->sep_, 0, rc1->coupling_ + rc1->fringe_, 0.0, rc1->res_);
    }
  }

  fprintf(fp, "DIST count %d width %g\n", cnt, w);

  for (uint ii = 0; ii < cnt; ii++) {
    table->get(ii)->writeRC(fp, bin);
  }

  fprintf(fp, "END DIST\n");
  return cnt;
}

uint extDistRCTable::writeDiagRules(FILE* fp,
                                    Ath__array1D<extDistRC*>* table,
                                    double w1,
                                    double w2,
                                    double s,
                                    bool bin)
{
  uint cnt = table->getCnt();

  fprintf(fp,
          "DIST count %d width %g diag_width %g diag_dist %g\n",
          cnt,
          w1,
          w2,
          s);
  for (uint ii = 0; ii < cnt; ii++) {
    table->get(ii)->writeRC(fp, bin);
  }

  fprintf(fp, "END DIST\n");
  return cnt;
}

uint extDistRCTable::writeRules(FILE* fp, double w, bool compute, bool bin)
{
  if (compute) {
    return writeRules(fp, computeTable_, w, bin);
  }
  return writeRules(fp, measureTable_, w, bin);
}

uint extDistRCTable::writeDiagRules(FILE* fp,
                                    double w1,
                                    double w2,
                                    double s,
                                    bool compute,
                                    bool bin)
{
  if (compute) {
    return writeDiagRules(fp, computeTable_, w1, w2, s, bin);
  }
  return writeDiagRules(fp, measureTable_, w1, w2, s, bin);
}

uint extMetRCTable::readRCstats(Ath__parser* parser)
{
  uint cnt = 0;

  extMeasure m(logger_);

  while (parser->parseNextLine() > 0) {
    cnt++;

    m._overUnder = false;
    m._over = false;

    if ((parser->isKeyword(2, "OVER")) && (parser->isKeyword(4, "UNDER"))) {
      m._met = parser->getInt(1);
      m._overMet = parser->getInt(5);
      m._underMet = parser->getInt(3);
      m._overUnder = true;

      m._w_m = parser->getDouble(7);
      m._w_nm = lround(m._w_m * 1000);

      m._s_m = parser->getDouble(9);
      m._s_nm = lround(m._s_m * 1000);

      extDistRC* rc = _rcPoolPtr->alloc();

      rc->set(m._s_nm,
              parser->getDouble(10),
              parser->getDouble(11),
              0.0,
              parser->getDouble(12));

      m._tmpRC = rc;
    } else if (parser->isKeyword(2, "UNDER")
               || parser->isKeyword(2, "DIAGUNDER")) {
      m._met = parser->getInt(1);
      m._overMet = parser->getInt(4);
      m._underMet = -1;

      m._w_m = parser->getDouble(6);
      m._w_nm = lround(m._w_m * 1000);

      m._s_m = parser->getDouble(8);
      m._s_nm = lround(m._s_m * 1000);

      extDistRC* rc = _rcPoolPtr->alloc();

      rc->set(m._s_nm,
              parser->getDouble(9),
              parser->getDouble(10),
              0.0,
              parser->getDouble(11));

      m._tmpRC = rc;
    } else if (parser->isKeyword(2, "OVER")) {
      m._met = parser->getInt(1);
      m._underMet = parser->getInt(4);
      m._overMet = -1;
      m._over = true;

      m._w_m = parser->getDouble(6);
      m._w_nm = lround(m._w_m * 1000);

      m._s_m = parser->getDouble(8);
      m._s_nm = lround(m._s_m * 1000);

      extDistRC* rc = _rcPoolPtr->alloc();

      rc->set(m._s_nm,
              parser->getDouble(9),
              parser->getDouble(10),
              0.0,
              parser->getDouble(11));

      m._tmpRC = rc;
    }
    addRCw(&m);
  }
  mkWidthAndSpaceMappings();
  return cnt;
}

uint extDistRCTable::readRules_res2(Ath__parser* parser,
                                    AthPool<extDistRC>* rcPool,
                                    bool compute,
                                    bool bin,
                                    bool ignore,
                                    double dbFactor)
{
  parser->parseNextLine();
  uint cnt = parser->getInt(2);
  if (cnt < 32) {
    cnt = 32;
  }

  Ath__array1D<extDistRC*>* table = nullptr;
  if (!ignore) {
    table = new Ath__array1D<extDistRC*>(cnt);
  }

  if (!measureInR_) {
    delete measureTable_;
  }
  measureInR_ = false;

  Ath__array1D<extDistRC*>* table0 = new Ath__array1D<extDistRC*>(8);
  int cnt1 = 0;
  int kk = 0;
  extDistRC* rc0 = nullptr;
  while (parser->parseNextLine() > 0) {
    if (parser->isKeyword(0, "END")) {
      break;
    }

    if (ignore) {
      continue;
    }

    extDistRC* rc = rcPool->alloc();
    rc->readRC_res2(parser, dbFactor);
    table->add(rc);
    if (rc0 != nullptr && rc0->coupling_ != rc->coupling_) {
      measureTable_ = table0;
      if (table0->getCnt() > 1) {
        interpolate(4, -1, rcPool);
        computeTableR_[kk] = computeTable_;
      }
      measureTableR_[kk] = table0;
      kk++;

      table0 = new Ath__array1D<extDistRC*>(cnt1);
      cnt1 = 0;

      maxDist_ = rc0->sep_;
    }
    if (rc0 == nullptr) {
      table0->add(rc);
    } else if (cnt1 == 0) {
      table0->add(rc);
    } else if (rc0->res_ != rc->res_) {
      table0->add(rc);
    }
    cnt1++;
    rc0 = rc;
  }
  distCnt_ = kk + 1;
  measureTableR_[kk] = table0;
  measureTable_ = table;

  return cnt;
}

uint extDistRCTable::readRules(Ath__parser* parser,
                               AthPool<extDistRC>* rcPool,
                               bool compute,
                               bool bin,
                               bool ignore,
                               double dbFactor)
{
  parser->parseNextLine();
  uint cnt = parser->getInt(2);
  if (cnt < 32) {
    cnt = 32;
  }

  Ath__array1D<extDistRC*>* table = nullptr;
  if (!ignore) {
    table = new Ath__array1D<extDistRC*>(cnt);
  }

  while (parser->parseNextLine() > 0) {
    if (parser->isKeyword(0, "END")) {
      break;
    }

    if (ignore) {
      continue;
    }

    extDistRC* rc = rcPool->alloc();
    rc->readRC(parser, dbFactor);
    table->add(rc);
  }
  if (ignore) {
    return cnt;
  }

  if (!measureInR_) {
    delete measureTable_;
  }
  measureInR_ = false;

  measureTable_ = table;

  if (compute) {
    interpolate(4, -1, rcPool);
  }

  return cnt;
}

void extDistRCTable::ScaleRes(double SUB_MULT_RES,
                              Ath__array1D<extDistRC*>* table)
{
  uint cnt = table->getCnt();
  if (cnt == 0) {
    return;
  }

  extDistRC* rc_last = table->get(cnt - 1);

  for (uint jj = 0; jj < cnt; jj++) {
    extDistRC* rc = table->get(jj);
    double delta = rc->res_ - rc_last->res_;
    if (delta < 0) {
      delta = -delta;
    }
    if (delta > 0.000001) {
      continue;
    }

    rc->res_ *= SUB_MULT_RES;
  }
}

void extDistRCTable::makeComputeTable(uint maxDist, uint distUnit)
{
  unit_ = distUnit;  // in nm
  uint n = maxDist / distUnit;
  n = distUnit * (n / distUnit + 1);

  computeTable_ = new Ath__array1D<extDistRC*>(n + 1);
}

uint extDistRCTable::addMeasureRC(extDistRC* rc)
{
  return measureTable_->add(rc);
}

extDistRC* extDistRCTable::getRC_99()
{
  if (measureTable_ == nullptr) {
    return nullptr;
  }

  uint cnt = measureTable_->getCnt();
  if (cnt < 2) {
    return nullptr;
  }

  extDistRC* before_lastRC = measureTable_->get(cnt - 2);
  if (before_lastRC->sep_ == 99000) {
    return before_lastRC;
  }

  extDistRC* lastRC
      = measureTable_->getLast();  // assuming last is 100 equivalent to inf
  if (lastRC->sep_ == 99000) {
    return lastRC;
  }

  return nullptr;
}

extDistRC* extDistRCTable::getComputeRC(uint dist)
{
  if (measureTable_ == nullptr) {
    return nullptr;
  }

  if (measureTable_->getCnt() <= 0) {
    return nullptr;
  }

  extDistRC* firstRC = measureTable_->get(0);
  uint firstDist = firstRC->sep_;
  if (dist <= firstDist) {
    return firstRC;
  }

  if (measureTable_->getLast()->sep_ == 100000) {
    extDistRC* before_lastRC = measureTable_->getLast()
                               - 1;  // assuming last is 100 equivalent to inf
    uint lastDist = before_lastRC->sep_;

    if (lastDist == 99000) {
      before_lastRC = before_lastRC - 1;
    }

    lastDist = before_lastRC->sep_;
    if (dist >= lastDist) {    // send Inf dist
      if (dist == lastDist) {  // send Inf dist
        return before_lastRC;
      }
      if (dist <= 2 * lastDist) {  // send Inf dist

        uint cnt = measureTable_->getCnt();
        extDistRC* rc31 = measureTable_->geti(31);
        extDistRC* rc2 = measureTable_->get(cnt - 2);
        extDistRC* rc3 = measureTable_->get(cnt - 3);

        rc31->sep_ = dist;
        rc31->interpolate(dist, rc3, rc2);

        rc31->coupling_
            = (before_lastRC->coupling_ / dist) * before_lastRC->sep_;
        rc31->fringe_ = before_lastRC->fringe_;
        return rc31;
      }
      if (dist > lastDist) {  // send Inf dist
        return measureTable_->getLast();
      }
    }
  } else {
    extDistRC* before_lastRC
        = measureTable_->getLast();  // assuming last is 100 equivalent to inf
    uint lastDist = before_lastRC->sep_;
    if (dist >= lastDist - unit_ && lastDist > 0) {  // send Inf dist
      return measureTable_->getLast();
    }
  }

  uint n = dist / unit_;
  return computeTable_->geti(n);
}

uint extDistWidthRCTable::getWidthIndex(uint w)
{
  // To notify that the RC info for a particular pattern
  // is empty
  if (_lastWidth == -1) {
    return -1;
  }

  if ((int) w >= _lastWidth) {
    return _widthTable->getCnt() - 1;
  }

  int v = w - _firstWidth;
  if (v < 0) {
    return 0;
  }

  return _widthMapTable->geti(v / _modulo);
}

uint extDistWidthRCTable::getDiagWidthIndex(uint m, uint w)
{
  if (_lastDiagWidth == nullptr) {
    return -1;
  }

  if ((int) w >= _lastDiagWidth->geti(m)) {
    return _diagWidthTable[m]->getCnt() - 1;
  }

  int v = w - _firstDiagWidth->geti(m);
  if (v < 0) {
    return 0;
  }

  return _diagWidthMapTable[m]->geti(v / _modulo);
}

uint extDistWidthRCTable::getDiagDistIndex(uint m, uint s)
{
  if ((int) s >= _lastDiagDist->geti(m)) {
    return _diagDistTable[m]->getCnt() - 1;
  }

  int v = s - _firstDiagDist->geti(m);
  if (v < 0) {
    return 0;
  }

  return _diagDistMapTable[m]->geti(v / _modulo);
}

extDistWidthRCTable::extDistWidthRCTable(bool over,
                                         uint met,
                                         uint layerCnt,
                                         uint metCnt,
                                         uint maxWidthCnt,
                                         AthPool<extDistRC>* rcPool,
                                         bool OUREVERSEORDER)
    : _ouReadReverse(OUREVERSEORDER),
      _over(over),
      _layerCnt(layerCnt),
      _met(met),
      _modulo(4),
      _widthTableAllocFlag(false),
      _metCnt(metCnt),
      _widthCnt(maxWidthCnt),
      _rcPoolPtr(rcPool)
{
  _widthTable = new Ath__array1D<int>(maxWidthCnt);

  _lastWidth = std::numeric_limits<int>::max();

  _rcDistTable = new extDistRCTable**[_metCnt];
  uint jj;
  for (jj = 0; jj < _metCnt; jj++) {
    _rcDistTable[jj] = new extDistRCTable*[maxWidthCnt];
    for (uint ii = 0; ii < maxWidthCnt; ii++) {
      _rcDistTable[jj][ii] = new extDistRCTable(10);
    }
  }

  for (jj = 0; jj < diagDepth; jj++) {
    _diagWidthMapTable[jj] = nullptr;
    _diagDistMapTable[jj] = nullptr;
    _diagWidthTable[jj] = nullptr;
    _diagDistTable[jj] = nullptr;
  }
  _rc31 = rcPool->alloc();
}

void extDistWidthRCTable::createWidthMap()
{
  uint widthCnt = _widthTable->getCnt();
  if (widthCnt == 0) {
    return;
  }

  _firstWidth = _widthTable->get(0);
  _lastWidth = _widthTable->getLast();
  _modulo = 4;

  _widthTableAllocFlag = true;
  _widthMapTable = new Ath__array1D<uint>(10 * widthCnt);

  uint jj;
  for (jj = 0; jj < widthCnt - 1; jj++) {
    double v1 = _widthTable->get(jj);
    double v2 = _widthTable->get(jj + 1);

    int w1 = lround(v1);
    int w2 = lround(v2);

    for (int w = w1; w <= w2; w += _modulo) {
      if (w >= _lastWidth) {
        continue;
      }

      uint n = 0;
      int v = w - _firstWidth;
      if (v > 0) {
        n = v / _modulo;
      }

      _widthMapTable->set(n, jj);
    }
  }
}

void extDistWidthRCTable::makeWSmapping()
{
  createWidthMap();

  for (uint jj = 0; jj < _metCnt; jj++) {
    for (uint ii = 0; ii < _widthTable->getCnt(); ii++) {
      _rcDistTable[jj][ii]->interpolate(4, -1, _rcPoolPtr);
    }
  }
}

extDistWidthRCTable::extDistWidthRCTable(bool dummy,
                                         uint met,
                                         uint layerCnt,
                                         uint widthCnt,
                                         bool OUREVERSEORDER)
    : _ouReadReverse(OUREVERSEORDER),
      _layerCnt(layerCnt),
      _met(met),
      _widthTableAllocFlag(true),
      _metCnt(layerCnt),
      _widthCnt(widthCnt)
{
  _widthTable = new Ath__array1D<int>(widthCnt);
  for (uint ii = 0; ii < widthCnt; ii++) {
    _widthTable->add(0);
  }

  _rcDistTable = new extDistRCTable**[_metCnt];
  uint jj;
  for (jj = 0; jj < _metCnt; jj++) {
    _rcDistTable[jj] = new extDistRCTable*[widthCnt];
    for (uint ii = 0; ii < widthCnt; ii++) {
      _rcDistTable[jj][ii] = new extDistRCTable(1);
    }
  }
  for (jj = 0; jj < diagDepth; jj++) {
    _diagWidthMapTable[jj] = nullptr;
    _diagDistMapTable[jj] = nullptr;
    _diagWidthTable[jj] = nullptr;
    _diagDistTable[jj] = nullptr;
  }
}

extDistWidthRCTable::extDistWidthRCTable(bool over,
                                         uint met,
                                         uint layerCnt,
                                         uint metCnt,
                                         Ath__array1D<double>* widthTable,
                                         AthPool<extDistRC>* rcPool,
                                         bool OUREVERSEORDER,
                                         double dbFactor)
    : _ouReadReverse(OUREVERSEORDER),
      _over(over),
      _layerCnt(layerCnt),
      _met(met),
      _metCnt(layerCnt)
{
  for (uint jj = 0; jj < diagDepth; jj++) {
    _diagWidthMapTable[jj] = nullptr;
    _diagDistMapTable[jj] = nullptr;
    _diagWidthTable[jj] = nullptr;
    _diagDistTable[jj] = nullptr;
  }

  // dkf 09202024 skip width map table when  not knowing number of widths is not
  // know in advance after reading rules of different width width mapping should
  // be re-done before writting rules
  bool skip_width_map_table = widthTable == nullptr;

  if (!skip_width_map_table && widthTable->getCnt() == 0) {
    return;
  }
  _metCnt = metCnt;
  _widthTableAllocFlag = true;
  int widthCnt = 12;
  _widthTable = new Ath__array1D<int>(widthCnt);

  if (!skip_width_map_table) {
    widthCnt = widthTable->getCnt();
    for (uint ii = 0; ii < widthCnt; ii++) {
      int w = lround(dbFactor * 1000 * widthTable->get(ii));
      _widthTable->add(w);
    }
    if (widthCnt > 0) {
      _firstWidth = _widthTable->get(0);
      _lastWidth = _widthTable->get(widthCnt - 1);
    }
    _modulo = 4;

    _widthMapTable = new Ath__array1D<uint>(10 * widthCnt);

    uint jj;
    for (jj = 0; jj < widthCnt - 1; jj++) {
      double v1 = _widthTable->get(jj);
      double v2 = _widthTable->get(jj + 1);

      int w1 = lround(v1);
      int w2 = lround(v2);

      for (int w = w1; w <= w2; w += _modulo) {
        if (w >= _lastWidth) {
          continue;
        }

        uint n = 0;
        int v = w - _firstWidth;
        if (v > 0) {
          n = v / _modulo;
        }

        _widthMapTable->set(n, jj);
      }
    }
  }

  _widthCnt = widthCnt;
  _rcDistTable = new extDistRCTable**[_metCnt];
  for (uint jj = 0; jj < _metCnt; jj++) {
    _rcDistTable[jj] = new extDistRCTable*[widthCnt];
    for (uint ii = 0; ii < widthCnt; ii++) {
      _rcDistTable[jj][ii] = new extDistRCTable(10);
    }
  }
  _rcPoolPtr = rcPool;

  _rc31 = rcPool->alloc();
}

extDistWidthRCTable::extDistWidthRCTable(bool over,
                                         uint met,
                                         uint layerCnt,
                                         uint metCnt,
                                         Ath__array1D<double>* widthTable,
                                         int diagWidthCnt,
                                         int diagDistCnt,
                                         AthPool<extDistRC>* rcPool,
                                         bool OUREVERSEORDER,
                                         double dbFactor)
    : _ouReadReverse(OUREVERSEORDER),
      _over(over),
      _layerCnt(layerCnt),
      _met(met),
      _modulo(4),
      _widthTableAllocFlag(true),
      _metCnt(metCnt),
      _diagWidthCnt(diagWidthCnt),
      _diagDistCnt(diagDistCnt),
      _rcPoolPtr(rcPool)
{
  uint widthCnt = widthTable->getCnt();
  _widthTable = new Ath__array1D<int>(widthCnt);
  for (uint ii = 0; ii < widthCnt; ii++) {
    int w = lround(dbFactor * 1000 * widthTable->get(ii));
    _widthTable->add(w);
  }
  for (uint i = 0; i < layerCnt; i++) {
    _diagWidthTable[i] = new Ath__array1D<int>(diagWidthCnt);
    _diagDistTable[i] = new Ath__array1D<int>(diagDistCnt);
    _diagWidthMapTable[i] = new Ath__array1D<uint>(10 * diagWidthCnt);
    _diagDistMapTable[i] = new Ath__array1D<uint>(10 * diagDistCnt);
  }

  for (uint i = layerCnt; i < diagDepth; i++) {
    _diagWidthTable[i] = nullptr;
    _diagDistTable[i] = nullptr;
    _diagWidthMapTable[i] = nullptr;
    _diagDistMapTable[i] = nullptr;
  }

  _firstWidth = _widthTable->get(0);
  _lastWidth = _widthTable->get(widthCnt - 1);
  _firstDiagWidth = new Ath__array1D<int>(layerCnt);
  _lastDiagWidth = new Ath__array1D<int>(layerCnt);
  _firstDiagDist = new Ath__array1D<int>(layerCnt);
  _lastDiagDist = new Ath__array1D<int>(layerCnt);

  _widthMapTable = new Ath__array1D<uint>(10 * widthCnt);
  uint jj;
  for (jj = 0; jj < widthCnt - 1; jj++) {
    double v1 = _widthTable->get(jj);
    double v2 = _widthTable->get(jj + 1);

    int w1 = lround(v1);
    int w2 = lround(v2);

    for (int w = w1; w <= w2; w += _modulo) {
      if (w >= _lastWidth) {
        continue;
      }

      uint n = 0;
      int v = w - _firstWidth;
      if (v > 0) {
        n = v / _modulo;
      }

      _widthMapTable->set(n, jj);
    }
  }

  _widthCnt = widthCnt;
  _rcDiagDistTable = new extDistRCTable****[_metCnt];
  for (jj = 0; jj < _metCnt; jj++) {
    _rcDiagDistTable[jj] = new extDistRCTable***[widthCnt];
    for (uint ii = 0; ii < widthCnt; ii++) {
      _rcDiagDistTable[jj][ii] = new extDistRCTable**[diagWidthCnt];
      for (int kk = 0; kk < diagWidthCnt; kk++) {
        _rcDiagDistTable[jj][ii][kk] = new extDistRCTable*[diagDistCnt];
        for (int ll = 0; ll < diagDistCnt; ll++) {
          _rcDiagDistTable[jj][ii][kk][ll] = new extDistRCTable(10);
        }
      }
    }
  }

  _rc31 = rcPool->alloc();
}

void extDistWidthRCTable::setDiagUnderTables(
    uint met,
    Ath__array1D<double>* diagWidthTable,
    Ath__array1D<double>* diagDistTable,
    double dbFactor)
{
  uint diagWidthCnt = diagWidthTable->getCnt();
  _diagWidthTable[met]->resetCnt();
  uint ii;
  for (ii = 0; ii < diagWidthCnt; ii++) {
    int w = lround(dbFactor * 1000 * diagWidthTable->get(ii));
    _diagWidthTable[met]->add(w);
  }
  _firstDiagWidth->set(met, _diagWidthTable[met]->get(0));
  _lastDiagWidth->set(met, _diagWidthTable[met]->get(diagWidthCnt - 1));
  uint diagDistCnt = diagDistTable->getCnt();
  _diagDistTable[met]->resetCnt();
  for (ii = 0; ii < diagDistCnt; ii++) {
    int s = lround(dbFactor * 1000 * diagDistTable->get(ii));
    _diagDistTable[met]->add(s);
  }
  _firstDiagDist->set(met, _diagDistTable[met]->get(0));
  _lastDiagDist->set(met, _diagDistTable[met]->get(diagDistCnt - 1));
  uint jj;
  for (jj = 0; jj < diagWidthCnt - 1; jj++) {
    double v1 = _diagWidthTable[met]->get(jj);
    double v2 = _diagWidthTable[met]->get(jj + 1);

    int w1 = lround(v1);
    int w2 = lround(v2);

    for (int w = w1; w <= w2; w += _modulo) {
      if (w >= _lastDiagWidth->geti(met)) {
        continue;
      }

      uint n = 0;
      int v = w - _firstDiagWidth->geti(met);
      if (v > 0) {
        n = v / _modulo;
      }

      _diagWidthMapTable[met]->set(n, jj);
    }
  }
  for (jj = 0; jj < diagDistCnt - 1; jj++) {
    double v1 = _diagDistTable[met]->get(jj);
    double v2 = _diagDistTable[met]->get(jj + 1);

    int s1 = lround(v1);
    int s2 = lround(v2);

    for (int s = s1; s <= s2; s += _modulo) {
      if (s >= _lastDiagDist->geti(met)) {
        continue;
      }

      int d = (s2 - s1) / 2;

      uint n = 0;
      int v = s - _firstDiagDist->geti(met);
      if (v > 0) {
        n = v / _modulo;
      }

      if (v < s1 + d) {
        _diagDistMapTable[met]->set(n, jj);
      } else {
        _diagDistMapTable[met]->set(n, jj + 1);
      }
    }
  }
}

extDistWidthRCTable::~extDistWidthRCTable()
{
  uint ii, jj, kk, ll;
  if (_rcDistTable != nullptr) {
    for (jj = 0; jj < _metCnt; jj++) {
      for (ii = 0; ii < _widthCnt; ii++) {
        delete _rcDistTable[jj][ii];
      }
      delete[] _rcDistTable[jj];
    }
    delete[] _rcDistTable;
  }

  if (_rcDiagDistTable != nullptr) {
    for (jj = 0; jj < _metCnt; jj++) {
      for (ii = 0; ii < _widthCnt; ii++) {
        for (kk = 0; kk < _diagWidthCnt; kk++) {
          for (ll = 0; ll < _diagDistCnt; ll++) {
            delete _rcDiagDistTable[jj][ii][kk][ll];
          }
          delete[] _rcDiagDistTable[jj][ii][kk];
        }
        delete[] _rcDiagDistTable[jj][ii];
      }
      delete[] _rcDiagDistTable[jj];
    }
    delete[] _rcDiagDistTable;
  }

  delete _widthTable;
  delete _widthMapTable;
  delete _firstDiagWidth;
  delete _lastDiagWidth;
  delete _firstDiagDist;
  delete _lastDiagDist;
  for (uint i = 0; i < _layerCnt; i++) {
    if (_diagWidthTable[i] != nullptr) {
      delete _diagWidthTable[i];
    }
    if (_diagDistTable[i] != nullptr) {
      delete _diagDistTable[i];
    }
    if (_diagWidthMapTable[i] != nullptr) {
      delete _diagWidthMapTable[i];
    }
    if (_diagDistMapTable[i] != nullptr) {
      delete _diagDistMapTable[i];
    }
  }
}

uint extDistWidthRCTable::writeWidthTable(FILE* fp, bool bin)
{
  uint widthCnt = _widthTable->getCnt();
  fprintf(fp, "WIDTH Table %d entries: ", widthCnt);
  for (uint ii = 0; ii < widthCnt; ii++) {
    fprintf(fp, " %g", 0.001 * _widthTable->get(ii));
  }
  fprintf(fp, "\n");
  return widthCnt;
}

uint extDistWidthRCTable::writeDiagWidthTable(FILE* fp, uint met, bool bin)
{
  uint diagWidthCnt = _diagWidthTable[met]->getCnt();
  fprintf(fp, "DIAG_WIDTH Table %d entries: ", diagWidthCnt);
  for (uint ii = 0; ii < diagWidthCnt; ii++) {
    fprintf(fp, " %g", 0.001 * _diagWidthTable[met]->get(ii));
  }
  fprintf(fp, "\n");
  return diagWidthCnt;
}

void extDistWidthRCTable::writeDiagTablesCnt(FILE* fp, uint met, bool bin)
{
  uint diagWidthCnt = _diagWidthTable[met]->getCnt();
  uint diagDistCnt = _diagDistTable[met]->getCnt();
  fprintf(fp, "DIAG_WIDTH Table Count: %d\n", diagWidthCnt);
  fprintf(fp, "DIAG_DIST Table Count: %d\n", diagDistCnt);
}

uint extDistWidthRCTable::writeDiagDistTable(FILE* fp, uint met, bool bin)
{
  uint diagDistCnt = _diagDistTable[met]->getCnt();
  fprintf(fp, "DIAG_DIST Table %d entries: ", diagDistCnt);
  for (uint ii = 0; ii < diagDistCnt; ii++) {
    fprintf(fp, " %g", 0.001 * _diagDistTable[met]->get(ii));
  }
  fprintf(fp, "\n");
  return diagDistCnt;
}

uint extDistWidthRCTable::writeRulesOver(FILE* fp, bool bin)
{
  uint cnt = 0;
  fprintf(fp, "\nMetal %d OVER\n", _met);

  writeWidthTable(fp, bin);
  uint widthCnt = _widthTable->getCnt();

  for (uint ii = 0; ii < _met; ii++) {
    fprintf(fp, "\nMetal %d OVER %d\n", _met, ii);

    for (uint jj = 0; jj < widthCnt; jj++) {
      cnt += _rcDistTable[ii][jj]->writeRules(
          fp, 0.001 * _widthTable->get(jj), false, bin);
    }
  }
  return cnt;
}

uint extDistWidthRCTable::readMetalHeader(Ath__parser* parser,
                                          uint& met,
                                          const char* keyword,
                                          bool bin,
                                          bool ignore)
{
  if (!(parser->parseNextLine() > 0)) {
    return 0;
  }

  if (parser->isKeyword(0, "Metal") && (strcmp(parser->get(2), keyword) == 0)) {
    met = parser->getInt(1);
    return 1;
  }

  return 0;
}

uint extDistWidthRCTable::readRulesOver(Ath__parser* parser,
                                        uint widthCnt,
                                        bool bin,
                                        bool ignore,
                                        const char* OVER,
                                        double dbFactor)
{
  bool res = strcmp(OVER, "RESOVER") == 0;
  uint cnt = 0;
  for (uint ii = 0; ii < _met; ii++) {
    uint met = 0;
    if (readMetalHeader(parser, met, OVER, bin, ignore) <= 0) {
      return 0;
    }

    parser->getInt(3);

    for (uint jj = 0; jj < widthCnt; jj++) {
      if (res) {
        if (!ignore) {
          cnt += _rcDistTable[ii][jj]->readRules_res2(
              parser, _rcPoolPtr, true, bin, ignore, dbFactor);
        } else {
          cnt += _rcDistTable[0][0]->readRules_res2(
              parser, _rcPoolPtr, true, bin, ignore, dbFactor);
        }
      } else {
        if (!ignore) {
          cnt += _rcDistTable[ii][jj]->readRules(
              parser, _rcPoolPtr, true, bin, ignore, dbFactor);
        } else {
          cnt += _rcDistTable[0][0]->readRules(
              parser, _rcPoolPtr, true, bin, ignore, dbFactor);
        }
      }
    }
  }
  return cnt;
}

uint extDistWidthRCTable::readRulesUnder(Ath__parser* parser,
                                         uint widthCnt,
                                         bool bin,
                                         bool ignore,
                                         double dbFactor)
{
  uint cnt = 0;
  for (uint ii = _met + 1; ii < _layerCnt; ii++) {
    uint met = 0;
    if (readMetalHeader(parser, met, "UNDER", bin, ignore) <= 0) {
      return 0;
    }

    uint metIndex = getMetIndexUnder(ii);
    if (ignore) {
      metIndex = 0;
    }

    parser->getInt(3);

    for (uint jj = 0; jj < widthCnt; jj++) {
      cnt += _rcDistTable[metIndex][jj]->readRules(
          parser, _rcPoolPtr, true, bin, ignore, dbFactor);
    }
  }
  return cnt;
}

uint extDistWidthRCTable::readRulesDiagUnder(Ath__parser* parser,
                                             uint widthCnt,
                                             uint diagWidthCnt,
                                             uint diagDistCnt,
                                             bool bin,
                                             bool ignore,
                                             double dbFactor)
{
  uint cnt = 0;
  for (uint ii = _met + 1; ii < _met + 5 && ii < _layerCnt; ii++) {
    uint met = 0;
    if (readMetalHeader(parser, met, "DIAGUNDER", bin, ignore) <= 0) {
      return 0;
    }
    Ath__array1D<double>* dwTable = nullptr;
    Ath__array1D<double>* ddTable = nullptr;
    parser->parseNextLine();
    dwTable = parser->readDoubleArray("DIAG_WIDTH", 4);
    parser->parseNextLine();
    ddTable = parser->readDoubleArray("DIAG_DIST", 4);
    uint diagWidthCnt = dwTable->getCnt();
    uint diagDistCnt = ddTable->getCnt();
    uint metIndex = getMetIndexUnder(ii);

    if (!ignore) {
      setDiagUnderTables(metIndex, dwTable, ddTable);
    }

    parser->getInt(3);

    for (uint jj = 0; jj < widthCnt; jj++) {
      for (uint kk = 0; kk < diagWidthCnt; kk++) {
        for (uint ll = 0; ll < diagDistCnt; ll++) {
          if (!ignore) {
            cnt += _rcDiagDistTable[metIndex][jj][kk][ll]->readRules(
                parser, _rcPoolPtr, true, bin, ignore, dbFactor);
          } else {
            cnt += _rcDistTable[0][0]->readRules(
                parser, _rcPoolPtr, true, bin, ignore, dbFactor);
          }
        }
      }
    }
    delete dwTable;
    delete ddTable;
  }
  return cnt;
}

uint extDistWidthRCTable::readRulesDiagUnder(Ath__parser* parser,
                                             uint widthCnt,
                                             bool bin,
                                             bool ignore,
                                             double dbFactor)
{
  uint cnt = 0;
  for (uint ii = _met + 1; ii < _layerCnt; ii++) {
    uint met = 0;
    if (readMetalHeader(parser, met, "DIAGUNDER", bin, ignore) <= 0) {
      return 0;
    }

    uint metIndex = getMetIndexUnder(ii);
    parser->getInt(3);

    for (uint jj = 0; jj < widthCnt; jj++) {
      cnt += _rcDistTable[metIndex][jj]->readRules(
          parser, _rcPoolPtr, true, bin, ignore, dbFactor);
    }
  }
  return cnt;
}

uint extDistWidthRCTable::readRulesOverUnder(Ath__parser* parser,
                                             uint widthCnt,
                                             bool bin,
                                             bool ignore,
                                             double dbFactor)
{
  uint cnt = 0;
  for (uint u = 1; u < _met; u++) {
    for (uint o = _met + 1; o < _layerCnt; o++) {
      uint mOver = o;

      uint met = 0;
      if (readMetalHeader(parser, met, "OVER", bin, ignore) <= 0) {
        return 0;
      }

      if (_ouReadReverse) {
        mOver = parser->getInt(5);
      }

      uint mUnder = parser->getInt(3);

      // Commented out this code per Dimitris...
      // The variable mOver is already defined above...
      // uint mOver= parser->getInt(5);

      int metIndex = 0;
      if (!ignore) {
        metIndex = extRCModel::getMetIndexOverUnder(
            _met, mUnder, mOver, _layerCnt, _metCnt);
      }
      int mcnt = 0;
      for (uint jj = 0; jj < widthCnt; jj++) {
        if (!ignore) {
          mcnt += _rcDistTable[metIndex][jj]->readRules(
              parser, _rcPoolPtr, true, bin, ignore, dbFactor);
        } else {
          mcnt += _rcDistTable[0][0]->readRules(
              parser, _rcPoolPtr, true, bin, ignore, dbFactor);
        }
      }
      cnt += mcnt;
    }
  }
  return cnt;
}

uint extDistWidthRCTable::getMetIndexUnder(uint mOver)
{
  return mOver - _met - 1;
}

uint extDistWidthRCTable::writeRulesUnder(FILE* fp, bool bin)
{
  uint cnt = 0;
  fprintf(fp, "\nMetal %d UNDER\n", _met);

  writeWidthTable(fp, bin);
  uint widthCnt = _widthTable->getCnt();

  for (uint ii = _met + 1; ii < _layerCnt; ii++) {
    fprintf(fp, "\nMetal %d UNDER %d\n", _met, ii);

    uint metIndex = getMetIndexUnder(ii);

    for (uint jj = 0; jj < widthCnt; jj++) {
      cnt += _rcDistTable[metIndex][jj]->writeRules(
          fp, 0.001 * _widthTable->get(jj), false, bin);
    }
  }
  return cnt;
}
uint extDistWidthRCTable::writeRulesOverUnder(FILE* fp, bool bin)
{
  uint cnt = 0;
  fprintf(fp, "\nMetal %d OVERUNDER\n", _met);

  writeWidthTable(fp, bin);
  uint widthCnt = _widthTable->getCnt();

  for (uint mUnder = 1; mUnder < _met; mUnder++) {
    for (uint mOver = _met + 1; mOver < _layerCnt; mOver++) {
      fprintf(fp, "\nMetal %d OVER %d UNDER %d\n", _met, mUnder, mOver);

      int metIndex = extRCModel::getMetIndexOverUnder(
          _met, mUnder, mOver, _layerCnt, _metCnt);
      assert(metIndex >= 0);

      for (uint jj = 0; jj < widthCnt; jj++) {
        cnt += _rcDistTable[metIndex][jj]->writeRules(
            fp, 0.001 * _widthTable->get(jj), false, bin);
      }
    }
  }
  return cnt;
}
extMetRCTable::extMetRCTable(uint layerCnt,
                             AthPool<extDistRC>* rcPool,
                             Logger* logger,
                             bool OUREVERSEORDER)
{
  logger_ = logger;
  _layerCnt = layerCnt;
  _wireCnt = 2;

  _resOver = new extDistWidthRCTable*[layerCnt];
  _capOver = new extDistWidthRCTable*[layerCnt];
  _capDiagUnder = new extDistWidthRCTable*[layerCnt];
  _capUnder = new extDistWidthRCTable*[layerCnt];
  _capOverUnder = new extDistWidthRCTable*[layerCnt];

  _capOver_open = allocTable();
  _capUnder_open = allocTable();
  _capOverUnder_open = allocTable();

  for (uint ii = 0; ii < layerCnt; ii++) {
    _resOver[ii] = nullptr;
    _capOver[ii] = nullptr;
    _capDiagUnder[ii] = nullptr;
    _capUnder[ii] = nullptr;
    _capOverUnder[ii] = nullptr;
  }
  _rcPoolPtr = rcPool;
  _rate = -1000.0;

  _OUREVERSEORDER = OUREVERSEORDER;
}

extMetRCTable::~extMetRCTable()
{
  for (uint ii = 0; ii < _layerCnt; ii++) {
    delete _capUnder[ii];
    delete _capDiagUnder[ii];
    delete _resOver[ii];
    delete _capOver[ii];
    delete _capOverUnder[ii];

    for (uint jj = 0; jj < _wireCnt; jj++) {
      delete _capOver_open[ii][jj];
      delete _capUnder_open[ii][jj];
      delete _capOverUnder_open[ii][jj];
    }
    delete[] _capOver_open[ii];
    delete[] _capUnder_open[ii];
    delete[] _capOverUnder_open[ii];
  }
  delete[] _resOver;
  delete[] _capOver;
  delete[] _capDiagUnder;
  delete[] _capUnder;
  delete[] _capOverUnder;
  delete[] _capOver_open;
  delete[] _capUnder_open;
  delete[] _capOverUnder_open;
}

void extMetRCTable::allocDiagUnderTable(uint met,
                                        Ath__array1D<double>* wTable,
                                        int diagWidthCnt,
                                        int diagDistCnt,
                                        double dbFactor)
{
  delete _capDiagUnder[met];
  _capDiagUnder[met] = new extDistWidthRCTable(false,
                                               met,
                                               _layerCnt,
                                               _layerCnt - met - 1,
                                               wTable,
                                               diagWidthCnt,
                                               diagDistCnt,
                                               _rcPoolPtr,
                                               _OUREVERSEORDER,
                                               dbFactor);
}

void extMetRCTable::setDiagUnderTables(uint met,
                                       uint overMet,
                                       Ath__array1D<double>* diagWTable,
                                       Ath__array1D<double>* diagSTable,
                                       double dbFactor)
{
  _capDiagUnder[met]->setDiagUnderTables(
      overMet, diagWTable, diagSTable, dbFactor);
}

void extMetRCTable::allocDiagUnderTable(uint met,
                                        Ath__array1D<double>* wTable,
                                        double dbFactor)
{
  delete _capDiagUnder[met];
  _capDiagUnder[met] = new extDistWidthRCTable(false,
                                               met,
                                               _layerCnt,
                                               _layerCnt - met - 1,
                                               wTable,
                                               _rcPoolPtr,
                                               _OUREVERSEORDER,
                                               dbFactor);
}

void extMetRCTable::allocUnderTable(uint met,
                                    Ath__array1D<double>* wTable,
                                    double dbFactor)
{
  delete _capUnder[met];
  _capUnder[met] = new extDistWidthRCTable(false,
                                           met,
                                           _layerCnt,
                                           _layerCnt - met - 1,
                                           wTable,
                                           _rcPoolPtr,
                                           _OUREVERSEORDER,
                                           dbFactor);
}

void extMetRCTable::allocOverUnderTable(uint met,
                                        Ath__array1D<double>* wTable,
                                        double dbFactor)
{
  if (met < 2) {
    return;
  }

  int n = extRCModel::getMaxMetIndexOverUnder(met, _layerCnt);
  delete _capOverUnder[met];
  _capOverUnder[met] = new extDistWidthRCTable(false,
                                               met,
                                               _layerCnt,
                                               n + 1,
                                               wTable,
                                               _rcPoolPtr,
                                               _OUREVERSEORDER,
                                               dbFactor);
}

extRCTable::extRCTable(bool over, uint layerCnt)
{
  _maxCnt1 = layerCnt + 1;
  _inTable = new Ath__array1D<extDistRC*>**[_maxCnt1];

  if (over) {
    makeCapTableOver();
  } else {
    makeCapTableUnder();
  }
}

extRCTable::~extRCTable()
{
  for (uint jj = 1; jj < _maxCnt1; jj++) {
    const uint max = _over ? jj : _maxCnt1;

    for (uint kk = 0; kk < max; kk++) {
      delete _inTable[jj][kk];
    }
    delete[] _inTable[jj];
  }

  delete[] _inTable;
}

extDistRC* extRCTable::getCapOver(uint met, uint metUnder)
{
  return _inTable[met][metUnder]->get(0);
}

double extRCModel::getTotCapOverSub(uint met)
{
  extDistRC* rc = _capOver->getCapOver(met, 0);
  return rc->getFringe();
}

extDistRC* extDistRCTable::getRC_index(int n)
{
  if (n < 0) {
    return nullptr;
  }
  int cnt = measureTable_->getCnt();
  if (n >= cnt) {
    return nullptr;
  }
  return measureTable_->get(n);
}

extDistRC* extDistRCTable::getLastRC()
{
  int cnt = measureTable_->getCnt();
  return getRC_index(cnt - 1);
}

extDistRC* extDistRCTable::getRC(uint s, bool compute)
{
  if (compute) {
    return getComputeRC(s);
  }
  return nullptr;
}

extDistRC* extDistWidthRCTable::getFringeRC(uint mou, uint w, int index_dist)
{
  const int wIndex = getWidthIndex(w);
  if ((wIndex < 0) || (wIndex >= (int) _widthTable->getCnt())) {
    return nullptr;
  }

  if (mou >= _metCnt || wIndex >= (int) _widthTable->getCnt()
      || _rcDistTable[mou][wIndex] == nullptr) {
    return nullptr;
  }

  extDistRC* rc;
  if (index_dist < 0) {
    rc = _rcDistTable[mou][wIndex]->getLastRC();
  } else {
    rc = _rcDistTable[mou][wIndex]->getRC_index(index_dist);
  }
  return rc;
}

extDistRC* extDistWidthRCTable::getLastWidthFringeRC(uint mou)
{
  if (mou >= _metCnt) {
    return nullptr;
  }

  int wIndex = _widthTable->getCnt() - 1;

  if (wIndex >= (int) _widthTable->getCnt()
      || _rcDistTable[mou][wIndex] == nullptr) {
    return nullptr;
  }

  return _rcDistTable[mou][wIndex]->getLastRC();
}

extDistRC* extDistWidthRCTable::getRC(uint mou, uint w, uint s)
{
  int wIndex = getWidthIndex(w);
  if (wIndex < 0) {
    return nullptr;
  }

  return _rcDistTable[mou][wIndex]->getRC(s, true);
}

extDistRC* extDistWidthRCTable::getRC(uint mou,
                                      uint w,
                                      uint dw,
                                      uint ds,
                                      uint s)
{
  int wIndex = getWidthIndex(w);
  if (wIndex < 0) {
    return nullptr;
  }
  int dwIndex = getDiagWidthIndex(mou, dw);
  if (dwIndex < 0) {
    return nullptr;
  }
  int dsIndex = getDiagDistIndex(mou, ds);
  if (dsIndex < 0) {
    return nullptr;
  }
  return _rcDiagDistTable[mou][wIndex][dwIndex][dsIndex]->getRC(s, true);
}

extDistRC* extDistWidthRCTable::getRC_99(uint mou, uint w, uint dw, uint ds)
{
  int wIndex = getWidthIndex(w);
  if (wIndex < 0) {
    return nullptr;
  }

  int dwIndex = getDiagWidthIndex(mou, dw);
  if (dwIndex < 0) {
    return nullptr;
  }

  int dsIndex = getDiagDistIndex(mou, ds);
  if (dsIndex < 0) {
    return nullptr;
  }

  uint s2 = _diagDistTable[mou]->get(dsIndex);
  extDistRC* rc2 = _rcDiagDistTable[mou][wIndex][dwIndex][dsIndex]->getRC_99();
  if (dsIndex == 0) {
    return rc2;
  }

  if ((int) ds == _diagDistTable[mou]->get(dsIndex)) {
    return rc2;
  }

  _rc31->sep_ = ds;

  uint lastDist = _lastDiagDist->geti(mou);
  if (ds > lastDist) {  // extrapolate
    _rc31->fringe_ = (rc2->fringe_ / ds) * lastDist;

    return _rc31;
  }
  // interpolate;
  uint s1 = _diagDistTable[mou]->get(dsIndex - 1);

  if (ds <= (s1 - s2) / 4) {  // too close!
    return rc2;
  }

  extDistRC* rc1
      = _rcDiagDistTable[mou][wIndex][dwIndex][dsIndex - 1]->getRC_99();

  _rc31->fringe_ = lineSegment(ds, s1, s2, rc1->fringe_, rc2->fringe_);

  return _rc31;
}

double extRCModel::getFringeOver(uint met, uint mUnder, uint w, uint s)
{
  extDistRC* rc = _modelTable[_tmpDataRate]->_capOver[met]->getRC(mUnder, w, s);

  return rc->getFringe();
}

double extRCModel::getCouplingOver(uint met, uint mUnder, uint w, uint s)
{
  extDistRC* rc = _modelTable[_tmpDataRate]->_capOver[met]->getRC(mUnder, w, s);

  return rc->getCoupling();
}

extDistRC* extRCModel::getOverRC(extMeasure* m)
{
  if (_modelTable[_tmpDataRate] == nullptr
      || _modelTable[_tmpDataRate]->_capOver[m->_met] == nullptr) {
    return nullptr;
  }
  extDistRC* rc = _modelTable[_tmpDataRate]->_capOver[m->_met]->getRC(
      m->_underMet, m->_width, m->_dist);

  return rc;
}

extDistRC* extRCModel::getUnderRC(extMeasure* m)
{
  uint n = getUnderIndex(m);
  if (_modelTable[_tmpDataRate] == nullptr
      || _modelTable[_tmpDataRate]->_capUnder[m->_met] == nullptr) {
    return nullptr;
  }
  extDistRC* rc = _modelTable[_tmpDataRate]->_capUnder[m->_met]->getRC(
      n, m->_width, m->_dist);

  return rc;
}
extDistRC* extRCModel::getUnderRC(int met, int overMet, int width, int dist)
{
  uint n = overMet - met - 1;
  if (_modelTable[_tmpDataRate] == nullptr
      || _modelTable[_tmpDataRate]->_capUnder[met] == nullptr) {
    return nullptr;
  }
  extDistRC* rc
      = _modelTable[_tmpDataRate]->_capUnder[met]->getRC(n, width, dist);

  return rc;
}

extDistRC* extRCModel::getOverUnderRC(extMeasure* m)
{
  uint maxOverUnderIndex
      = _modelTable[_tmpDataRate]->_capOverUnder[m->_met]->_metCnt;
  uint n = getOverUnderIndex(m, maxOverUnderIndex);
  extDistRC* rc = _modelTable[_tmpDataRate]->_capOverUnder[m->_met]->getRC(
      n, m->_width, m->_dist);

  return rc;
}
extDistRC* extRCModel::getOverUnderRC(uint met,
                                      uint underMet,
                                      int overMet,
                                      int width,
                                      int dist)
{
  uint maxOverUnderIndex
      = _modelTable[_tmpDataRate]->_capOverUnder[met]->_metCnt;
  // uint n = getOverUnderIndex(m, maxOverUnderIndex);
  uint n = getMetIndexOverUnder(
      met, underMet, overMet, _layerCnt, maxOverUnderIndex);
  extDistRC* rc
      = _modelTable[_tmpDataRate]->_capOverUnder[met]->getRC(n, width, dist);

  return rc;
}

extDistRC* extRCModel::getOverFringeRC(uint met, uint underMet, uint width)
{
  if (met >= _layerCnt) {
    return nullptr;
  }

  extDistRC* rc
      = _modelTable[_tmpDataRate]->_capOver[met]->getFringeRC(underMet, width);

  return rc;
}

extDistRC* extMetRCTable::getOverFringeRC(extMeasure* m, int index_dist)
{
  if (m->_met >= (int) _layerCnt) {
    return nullptr;
  }

  extDistRC* rc
      = _capOver[m->_met]->getFringeRC(m->_underMet, m->_width, index_dist);

  return rc;
}

extDistRC* extMetRCTable::getOverFringeRC_last(int met, int width)
{
  if (met >= (int) _layerCnt) {
    return nullptr;
  }

  extDistRC* rc = _capOver[met]->getFringeRC(0, width, -1);

  return rc;
}

extDistRC* extRCModel::getOverFringeRC(extMeasure* m)
{
  if (m->_met >= (int) _layerCnt) {
    return nullptr;
  }

  extDistRC* rc = _modelTable[_tmpDataRate]->_capOver[m->_met]->getFringeRC(
      m->_underMet, m->_width);

  return rc;
}

extDistRC* extRCModel::getUnderFringeRC(extMeasure* m)
{
  uint n = getUnderIndex(m);
  if (_modelTable[_tmpDataRate] == nullptr
      || _modelTable[_tmpDataRate]->_capUnder[m->_met] == nullptr) {
    return nullptr;
  }
  extDistRC* rc = _modelTable[_tmpDataRate]->_capUnder[m->_met]->getFringeRC(
      n, m->_width);

  return rc;
}

extDistRC* extRCModel::getOverUnderFringeRC(extMeasure* m)
{
  uint maxCnt = _modelTable[_tmpDataRate]->_capOverUnder[m->_met]->_metCnt;
  uint n = getOverUnderIndex(m, maxCnt);
  if (_modelTable[_tmpDataRate] == nullptr
      || _modelTable[_tmpDataRate]->_capOverUnder[m->_met] == nullptr) {
    return nullptr;
  }
  extDistRC* rc
      = _modelTable[_tmpDataRate]->_capOverUnder[m->_met]->getFringeRC(
          n, m->_width);

  return rc;
}

extDistRC* extMeasure::getOverUnderFringeRC(extMetRCTable* rcModel)
{
  uint maxCnt = _currentModel->getMaxCnt(_met);
  int n = extRCModel::getMetIndexOverUnder(
      _met, _underMet, _overMet, _layerCnt, maxCnt);

  if (rcModel == nullptr || rcModel->_capOverUnder[_met] == nullptr) {
    return nullptr;
  }

  extDistRC* rc = rcModel->_capOverUnder[_met]->getFringeRC(n, _width);

  return rc;
}

extDistRC* extMeasure::getOverUnderRC(extMetRCTable* rcModel)
{
  uint maxCnt = _currentModel->getMaxCnt(_met);
  int n = extRCModel::getMetIndexOverUnder(
      _met, _underMet, _overMet, _layerCnt, maxCnt);

  extDistRC* rc = nullptr;
  if (_dist < 0) {
    rc = rcModel->_capOverUnder[_met]->getFringeRC(n, _width);
  } else {
    rc = rcModel->_capOverUnder[_met]->getRC(n, _width, _dist);
  }

  return rc;
}

extDistRC* extMeasure::getOverRC(extMetRCTable* rcModel)
{
  if (_met >= (int) _layerCnt) {
    return nullptr;
  }

  extDistRC* rc = nullptr;
  if (_dist < 0) {
    rc = rcModel->_capOver[_met]->getFringeRC(_underMet, _width);
  } else {
    rc = rcModel->_capOver[_met]->getRC(_underMet, _width, _dist);
  }

  return rc;
}

uint extMeasure::getUnderIndex(uint overMet)
{
  return overMet - _met - 1;
}

uint extMeasure::getUnderIndex()
{
  return _overMet - _met - 1;
}

extDistRC* extMeasure::getUnderLastWidthDistRC(extMetRCTable* rcModel,
                                               uint overMet)
{
  if (rcModel->_capUnder[_met] == nullptr) {
    return nullptr;
  }

  uint n = getUnderIndex(overMet);

  return rcModel->_capUnder[_met]->getLastWidthFringeRC(n);
}

extDistRC* extMeasure::getUnderRC(extMetRCTable* rcModel)
{
  if (rcModel->_capUnder[_met] == nullptr) {
    return nullptr;
  }

  uint n = getUnderIndex();

  extDistRC* rc = nullptr;
  if (_dist < 0) {
    rc = rcModel->_capUnder[_met]->getFringeRC(n, _width);
  } else {
    rc = rcModel->_capUnder[_met]->getRC(n, _width, _dist);
  }

  return rc;
}

extDistRC* extMeasure::getVerticalUnderRC(extMetRCTable* rcModel,
                                          uint diagDist,
                                          uint tgtWidth,
                                          uint overMet)
{
  if (rcModel->_capDiagUnder[_met] == nullptr) {
    return getUnderRC(rcModel);  // DELETE
    return nullptr;
  }

  uint n = getUnderIndex(overMet);

  extDistRC* rc
      = rcModel->_capDiagUnder[_met]->getRC_99(n, _width, tgtWidth, diagDist);

  return rc;
}

double extMeasure::getDiagUnderCC(extMetRCTable* rcModel,
                                  uint dist,
                                  uint overMet)
{
  if (rcModel->_capDiagUnder[_met] == nullptr) {
    return 0.0;
  }

  uint n = getUnderIndex(overMet);

  extDistRC* rc = rcModel->_capDiagUnder[_met]->getRC(n, _width, dist);

  if (rc != nullptr) {
    if (IsDebugNet()) {
      int dbUnit = _extMain->_block->getDbUnitsPerMicron();
      rc->printDebugRC_diag(_met, overMet, 0, _width, dist, dbUnit, logger_);
    }
    return rc->fringe_;
  }
  return 0.0;
}

double extMeasure::getDiagUnderCC(extMetRCTable* rcModel,
                                  uint diagWidth,
                                  uint diagDist,
                                  uint overMet)
{
  if (rcModel->_capDiagUnder[_met] == nullptr) {
    return 0.0;
  }

  uint n = getUnderIndex(overMet);

  extDistRC* rc = rcModel->_capDiagUnder[_met]->getRC(
      n, _width, diagWidth, diagDist, _dist);

  if (rc != nullptr) {
    return rc->fringe_;
  }
  return 0.0;
}

extDistRC* extMeasure::getDiagUnderCC2(extMetRCTable* rcModel,
                                       uint diagWidth,
                                       uint diagDist,
                                       uint overMet)
{
  if (rcModel->_capDiagUnder[_met] == nullptr) {
    return nullptr;
  }

  uint n = getUnderIndex(overMet);

  extDistRC* rc = rcModel->_capDiagUnder[_met]->getRC(
      n, _width, diagWidth, diagDist, _dist);

  if (rc == nullptr) {
    return nullptr;
  }
  return rc;
}

double extRCModel::getRes(uint met)
{
  if (met > 13) {
    return 0;
  }

  extDistRC* rc = _capOver->getCapOver(met, 0);
  if (rc == nullptr) {
    return 0;
  }

  return rc->getRes();
}

uint extRCTable::addCapOver(uint met, uint metUnder, extDistRC* rc)
{
  return _inTable[met][metUnder]->add(rc);
}

extRCModel::extRCModel(uint layerCnt, const char* name, Logger* logger)
{
  logger_ = logger;
  _layerCnt = layerCnt;
  strcpy(_name, name);
  _resOver = new extRCTable(true, layerCnt);
  _capOver = new extRCTable(true, layerCnt);
  _capUnder = new extRCTable(false, layerCnt);
  _capDiagUnder = new extRCTable(false, layerCnt);
  _rcPoolPtr = new AthPool<extDistRC>(1024);
  _process = nullptr;
  _maxMinFlag = false;

  _wireDirName = new char[2048];
  _topDir = new char[1024];
  _patternName = new char[1024];
  _parser = new Ath__parser(logger_);
  _solverFileName = new char[1024];
  _wireFileName = new char[1024];
  _capLogFP = nullptr;
  _logFP = nullptr;

  _readCapLog = false;
  _commentFlag = false;

  _modelCnt = 0;
  _dataRateTable = nullptr;
  _modelTable = nullptr;
  _tmpDataRate = 0;
  _extMain = nullptr;
  _ruleFileName = nullptr;
  _diagModel = 0;
  _verticalDiag = false;
  _metLevel = 0;
}

extRCModel::extRCModel(const char* name, Logger* logger)
{
  logger_ = logger;
  _layerCnt = 0;
  strcpy(_name, name);
  _resOver = nullptr;
  _capOver = nullptr;
  _capUnder = nullptr;
  _capDiagUnder = nullptr;
  _rcPoolPtr = new AthPool<extDistRC>(1024);
  _process = nullptr;
  _maxMinFlag = false;

  _wireDirName = new char[2048];
  _topDir = new char[1024];
  _patternName = new char[1024];
  _parser = new Ath__parser(logger_);
  _solverFileName = new char[1024];
  _wireFileName = new char[1024];
  _capLogFP = nullptr;
  _logFP = nullptr;

  _readCapLog = false;
  _commentFlag = false;

  _modelCnt = 0;
  _modelTable = nullptr;
  _tmpDataRate = 0;

  _noVariationIndex = -1;

  _extMain = nullptr;
  _ruleFileName = nullptr;
  _diagModel = 0;
  _verticalDiag = false;
  _metLevel = 0;
}

extRCModel::~extRCModel()
{
  free(_ruleFileName);
  delete _resOver;
  delete _capOver;
  delete _capUnder;
  delete _capDiagUnder;
  delete _rcPoolPtr;

  delete[] _wireDirName;
  delete[] _topDir;
  delete[] _patternName;
  delete _parser;
  delete[] _solverFileName;
  delete[] _wireFileName;

  for (uint ii = 0; ii < _modelCnt; ii++) {
    delete _modelTable[ii];
  }

  delete[] _modelTable;
  delete _dataRateTable;
}

void extRCModel::setExtMain(extMain* x)
{
  _extMain = x;
}

extProcess* extRCModel::getProcess()
{
  return _process;
}

void extRCModel::setProcess(extProcess* p)
{
  _process = p;
}

// extMetRCTable holds one RC model per process corner
void extRCModel::createModelTable(uint n, uint layerCnt)
{
  _layerCnt = layerCnt;
  _modelCnt = n;

  _dataRateTable = new Ath__array1D<double>(_modelCnt);
  _modelTable = new extMetRCTable*[_modelCnt];
  for (uint jj = 0; jj < _modelCnt; jj++) {
    _modelTable[jj]
        = new extMetRCTable(_layerCnt, _rcPoolPtr, logger_, _OUREVERSEORDER);
  }
}

void extRCModel::setDataRateTable(uint met)
{
  if (_process == nullptr) {
    return;
  }
  _maxMinFlag = _process->getMaxMinFlag();
  bool thickVarFlag = _process->getThickVarFlag();
  extVariation* xvar = _process->getVariation(met);

  if (xvar != nullptr) {
    Ath__array1D<double>* dTable = xvar->getDataRateTable();

    createModelTable(dTable->getCnt() + 1, _layerCnt);

    _dataRateTable->add(0.0);
    for (uint ii = 0; ii < dTable->getCnt(); ii++) {
      _dataRateTable->add(dTable->get(ii));
    }

  } else if (_maxMinFlag) {
    createModelTable(3, _layerCnt);
    for (uint i = 0; i < 3; i++) {
      _dataRateTable->add(i);
    }
  } else if (thickVarFlag) {
    Ath__array1D<double>* dTable = _process->getDataRateTable(1);
    createModelTable(dTable->getCnt(), _layerCnt);
    for (uint ii = 0; ii < dTable->getCnt(); ii++) {
      _dataRateTable->add(dTable->get(ii));
    }
  } else {
    createModelTable(1, _layerCnt);
    _dataRateTable->add(0.0);
  }
  _tmpDataRate = 0;
}

uint extRCModel::addLefTotRC(uint met, uint underMet, double fr, double r)
{
  extDistRC* rc = _rcPoolPtr->alloc();
  rc->set(0, 0.0, fr, 0.0, r);

  uint n = _capOver->addCapOver(met, underMet, rc);
  return n;
}

uint extRCModel::addCapOver(uint met,
                            uint underMet,
                            uint d,
                            double cc,
                            double fr,
                            double a,
                            double r)
{
  extDistRC* rc = _rcPoolPtr->alloc();
  rc->set(d, cc, fr, a, r);

  uint n = _capOver->addCapOver(met, underMet, rc);
  return n;
}

extMeasure::extMeasure(utl::Logger* logger)
    : _create_net_util(logger), logger_(logger)
{
  _met = -1;
  _underMet = -1;
  _overMet = -1;

  _w_m = 0.0;
  _s_m = 0.0;
  _w_nm = 0;
  _s_nm = 0;
  _r = 0.0;
  _t = 0.0;
  _h = 0.0;
  _w2_m = 0.0;
  _s2_m = 0.0;
  _w2_nm = 0;
  _s2_nm = 0;

  _topWidth = 0.0;
  _botWidth = 0.0;
  _teff = 0.0;
  _heff = 0.0;
  _seff = 0.0;

  _varFlag = false;
  _diag = false;
  _verticalDiag = false;
  _plate = false;
  _metExtFlag = false;

  for (auto& ii : _rc) {
    ii = new extDistRC();
    ii->setLogger(logger_);
  }

  _tmpRC = _rc[0];

  _capTable = nullptr;
  _ll[0] = 0;
  _ll[1] = 0;
  _ur[0] = 0;
  _ur[1] = 0;
  /*
    _maxCapNodeCnt = 100;
    for (int n = 0; n < (int) _maxCapNodeCnt; n++) {
      for (int k = 0; k < (int) _maxCapNodeCnt; k++) {
        _capMatrix[n][k] = 0.0;
      }
    }
    */
  _extMain = nullptr;

  _2dBoxPool = new AthPool<ext2dBox>(1024);

  _lenOUPool = nullptr;
  _lenOUtable = nullptr;

  _totCCcnt = 0;
  _totSmallCCcnt = 0;
  _totBigCCcnt = 0;
  _totSignalSegCnt = 0;
  _totSegCnt = 0;

  _tmpDstTable = new Ath__array1D<SEQ*>(32);
  _tmpSrcTable = new Ath__array1D<SEQ*>(32);
  _diagTable = new Ath__array1D<SEQ*>(32);
  _tmpTable = new Ath__array1D<SEQ*>(32);
  _ouTable = new Ath__array1D<SEQ*>(32);
  _overTable = new Ath__array1D<SEQ*>(32);
  _underTable = new Ath__array1D<SEQ*>(32);

  _seqPool = new AthPool<SEQ>(1024);

  _dgContextFile = nullptr;
  _diagFlow = false;
  _rotatedGs = false;
  _sameNetFlag = false;
}

void extMeasure::allocOUpool()
{
  _lenOUPool = new AthPool<extLenOU>(128);
  _lenOUtable = new Ath__array1D<extLenOU*>(128);
}

extMeasure::~extMeasure()
{
  for (auto& ii : _rc) {
    delete ii;
  }

  delete _tmpDstTable;
  delete _tmpSrcTable;
  delete _diagTable;
  delete _tmpTable;
  delete _ouTable;
  delete _overTable;
  delete _underTable;

  delete _seqPool;

  delete _2dBoxPool;
  if (_lenOUPool != nullptr) {
    delete _lenOUPool;
    delete _lenOUtable;
  }
}

void extMeasure::setMets(int m, int u, int o)
{
  _met = m;
  _underMet = u;
  _overMet = o;
  _over = false;
  _overUnder = false;
  if ((u > 0) && (o > 0)) {
    _overUnder = true;
  } else if ((u >= 0) && (o < 0)) {
    _over = true;
  }
}

void extMeasure::setTargetParams(double w,
                                 double s,
                                 double r,
                                 double t,
                                 double h,
                                 double w2,
                                 double s2)
{
  _w_m = w;
  _s_m = s;
  _w_nm = lround(1000 * w);
  _s_nm = lround(1000 * s);
  _r = r;
  _t = t;
  _h = h;
  if (w2 > 0.0) {
    _w2_m = w2;
    _w2_nm = lround(1000 * w2);
  } else {
    _w2_m = _w_m;
    _w2_nm = _w_nm;
  }
  if (s2 > 0.0 || (s2 == 0.0 && _diag)) {
    long int n2 = _s2_nm = lround(1000 * s2);
    n2 = (n2 / 10) * 10;

    _s2_m = 0.001 * n2;
    // _s2_nm = lround(1000 * s2);
    _s2_nm = n2;
  } else {
    _s2_m = _s_m;
    _s2_nm = _s_nm;
  }
}

void extMeasure::setEffParams(double wTop, double wBot, double teff)
{
  _topWidth = wTop;
  _botWidth = wBot;
  _teff = teff;
  _heff = _h;
  if (!_metExtFlag && _s_m != 99) {
    _seff = _w_m + _s_m - wTop;
  } else {
    _seff = _s_m;
  }
}

extDistRC* extMeasure::addRC(extDistRC* rcUnit, uint len, uint jj)
{
  if (rcUnit == nullptr) {
    return nullptr;
  }
  int dbUnit = _extMain->_block->getDbUnitsPerMicron();
  if (IsDebugNet()) {
    rcUnit->printDebugRC(
        _met, _overMet, _underMet, _width, _dist, dbUnit, logger_);
  }

  if (_sameNetFlag) {  // TO OPTIMIZE
    _rc[jj]->fringe_ += 0.5 * rcUnit->fringe_ * len;
  } else {
    _rc[jj]->fringe_ += rcUnit->fringe_ * len;

    if (_dist > 0) {  // dist based
      _rc[jj]->coupling_ += rcUnit->coupling_ * len;
    }
  }

  _rc[jj]->res_ += rcUnit->res_ * len;
  if (IsDebugNet()) {
    _rc[jj]->printDebugRC_sum(len, dbUnit, logger_);
  }
  return rcUnit;
}

extDistRC* extMeasure::computeOverUnderRC(uint len)
{
  extDistRC* rcUnit = nullptr;

  for (uint ii = 0; ii < _metRCTable.getCnt(); ii++) {
    extMetRCTable* rcModel = _metRCTable.get(ii);

    rcUnit = getOverUnderRC(rcModel);

    addRC(rcUnit, len, ii);
  }
  return rcUnit;
}

extDistRC* extMeasure::computeOverRC(uint len)
{
  extDistRC* rcUnit = nullptr;

  for (uint ii = 0; ii < _metRCTable.getCnt(); ii++) {
    extMetRCTable* rcModel = _metRCTable.get(ii);

    rcUnit = getOverRC(rcModel);

    addRC(rcUnit, len, ii);
  }
  return rcUnit;
}

extDistRC* extMeasure::computeR(uint len, double* valTable)
{
  extDistRC* rcUnit = nullptr;

  for (uint ii = 0; ii < _metRCTable.getCnt(); ii++) {
    extMetRCTable* rcModel = _metRCTable.get(ii);

    rcUnit = getOverRC(rcModel);
    if (rcUnit != nullptr) {
      _rc[ii]->res_ += rcUnit->res_ * len;
    }
  }
  return rcUnit;
}

extDistRC* extMeasure::computeUnderRC(uint len)
{
  extDistRC* rcUnit = nullptr;

  for (uint ii = 0; ii < _metRCTable.getCnt(); ii++) {
    extMetRCTable* rcModel = _metRCTable.get(ii);

    rcUnit = getUnderRC(rcModel);

    addRC(rcUnit, len, ii);
  }
  return rcUnit;
}

void extMeasure::printMets(FILE* fp)
{
  if (_overUnder) {
    fprintf(fp, "M%d over M%d under M%d ", _met, _underMet, _overMet);
  } else if (_over) {
    if (_diag) {
      fprintf(fp, "M%d over diag M%d ", _met, _underMet);
    } else {
      fprintf(fp, "M%d over M%d ", _met, _underMet);
    }
  } else {
    if (_diag) {
      fprintf(fp, "M%d under diag M%d ", _met, _overMet);
    } else {
      fprintf(fp, "M%d under M%d ", _met, _overMet);
    }
  }
}

void extMeasure::printStats(FILE* fp)
{
  fprintf(fp,
          "<==> w= %g[%g %g]  s= %g[%g]  th= %g[%g]  h= %g[%g]  r= %g",
          _w_m,
          _topWidth,
          _botWidth,
          _s_m,
          _seff,
          _t,
          _teff,
          _h,
          _heff,
          _r);
}

FILE* extRCModel::openFile(const char* topDir,
                           const char* name,
                           const char* suffix,
                           const char* permissions)
{
  char filename[2048];

  filename[0] = '\0';
  if (topDir != nullptr) {
    sprintf(filename, "%s/", topDir);
  }
  strcat(filename, name);
  if (suffix != nullptr) {
    strcat(filename, suffix);
  }

  FILE* fp = fopen(filename, permissions);
  if (fp == nullptr) {
    logger_->info(RCX,
                  486,
                  "Cannot open file {} with permissions {}",
                  filename,
                  permissions);
    return nullptr;
  }
  return fp;
}

void extRCModel::mkFileNames(extMeasure* m, char* wiresNameSuffix)
{
  char overUnder[128];

  if ((m->_overMet > 0) && (m->_underMet > 0)) {
    sprintf(overUnder, "M%doM%duM%d", m->_met, m->_underMet, m->_overMet);

  } else if (m->_overMet > 0) {
    if (m->_diag) {
      sprintf(overUnder, "M%dduM%d", m->_met, m->_overMet);
    } else {
      sprintf(overUnder, "M%duM%d", m->_met, m->_overMet);
    }

  } else if (m->_underMet >= 0) {
    sprintf(overUnder, "M%doM%d", m->_met, m->_underMet);

  } else {
    sprintf(overUnder, "Uknown");
  }

  double w = m->_w_m;
  double s = m->_s_m;
  double w2 = m->_w2_m;
  double s2 = m->_s2_m;

  sprintf(_wireDirName,
          "%s/%s/%s/W%g_W%g/S%g_S%g_L%d",
          _topDir,
          _patternName,
          overUnder,
          w,
          w2,
          s,
          s2,
          m->_len);

  if (wiresNameSuffix != nullptr) {
    sprintf(_wireFileName, "%s.%s", "wires", wiresNameSuffix);
  } else {
    sprintf(_wireFileName, "%s", "wires");
  }

  fprintf(_logFP, "PATTERN %s\n\n", _wireDirName);
  fflush(_logFP);
}

double get_nm(extMeasure* m, double n)
{
  if (n == 0) {
    return 0;
  }
  double a = 1000 * 1000.0 * n / m->_dbunit;
  return a;
}

int get_nm(int n, int units)
{
  if (n == 0) {
    return 0;
  }
  int a = (1000 * n) / units;
  return a;
}

void extRCModel::mkNet_prefix(extMeasure* m, const char* wiresNameSuffix)
{
  char overUnder[128];

  if ((m->_overMet > 0) && (m->_underMet > 0)) {
    sprintf(overUnder, "M%doM%duM%d", m->_met, m->_underMet, m->_overMet);

  } else if (m->_overMet > 0) {
    if (m->_diag) {
      sprintf(overUnder, "M%duuM%d", m->_met, m->_overMet);
    } else {
      sprintf(overUnder, "M%duM%d", m->_met, m->_overMet);
    }

  } else if (m->_underMet >= 0) {
    if (m->_diag) {
      sprintf(overUnder, "M%duuM%d", m->_underMet, m->_met);
    } else {
      sprintf(overUnder, "M%doM%d", m->_met, m->_underMet);
    }
  } else {
    sprintf(overUnder, "Unknown");
  }

  sprintf(_wireDirName,
          "%s_%s_W%gW%g_S%05dS%05d",
          _patternName,
          overUnder,
          get_nm(m, m->_w_m),
          get_nm(m, m->_w2_m),
          get_nm(m->_s_nm, m->_dbunit),
          get_nm(m->_s2_nm, m->_dbunit));

  if (wiresNameSuffix != nullptr) {
    sprintf(_wireFileName, "%s.%s", "wires", wiresNameSuffix);
  } else {
    sprintf(_wireFileName, "%s", "wires");
  }

  // fprintf(_logFP, "pattern Dir %s\n\n", _wireDirName);
  fflush(_logFP);
}

FILE* extRCModel::mkPatternFile()
{
  _parser->mkDirTree(_wireDirName, "/");

  FILE* fp = openFile(_wireDirName, _wireFileName, nullptr, "w");
  if (fp == nullptr) {
    return nullptr;
  }

  fprintf(fp, "PATTERN %s\n\n", _wireDirName);
  if (strcmp("TYP/Under3/M6uM7/W0.42_W0.42/S0.84_S0.84", _wireDirName) == 0) {
    fprintf(stdout, "%s\n", _wireDirName);
  }

  return fp;
}

FILE* extRCModel::openSolverFile()
{
  FILE* fp = openFile(_wireDirName, _wireFileName, ".out", "r");
  if (fp != nullptr) {
    _parser->setInputFP(fp);
  }

  return fp;
}

bool extRCModel::openCapLogFile()
{
  if (!_readSolver && !_runSolver) {
    _readCapLog = false;
    return true;
  }

  _readCapLog = false;
  if (_readSolver && !_runSolver) {
    _readCapLog = true;
  }

  const char* capLog = "caps.log";

  char buff[1024];
  sprintf(buff, "%s/%s", _topDir, _patternName);
  _parser->mkDirTree(buff, "/");

  FILE* fp = openFile(buff, capLog, nullptr, "r");

  if (fp == nullptr) {  // no previous run
    _capLogFP = openFile(buff, capLog, nullptr, "w");
    _parser->setInputFP(_capLogFP);
    return false;
  }
  fclose(fp);

  FILE* fp1 = nullptr;
  if (_readCapLog) {
    fp1 = openFile(buff, capLog, nullptr, "r");
    _capLogFP = openFile(buff, capLog, "out", "a");
  } else if (_metLevel > 0) {
    _capLogFP = openFile(buff, capLog, nullptr, "a");
  } else {
    try {
      std::filesystem::path path0(_topDir);
      path0 += _patternName;
      path0 += capLog;
      std::filesystem::path path1(_topDir);
      path0 += _patternName;
      path0 += std::string(capLog) + ".in";

      std::filesystem::rename(path0, path1);
    } catch (const std::filesystem::filesystem_error&) {
      logger_->error(
          RCX, 489, "mv failed: {}/{}/{}", _topDir, _patternName, capLog);
    }

    _capLogFP = openFile(buff, capLog, nullptr, "w");

    fp1 = openFile(buff, capLog, ".in", "r");
  }
  if (fp1 == nullptr) {
    return false;
  }

  _parser->setInputFP(fp1);

  _readCapLog = true;
  return true;
}

void extRCModel::closeCapLogFile()
{
  if (_readCapLog) {
    fclose(_capLogFP);
  }
}

void extRCModel::writeWires2(FILE* fp, extMeasure* measure, uint wireCnt)
{
  extMasterConductor* m = _process->getMasterConductor(measure->_met);
  double pitch = measure->_topWidth + measure->_seff;
  double min_pitch = 0.001 * (measure->_minWidth + measure->_minSpace);

  uint n = wireCnt / 2;  // ASSUME odd number of wires, 2 will also work
  double orig = 0.0;

  // assume origin = (0,0)
  double x = -min_pitch * (n - 1) - pitch - 0.5 * measure->_topWidth + orig;
  for (uint ii = 0; ii < n - 1; ii++) {
    m->writeRaphaelPoly(fp, ii + 1, 0.001 * measure->_minWidth, x, 0.0);
    x += min_pitch;
  }
  x += 0.5 * measure->_topWidth;
  m->writeRaphaelPoly(fp, n, x, 0.0);

  m->writeRaphaelPoly(fp, n + 1, orig, 1.0);
  x = 0.5 * measure->_w_m + measure->_s2_m;
  m->writeRaphaelPoly(fp, n + 2, measure->_w2_m, x, 0.0);

  x = orig + 0.5 * measure->_topWidth + measure->_w2_m + measure->_s2_m
      + 0.001 * measure->_minSpace;
  for (uint jj = n + 2; jj < wireCnt; jj++) {
    m->writeRaphaelPoly(fp, jj + 1, 0.001 * measure->_minWidth, x, 0.0);
    x += min_pitch;
  }
}

int extRCModel::writeBenchWires(FILE* fp, extMeasure* measure)
{
  uint grid_gap_cnt = 20;

  int bboxLL[2];
  bboxLL[measure->_dir] = measure->_ur[measure->_dir];
  bboxLL[!measure->_dir] = measure->_ll[!measure->_dir];

  extMasterConductor* m = _process->getMasterConductor(measure->_met);

  int n
      = measure->_wireCnt / 2;  // ASSUME odd number of wires, 2 will also work

  double pitchUp_print;
  pitchUp_print = measure->_topWidth + measure->_seff;
  double pitch_print = 0.001 * (measure->_minWidth + measure->_minSpace);

  uint w_layout = measure->_minWidth;
  uint s_layout = measure->_minSpace;

  double x = -(measure->_topWidth * 0.5 + pitchUp_print + pitch_print);

  measure->clean2dBoxTable(measure->_met, false);

  double x_tmp[50];
  uint netIdTable[50];
  uint idCnt = 1;
  int ii;
  for (ii = 0; ii < n - 1; ii++) {
    netIdTable[idCnt]
        = measure->createNetSingleWire(_wireDirName, idCnt, w_layout, s_layout);
    idCnt++;
    x_tmp[ii] = x;
    x -= pitch_print;
  }

  double X[50];
  ii--;
  int cnt = 0;
  for (; ii >= 0; ii--) {
    X[cnt++] = x_tmp[ii];
  }

  uint WW = measure->_w_nm;
  uint SS1;
  SS1 = measure->_s_nm;
  uint WW2 = measure->_w2_nm;
  uint SS2 = measure->_s2_nm;

  X[cnt++] = -pitchUp_print;
  int mid = cnt;
  netIdTable[idCnt]
      = measure->createNetSingleWire(_wireDirName, idCnt, WW, s_layout);
  idCnt++;

  X[cnt++] = 0.0;
  netIdTable[idCnt]
      = measure->createNetSingleWire(_wireDirName, idCnt, WW, SS1);
  idCnt++;
  uint base = measure->_ll[measure->_dir] + WW / 2;

  X[cnt++] = (SS2 + WW * 0.5) * 0.001;
  netIdTable[idCnt]
      = measure->createNetSingleWire(_wireDirName, idCnt, WW2, SS2);
  idCnt++;

  x = measure->_topWidth * 0.5 + 0.001 * (WW2 + SS2 + measure->_minSpace);
  for (int jj = 0; jj < n - 1; jj++) {
    X[cnt++] = x;
    x += pitch_print;
    netIdTable[idCnt]
        = measure->createNetSingleWire(_wireDirName, idCnt, w_layout, s_layout);
    idCnt++;
  }

  for (ii = 0; ii < cnt; ii++) {
    if (ii == mid) {
      m->writeRaphaelPoly(fp, netIdTable[ii + 1], X[ii], 1.0);
    } else if (ii == mid - 1) {
      m->writeRaphaelPoly(fp, netIdTable[ii + 1], X[ii], 0.0);
    } else if (ii == mid + 1) {
      m->writeRaphaelPoly(fp, netIdTable[ii + 1], 0.001 * WW2, X[ii], 0.0);
    } else {
      m->writeRaphaelPoly(
          fp, netIdTable[ii + 1], 0.001 * measure->_minWidth, X[ii], 0.0);
    }
  }

  if (measure->_diag) {
    int met = 0;
    if (measure->_overMet > 0) {
      met = measure->_overMet;
    } else if (measure->_underMet > 0) {
      met = measure->_underMet;
    }

    m = _process->getMasterConductor(met);
    double minWidth = _process->getConductor(met)->_min_width;
    double minSpace = _process->getConductor(met)->_min_spacing;
    double min_pitch = minWidth + minSpace;
    measure->clean2dBoxTable(met, false);
    int i;
    uint begin
        = base - lround(measure->_seff * 1000) + lround(minWidth * 1000) / 2;
    for (i = 0; i < n + 1; i++) {
      netIdTable[idCnt]
          = measure->createDiagNetSingleWire(_wireDirName,
                                             idCnt,
                                             begin,
                                             lround(1000 * minWidth),
                                             lround(1000 * minSpace),
                                             measure->_dir);
      begin -= lround(min_pitch * 1000);
      idCnt++;
    }

    measure->_ur[measure->_dir] += grid_gap_cnt * (w_layout + s_layout);

    fprintf(fp, "\nOPTIONS SET_GRID=10000;\n\n");
    fprintf(fp, "POTENTIAL");
    fprintf(fp, " \n");

    return cnt;
  }

  int bboxUR[2] = {measure->_ur[0], measure->_ur[1]};

  double pitchMult = 1.0;

  measure->clean2dBoxTable(measure->_underMet, true);
  measure->createContextNets(
      _wireDirName, bboxLL, bboxUR, measure->_underMet, pitchMult);

  measure->clean2dBoxTable(measure->_overMet, true);
  measure->createContextNets(
      _wireDirName, bboxLL, bboxUR, measure->_overMet, pitchMult);

  int main_xlo, main_ylo, main_xhi, main_yhi;
  measure->getBox(measure->_met, false, main_xlo, main_ylo, main_xhi, main_yhi);
  measure->_ur[measure->_dir] += grid_gap_cnt * (w_layout + s_layout);

  fprintf(fp, "\nOPTIONS SET_GRID=10000;\n\n");
  fprintf(fp, "POTENTIAL");
  fprintf(fp, " \n");

  return cnt;
}

void extRCModel::setOptions(const char* topDir,
                            const char* pattern,
                            bool writeFiles)
{
  _logFP = openFile("./", "rulesGen", ".log", "w");
  _filesFP = openFile("./", "patternFiles.", pattern, "w");
  // strcpy(_topDir, topDir);
  strcpy(_topDir, pattern);
  strcpy(_patternName, pattern);

  _writeFiles = true;
  _readSolver = true;
  _runSolver = true;
}

void extRCModel::setOptions(const char* topDir, const char* pattern)
{
  _logFP = openFile("./", "rulesGen", ".log", "w");
  strcpy(_topDir, topDir);
  strcpy(_patternName, pattern);

  _writeFiles = true;
  _readSolver = true;
  _runSolver = true;
}

void extRCModel::closeFiles()
{
  fflush(_logFP);

  if (_logFP != nullptr) {
    fclose(_logFP);
  }
  fflush(_filesFP);
  if (_filesFP != nullptr) {
    fclose(_filesFP);
  }
}

void extRCModel::cleanFiles()
{
  try {
    std::filesystem::remove_all(_wireDirName);
  } catch (const std::filesystem::filesystem_error& err) {
    logger_->error(RCX, 491, "rm failed on {}: {}", _wireDirName, err.what());
  }
}

int extRCModel::getOverUnderIndex(extMeasure* m, uint maxCnt)
{
  return getMetIndexOverUnder(
      m->_met, m->_underMet, m->_overMet, _layerCnt, maxCnt);
}

uint extRCModel::getUnderIndex(extMeasure* m)
{
  return m->_overMet - m->_met - 1;
}

void extDistWidthRCTable::addRCw(uint n, uint w, extDistRC* rc)
{
  int wIndex = _widthTable->findIndex(w);
  if (wIndex < 0) {
    wIndex = _widthTable->add(w);
  }

  _rcDistTable[n][wIndex]->addMeasureRC(rc);
}
/* DKF  10032024
void extMetRCTable::addRCw(extMeasure* m)
{
  extDistWidthRCTable* table = nullptr;
  int n;
  if (m->_overUnder) {
    n = extRCModel::getMetIndexOverUnder(m->_met,
                             m->_underMet,
                             m->_overMet,
                             _layerCnt,
                             _capOverUnder[m->_met]->_metCnt);
    assert(n >= 0);
    table = _capOverUnder[m->_met];
  } else if (m->_over) {
    n = m->_underMet;
    if (m->_res) {
      table = _resOver[m->_met];
    } else {
      table = _capOver[m->_met];
    }
  } else if (m->_diag) {
    n = m->getUnderIndex();
    if (m->_diagModel == 1) {  // TODO 0620 : diagModel=2
      table = _capDiagUnder[m->_met];
    }
  } else {
    n = m->getUnderIndex();
    table = _capUnder[m->_met];
  }
  if (table != nullptr) {
    table->addRCw(n, m->_w_nm, m->_tmpRC);
  }
}
*/
void extMetRCTable::addRCw(extMeasure* m)
{
  extDistWidthRCTable* table = nullptr;
  int n;
  if (m->_overUnder) {
    n = extRCModel::getMetIndexOverUnder(m->_met,
                                         m->_underMet,
                                         m->_overMet,
                                         _layerCnt,
                                         _capOverUnder[m->_met]->_metCnt);
    assert(n >= 0);
    if (m->_open) {
      table = _capOverUnder_open[m->_met][0];
    } else if (m->_over1) {
      table = _capOverUnder_open[m->_met][1];
    } else {
      table = _capOverUnder[m->_met];
    }
  } else if (m->_over) {
    n = m->_underMet;
    if (m->_res) {
      table = _resOver[m->_met];
    } else if (m->_open) {
      table = _capOver_open[m->_met][0];
    } else if (m->_over1) {
      table = _capOver_open[m->_met][1];
    } else {
      table = _capOver[m->_met];
    }
  } else if (m->_diag) {
    n = m->getUnderIndex();
    if (m->_diagModel == 1) {  // TODO 0620 : diagModel=2
      table = _capDiagUnder[m->_met];
    }
  } else {
    n = m->getUnderIndex();
    if (m->_open) {
      table = _capUnder_open[m->_met][0];
    } else if (m->_over1) {
      table = _capUnder_open[m->_met][1];
    } else {
      table = _capUnder[m->_met];
    }
  }
  if (table != nullptr) {
    table->addRCw(n, m->_w_nm, m->_tmpRC);
  }
}

void extRCModel::addRC(extMeasure* m)
{
  if (m->_overUnder) {
    uint maxCnt = _modelTable[m->_rIndex]->_capOverUnder[m->_met]->_metCnt;
    uint n = getOverUnderIndex(m, maxCnt);
    _modelTable[m->_rIndex]
        ->_capOverUnder[m->_met]
        ->_rcDistTable[n][m->_wIndex]
        ->addMeasureRC(m->_tmpRC);
  } else if (m->_over) {
    if (m->_res) {
      _modelTable[m->_rIndex]
          ->_resOver[m->_met]
          ->_rcDistTable[m->_underMet][m->_wIndex]
          ->addMeasureRC(m->_tmpRC);
    } else {
      _modelTable[m->_rIndex]
          ->_capOver[m->_met]
          ->_rcDistTable[m->_underMet][m->_wIndex]
          ->addMeasureRC(m->_tmpRC);
    }
  } else if (m->_diag) {
    uint n = getUnderIndex(m);
    if (_diagModel == 2) {
      _modelTable[m->_rIndex]
          ->_capDiagUnder[m->_met]
          ->_rcDiagDistTable[n][m->_wIndex][m->_dwIndex][m->_dsIndex]
          ->addMeasureRC(m->_tmpRC);
    } else {
      _modelTable[m->_rIndex]
          ->_capDiagUnder[m->_met]
          ->_rcDistTable[n][m->_wIndex]
          ->addMeasureRC(m->_tmpRC);
    }
  } else {
    uint n = getUnderIndex(m);
    _modelTable[m->_rIndex]
        ->_capUnder[m->_met]
        ->_rcDistTable[n][m->_wIndex]
        ->addMeasureRC(m->_tmpRC);
  }
}

void extMetRCTable::mkWidthAndSpaceMappings()
{
  for (uint ii = 1; ii < _layerCnt; ii++) {
    if (_capOver[ii] != nullptr) {
      _capOver[ii]->makeWSmapping();
    } else {
      logger_->info(RCX, 72, "Can't find <OVER> rules for {}", ii);
    }

    if (_resOver[ii] != nullptr) {
      _resOver[ii]->makeWSmapping();
    } else {
      logger_->info(RCX, 358, "Can't find <RESOVER> Res rules for {}", ii);
    }

    if (_capUnder[ii] != nullptr) {
      _capUnder[ii]->makeWSmapping();
    } else {
      logger_->info(RCX, 216, "Can't find <UNDER> rules for {}", ii);
    }

    if (_capOverUnder[ii] != nullptr) {
      _capOverUnder[ii]->makeWSmapping();
    } else if ((ii > 1) && (ii < _layerCnt - 1)) {
      logger_->info(RCX, 217, "Can't find <OVERUNDER> rules for {}", ii);
    }
  }
}

void extRCModel::writeRules(char* name, bool binary)
{
  bool writeRes = true;
  FILE* fp = fopen(name, "w");

  fprintf(fp, "Extraction Rules for OpenRCX\n\n");
  if (_diag || _diagModel > 0) {
    if (_diagModel == 1) {
      fprintf(fp, "DIAGMODEL ON\n\n");
    } else if (_diagModel == 2) {
      fprintf(fp, "DIAGMODEL TRUE\n\n");
    }
  }

  fprintf(fp, "LayerCount %d\n", _layerCnt - 1);
  fprintf(fp, "DensityRate %d ", _modelCnt);
  for (uint kk = 0; kk < _modelCnt; kk++) {
    fprintf(fp, " %g", _dataRateTable->get(kk));
  }
  fprintf(fp, "\n");

  for (uint m = 0; m < _modelCnt; m++) {
    fprintf(fp, "\nDensityModel %d\n", m);

    for (uint ii = 1; ii < _layerCnt; ii++) {
      if (writeRes) {
        if (_modelTable[m]->_resOver[ii] != nullptr) {
          _modelTable[m]->_resOver[ii]->writeRulesOver_res(fp, binary);
        } else if ((m > 0) && (_modelTable[0]->_resOver[ii] != nullptr)) {
          _modelTable[0]->_resOver[ii]->writeRulesOver_res(fp, binary);
        } else if (m == 0) {
          logger_->info(RCX,
                        413,
                        "Cannot write <OVER> Res rules for <DensityModel> {} "
                        "and layer {}",
                        m,
                        ii);
        }
      }
      if (_modelTable[m]->_capOver[ii] != nullptr) {
        _modelTable[m]->_capOver[ii]->writeRulesOver(fp, binary);
      } else if ((m > 0) && (_modelTable[0]->_capOver[ii] != nullptr)) {
        _modelTable[0]->_capOver[ii]->writeRulesOver(fp, binary);
      } else if (m == 0) {
        logger_->info(
            RCX,
            412,
            "Cannot write <OVER> rules for <DensityModel> {} and layer {}",
            m,
            ii);
      }

      if (_modelTable[m]->_capUnder[ii] != nullptr) {
        _modelTable[m]->_capUnder[ii]->writeRulesUnder(fp, binary);
      } else if ((m > 0) && (_modelTable[0]->_capUnder[ii] != nullptr)) {
        _modelTable[0]->_capUnder[ii]->writeRulesUnder(fp, binary);
      } else if (m == 0) {
        logger_->info(
            RCX,
            219,
            "Cannot write <UNDER> rules for <DensityModel> {} and layer {}",
            m,
            ii);
      }

      if (_modelTable[m]->_capDiagUnder[ii] != nullptr) {
        if (_diagModel == 1) {
          _modelTable[m]->_capDiagUnder[ii]->writeRulesDiagUnder(fp, binary);
        }
        if (_diagModel == 2) {
          _modelTable[m]->_capDiagUnder[ii]->writeRulesDiagUnder2(fp, binary);
        }
      } else if ((m > 0) && (_modelTable[0]->_capDiagUnder[ii] != nullptr)) {
        if (_diagModel == 1) {
          _modelTable[0]->_capDiagUnder[ii]->writeRulesDiagUnder(fp, binary);
        }
        if (_diagModel == 2) {
          _modelTable[0]->_capDiagUnder[ii]->writeRulesDiagUnder2(fp, binary);
        }
      } else if (m == 0) {
        logger_->info(
            RCX,
            249,
            "Cannot write <DIAGUNDER> rules for <DensityModel> {} and layer {}",
            m,
            ii);
      }

      if (_modelTable[m]->_capOverUnder[ii] != nullptr) {
        _modelTable[m]->_capOverUnder[ii]->writeRulesOverUnder(fp, binary);
      } else if ((m > 0) && (_modelTable[0]->_capOverUnder[ii] != nullptr)) {
        _modelTable[0]->_capOverUnder[ii]->writeRulesOverUnder(fp, binary);
      } else if ((m == 0) && (ii > 1) && (ii < _layerCnt - 1)) {
        logger_->info(
            RCX,
            221,
            "Cannot write <OVERUNDER> rules for <DensityModel> {} and layer {}",
            m,
            ii);
      }
    }
    fprintf(fp, "END DensityModel %d\n", m);
  }
  fclose(fp);
}

uint extRCModel::readMetalHeader(Ath__parser* parser,
                                 uint& met,
                                 const char* keyword,
                                 bool bin,
                                 bool ignore)
{
  if (parser->isKeyword(0, "END")
      && (strcmp(parser->get(1), "DensityModel") == 0)) {
    return 0;
  }

  if (!(parser->parseNextLine() > 0)) {
    return 0;
  }

  if (parser->isKeyword(0, "Metal") && (strcmp(parser->get(2), keyword) == 0)) {
    met = parser->getInt(1);
    return 1;
  }

  return 0;
}

/*
void extMetRCTable::allocateInitialTables(uint layerCnt,
                                          uint widthCnt,
                                          bool over,
                                          bool under,
                                          bool diag)
{
  for (uint met = 1; met < _layerCnt; met++) {
    if (over && under && (met > 1) && (met < _layerCnt - 1)) {
      int n = extRCModel::getMaxMetIndexOverUnder(met, layerCnt);
      _capOverUnder[met] = new extDistWidthRCTable(
          false, met, layerCnt, n + 1, widthCnt, _rcPoolPtr);
    }
    if (over) {
      _capOver[met] = new extDistWidthRCTable(
          true, met, layerCnt, met, widthCnt, _rcPoolPtr);
      _resOver[met] = new extDistWidthRCTable(
          true, met, layerCnt, met, widthCnt, _rcPoolPtr);
    }
    if (under) {
      _capUnder[met] = new extDistWidthRCTable(
          false, met, layerCnt, _layerCnt - met - 1, widthCnt, _rcPoolPtr);
    }
    if (diag) {
      _capDiagUnder[met] = new extDistWidthRCTable(
          false, met, layerCnt, _layerCnt - met - 1, widthCnt, _rcPoolPtr);
    }
  }
}
*/
void extMetRCTable::allocateInitialTables(uint widthCnt,
                                          bool over,
                                          bool under,
                                          bool diag)
{
  for (uint met = 1; met < _layerCnt; met++) {
    if (over && under && (met > 1) && (met < _layerCnt - 1)) {
      int n = extRCModel::getMaxMetIndexOverUnder(met, _layerCnt);
      _capOverUnder[met] = new extDistWidthRCTable(
          false, met, _layerCnt, n + 1, widthCnt, _rcPoolPtr, _OUREVERSEORDER);
      for (uint jj = 0; jj < _wireCnt; jj++) {
        _capOverUnder_open[met][jj] = new extDistWidthRCTable(false,
                                                              met,
                                                              _layerCnt,
                                                              n + 1,
                                                              widthCnt,
                                                              _rcPoolPtr,
                                                              _OUREVERSEORDER);
      }
    }
    if (over) {
      _capOver[met] = new extDistWidthRCTable(
          true, met, _layerCnt, met, widthCnt, _rcPoolPtr, _OUREVERSEORDER);
      _resOver[met] = new extDistWidthRCTable(
          true, met, _layerCnt, met, widthCnt, _rcPoolPtr, _OUREVERSEORDER);
      for (uint jj = 0; jj < _wireCnt; jj++) {
        _capOver_open[met][jj] = new extDistWidthRCTable(
            true, met, _layerCnt, met, widthCnt, _rcPoolPtr, _OUREVERSEORDER);
      }
    }
    if (under) {
      _capUnder[met] = new extDistWidthRCTable(false,
                                               met,
                                               _layerCnt,
                                               _layerCnt - met - 1,
                                               widthCnt,
                                               _rcPoolPtr,
                                               _OUREVERSEORDER);
      for (uint jj = 0; jj < _wireCnt; jj++) {
        _capUnder_open[met][jj] = new extDistWidthRCTable(false,
                                                          met,
                                                          _layerCnt,
                                                          _layerCnt - met - 1,
                                                          widthCnt,
                                                          _rcPoolPtr,
                                                          _OUREVERSEORDER);
      }
    }
    if (diag) {
      _capDiagUnder[met] = new extDistWidthRCTable(false,
                                                   met,
                                                   _layerCnt,
                                                   _layerCnt - met - 1,
                                                   widthCnt,
                                                   _rcPoolPtr,
                                                   _OUREVERSEORDER);
    }
  }
}

Ath__array1D<double>* extRCModel::readHeaderAndWidth(Ath__parser* parser,
                                                     uint& met,
                                                     const char* ouKey,
                                                     const char* wKey,
                                                     bool bin,
                                                     bool ignore)
{
  if (readMetalHeader(parser, met, ouKey, bin, ignore) <= 0) {
    return nullptr;
  }

  if (!(parser->parseNextLine() > 0)) {
    return nullptr;
  }

  return parser->readDoubleArray("WIDTH", 4);
}

uint extRCModel::readRules(Ath__parser* parser,
                           uint m,
                           uint ii,
                           const char* ouKey,
                           const char* wKey,
                           bool over,
                           bool under,
                           bool bin,
                           bool diag,
                           bool ignore,
                           double dbFactor)
{
  uint cnt = 0;
  uint met = 0;
  Ath__array1D<double>* wTable
      = readHeaderAndWidth(parser, met, ouKey, wKey, bin, false);

  if (wTable == nullptr) {
    return 0;
  }

  uint widthCnt = wTable->getCnt();

  extDistWidthRCTable* dummy = nullptr;
  if (ignore) {
    dummy = new extDistWidthRCTable(
        true, met, _layerCnt, widthCnt, _OUREVERSEORDER);
  }

  uint diagWidthCnt = 0;
  uint diagDistCnt = 0;

  if (diag && strcmp(ouKey, "DIAGUNDER") == 0 && _diagModel == 2) {
    parser->parseNextLine();
    if (parser->isKeyword(0, "DIAG_WIDTH")) {
      diagWidthCnt = parser->getInt(3);
    }
    parser->parseNextLine();
    if (parser->isKeyword(0, "DIAG_DIST")) {
      diagDistCnt = parser->getInt(3);
    }
  }

  if (over && under && (met > 1)) {
    if (!ignore) {
      _modelTable[m]->allocOverUnderTable(met, wTable, dbFactor);
      _modelTable[m]->_capOverUnder[met]->readRulesOverUnder(
          parser, widthCnt, bin, ignore, dbFactor);
    } else {
      dummy->readRulesOverUnder(parser, widthCnt, bin, ignore, dbFactor);
    }
  } else if (over) {
    if (strcmp(ouKey, "OVER") == 0) {
      if (!ignore) {
        _modelTable[m]->_capOver[met]->readRulesOver(
            parser, widthCnt, bin, ignore, "OVER", dbFactor);
      } else {
        dummy->readRulesOver(parser, widthCnt, bin, ignore, "OVER", dbFactor);
      }
    } else {  // RESOVER
      if (!ignore) {
        _modelTable[m]->allocOverTable(met, wTable, dbFactor);
        _modelTable[m]->_resOver[met]->readRulesOver(
            parser, widthCnt, bin, ignore, "RESOVER", dbFactor);
      } else {
        dummy->readRulesOver(
            parser, widthCnt, bin, ignore, "RESOVER", dbFactor);
      }
    }
  } else if (under) {
    if (!ignore) {
      _modelTable[m]->allocUnderTable(met, wTable, dbFactor);
      _modelTable[m]->_capUnder[met]->readRulesUnder(
          parser, widthCnt, bin, ignore, dbFactor);
    } else {
      dummy->readRulesUnder(parser, widthCnt, bin, ignore, dbFactor);
    }
  } else if (diag) {
    if (!ignore && _diagModel == 2) {
      _modelTable[m]->allocDiagUnderTable(
          met, wTable, diagWidthCnt, diagDistCnt, dbFactor);
      _modelTable[m]->_capDiagUnder[met]->readRulesDiagUnder(
          parser, widthCnt, diagWidthCnt, diagDistCnt, bin, ignore, dbFactor);
    } else if (!ignore && _diagModel == 1) {
      _modelTable[m]->allocDiagUnderTable(met, wTable, dbFactor);
      _modelTable[m]->_capDiagUnder[met]->readRulesDiagUnder(
          parser, widthCnt, bin, ignore, dbFactor);
    } else if (ignore) {
      if (_diagModel == 2) {
        dummy->readRulesDiagUnder(
            parser, widthCnt, diagWidthCnt, diagDistCnt, bin, ignore, dbFactor);
      } else if (_diagModel == 1) {
        dummy->readRulesDiagUnder(parser, widthCnt, bin, ignore, dbFactor);
      }
    }
  }
  if (ignore) {
    delete dummy;
  }

  delete wTable;

  return cnt;
}

bool extRCModel::readRules_v1(char* name,
                              bool bin,
                              bool over,
                              bool under,
                              bool overUnder,
                              bool diag,
                              uint cornerCnt,
                              const uint* cornerTable,
                              double dbFactor)
{
  _OUREVERSEORDER = false;
  diag = false;
  free(_ruleFileName);
  _ruleFileName = strdup(name);
  Ath__parser parser(logger_);
  parser.addSeparator("\r");
  parser.openFile(name);
  while (parser.parseNextLine() > 0) {
    if (parser.isKeyword(0, "OUREVERSEORDER")) {
      if (strcmp(parser.get(1), "ON") == 0) {
        _OUREVERSEORDER = true;
      }
    }
    if (parser.isKeyword(0, "DIAGMODEL")) {
      if (strcmp(parser.get(1), "ON") == 0) {
        _diagModel = 1;
        diag = true;
      } else if (strcmp(parser.get(1), "TRUE") == 0) {
        _diagModel = 2;
        diag = true;
      }
      continue;
    }

    if (parser.isKeyword(0, "rcStats")) {
      _layerCnt = parser.getInt(2);
      createModelTable(1, _layerCnt);
      for (uint kk = 0; kk < _modelCnt; kk++) {
        _dataRateTable->add(0.0);
      }

      _modelTable[0]->allocateInitialTables(10, true, true, true);

      _modelTable[0]->readRCstats(&parser);

      continue;
    }
    if (parser.isKeyword(0, "Layer")) {
      _layerCnt = parser.getInt(2);
      continue;
    }
    if (parser.isKeyword(0, "LayerCount")) {
      _layerCnt = parser.getInt(1) + 1;
      _verticalDiag = true;
      continue;
    }
    if (parser.isKeyword(0, "DensityRate")) {
      uint rulesFileModelCnt = parser.getInt(1);
      if (cornerCnt > 0) {
        if ((rulesFileModelCnt > 0) && (rulesFileModelCnt < cornerCnt)) {
          logger_->warn(
              RCX,
              222,
              "There were {} extraction models defined but only {} exists "
              "in the extraction rules file {}",
              cornerCnt,
              rulesFileModelCnt,
              name);
          return false;
        }
        createModelTable(cornerCnt, _layerCnt);

        for (uint jj = 0; jj < cornerCnt; jj++) {
          uint modelIndex = cornerTable[jj];

          uint kk;
          for (kk = 0; kk < rulesFileModelCnt; kk++) {
            if (modelIndex != kk) {
              continue;
            }
            _dataRateTable->add(parser.getDouble(kk + 2));
            break;
          }
          if (kk == rulesFileModelCnt) {
            logger_->warn(RCX,
                          223,
                          "Cannot find model index {} in extRules file {}",
                          modelIndex,
                          name);
            return false;
          }
        }
      } else {
        createModelTable(1, _layerCnt);

        for (uint kk = 0; kk < _modelCnt; kk++) {
          _dataRateTable->add(parser.getDouble(kk + 2));
        }
        for (uint ii = 0; ii < _modelCnt; ii++) {
          _modelTable[ii]->_rate = _dataRateTable->get(ii);
        }
      }
      continue;
    }

    if (parser.isKeyword(0, "DensityModel")) {
      uint m = parser.getInt(1);
      uint modelIndex = m;
      bool skipModel = false;
      if (cornerCnt > 0) {
        uint jj = 0;
        for (; jj < cornerCnt; jj++) {
          if (m == cornerTable[jj]) {
            break;
          }
        }
        if (jj == cornerCnt) {
          skipModel = true;
          modelIndex = 0;
        } else {
          skipModel = false;
          modelIndex = jj;
        }
      } else {
        if (modelIndex) {
          skipModel = true;
        }
      }

      bool res_skipModel = false;

      for (uint ii = 1; ii < _layerCnt; ii++) {
        if (!res_skipModel) {
          readRules(&parser,
                    modelIndex,
                    ii,
                    "RESOVER",
                    "WIDTH",
                    over,
                    false,
                    bin,
                    false,
                    res_skipModel,
                    dbFactor);
        }
        readRules(&parser,
                  modelIndex,
                  ii,
                  "OVER",
                  "WIDTH",
                  over,
                  false,
                  bin,
                  false,
                  skipModel,
                  dbFactor);
        if (ii < _layerCnt - 1) {
          readRules(&parser,
                    modelIndex,
                    ii,
                    "UNDER",
                    "WIDTH",
                    false,
                    under,
                    bin,
                    false,
                    skipModel,
                    dbFactor);
          if (diag) {
            readRules(&parser,
                      modelIndex,
                      ii,
                      "DIAGUNDER",
                      "WIDTH",
                      false,
                      false,
                      bin,
                      diag,
                      skipModel,
                      dbFactor);
          }
        }

        if ((ii > 1) && (ii < _layerCnt - 1)) {
          readRules(&parser,
                    modelIndex,
                    ii,
                    "OVERUNDER",
                    "WIDTH",
                    overUnder,
                    overUnder,
                    bin,
                    false,
                    skipModel,
                    dbFactor);
        }
      }
      parser.parseNextLine();
    }
  }
  return true;
}

double extRCModel::measureResistance(extMeasure* m,
                                     double ro,
                                     double top_widthR,
                                     double bot_widthR,
                                     double thicknessR)
{
  double r = ro / ((top_widthR + bot_widthR) * thicknessR * 0.5);
  return r;
}
bool extRCModel::measurePatternVar(extMeasure* m,
                                   double top_width,
                                   double bot_width,
                                   double thickness,
                                   uint wireCnt,
                                   char* wiresNameSuffix,
                                   double res)
{
  if (m->_simVersion > 0) {
    return measurePatternVar_3D(
        m, top_width, bot_width, thickness, wireCnt, wiresNameSuffix, res);
  }

  m->setEffParams(top_width, bot_width, thickness);
  double thicknessChange
      = _process->adjustMasterLayersForHeight(m->_met, thickness);

  // _process->getMasterConductor(m->_met)->reset(m->_heff, top_width,
  // bot_width, thickness);
  _process->getMasterConductor(m->_met)->resetWidth(top_width, bot_width);

  mkFileNames(m, wiresNameSuffix);

  printCommentLine('$', m);
  fprintf(_logFP, "%s\n", _commentLine);
  fprintf(_logFP, "%c %g thicknessChange\n", '$', thicknessChange);
  fflush(_logFP);

  if (_writeFiles) {
    FILE* wfp = mkPatternFile();

    if (wfp == nullptr) {
      return false;  // should be an exception!! and return!
    }

    _process->adjustMasterDielectricsForHeight(m->_met, thicknessChange);

    if (_commentFlag) {
      fprintf(wfp, "%s\n", _commentLine);
    }

    if (m->_benchFlag) {
      writeBenchWires(wfp, m);
    } else {
      fprintf(_filesFP, "%s/wires\n", _wireDirName);
    }

    fclose(wfp);
  }
  return true;
}
void extRCModel::printCommentLine(char commentChar, extMeasure* m)
{
  sprintf(_commentLine,
          "%c %s w= %g s= %g r= %g\n\n%c %s %6.3f %6.3f top_width\n%c %s %6.3f "
          "%g bot_width\n%c %s %6.3f %6.3f spacing\n%c %s %6.3f %6.3f height "
          "\n%c %s %6.3f %6.3f thicknes\n",
          commentChar,
          "Layout params",
          m->_w_m,
          m->_s_m,
          m->_r,
          commentChar,
          "Layout/Eff",
          m->_w_m,
          m->_topWidth,
          commentChar,
          "Layout/Eff",
          m->_w_m,
          m->_botWidth,
          commentChar,
          "Layout/Eff",
          m->_s_m,
          m->_seff,
          commentChar,
          "Layout/Eff",
          m->_h,
          m->_heff,
          commentChar,
          "Layout/Eff",
          m->_t,
          m->_teff);
  _commentFlag = true;
}

void extRCModel::getDiagTables(extMeasure* m, uint widthCnt, uint spaceCnt)
{
  Ath__array1D<double>* diagSTable0 = nullptr;
  Ath__array1D<double>* diagWTable0 = nullptr;
  diagSTable0 = _process->getDiagSpaceTable(m->_overMet);
  diagWTable0 = _process->getWidthTable(m->_overMet);
  m->_diagWidthTable0.resetCnt();
  if (diagWTable0) {
    for (uint wIndex = 0;
         (wIndex < diagWTable0->getCnt()) && (wIndex < widthCnt);
         wIndex++) {
      double w = diagWTable0->get(wIndex);
      m->_diagWidthTable0.add(w);
    }
  }
  m->_diagSpaceTable0.resetCnt();
  if (diagSTable0) {
    for (uint dsIndex = 0;
         (dsIndex < diagSTable0->getCnt()) && (dsIndex < spaceCnt);
         dsIndex++) {
      double ds = diagSTable0->get(dsIndex);
      m->_diagSpaceTable0.add(ds);
    }
  }
}

void extRCModel::computeTables(extMeasure* m,
                               uint wireCnt,
                               uint widthCnt,
                               uint spaceCnt,
                               uint dCnt)
{
  extVariation* xvar = nullptr;
  if (!_maxMinFlag) {
    xvar = _process->getVariation(m->_met);
  }

  m->_thickVarFlag = _process->getThickVarFlag();

  Ath__array1D<double>* wTable = nullptr;
  Ath__array1D<double>* sTable = nullptr;
  Ath__array1D<double>* dTable = nullptr;
  Ath__array1D<double>* pTable = nullptr;
  Ath__array1D<double>* wTable0 = nullptr;
  Ath__array1D<double>* sTable0 = nullptr;
  Ath__array1D<double>* diagSTable0 = nullptr;
  Ath__array1D<double>* diagWTable0 = nullptr;
  if (xvar != nullptr) {
    wTable = xvar->getWidthTable();
    sTable = xvar->getSpaceTable();
    dTable = xvar->getDataRateTable();
    pTable = xvar->getPTable();
    wTable0 = _process->getWidthTable(m->_met);
    sTable0 = _process->getSpaceTable(m->_met);
    if (_diagModel == 2 && m->_overMet < (int) _layerCnt) {
      diagSTable0 = _process->getDiagSpaceTable(m->_overMet);
      diagWTable0 = _process->getWidthTable(m->_overMet);
    } else {
      diagSTable0 = _process->getDiagSpaceTable(m->_met);
    }
  } else {  // no variation
    wTable = _process->getWidthTable(m->_met);
    sTable = _process->getSpaceTable(m->_met);
    dTable = _process->getDataRateTable(m->_met);
    wTable0 = _process->getWidthTable(m->_met);
    sTable0 = _process->getSpaceTable(m->_met);
    if (_diagModel == 2 && m->_overMet < (int) _layerCnt) {
      diagSTable0 = _process->getDiagSpaceTable(m->_overMet);
      diagWTable0 = _process->getWidthTable(m->_overMet);
    } else {
      diagSTable0 = _process->getDiagSpaceTable(m->_met);
    }
    if (_maxMinFlag) {
      for (uint i = 1; i < 3; i++) {
        dTable->add(i);
      }
    }
  }
  m->_widthTable.resetCnt();
  for (uint wIndex = 0; (wIndex < wTable->getCnt()) && (wIndex < widthCnt);
       wIndex++) {
    double w = wTable->get(wIndex);  // layout
    m->_widthTable.add(w);
  }
  if (_diagModel == 2 && m->_overMet < (int) _layerCnt) {
    m->_diagWidthTable0.resetCnt();
    if (diagWTable0) {
      for (uint wIndex = 0;
           (wIndex < diagWTable0->getCnt()) && (wIndex < widthCnt);
           wIndex++) {
        double w = diagWTable0->get(wIndex);
        m->_diagWidthTable0.add(w);
      }
    }
  }
  m->_spaceTable.resetCnt();
  if (m->_diagModel == 1) {
    m->_spaceTable.add(0.0);
  }
  for (uint sIndex = 0; (sIndex < sTable->getCnt()) && (sIndex < spaceCnt);
       sIndex++) {
    double s = sTable->get(sIndex);  // layout
    m->_spaceTable.add(s);
  }
  if (m->_diagModel == 2) {
    m->_spaceTable.add(99);
  }
  double lastSpacing = 2 * sTable->getLast();
  m->_spaceTable.add(lastSpacing);
  // m->_spaceTable.add(100); // DKF 7/27/24
  m->_diagSpaceTable0.resetCnt();
  if (diagSTable0) {
    for (uint dsIndex = 0;
         (dsIndex < diagSTable0->getCnt()) && (dsIndex < spaceCnt);
         dsIndex++) {
      double ds = diagSTable0->get(dsIndex);
      m->_diagSpaceTable0.add(ds);
    }
  }
  m->_dataTable.resetCnt();
  if (!_maxMinFlag && xvar != nullptr) {
    m->_dataTable.add(0.0);
  }
  m->_widthTable0.resetCnt();
  for (uint wIndex1 = 0; (wIndex1 < wTable0->getCnt()) && (wIndex1 < widthCnt);
       wIndex1++) {
    double w = wTable0->get(wIndex1);
    m->_widthTable0.add(w);
  }
  m->_spaceTable0.resetCnt();
  if (m->_diagModel == 1) {
    m->_spaceTable0.add(0.0);
  }
  for (uint sIndex1 = 0; (sIndex1 < sTable0->getCnt()) && (sIndex1 < spaceCnt);
       sIndex1++) {
    double s = sTable0->get(sIndex1);
    m->_spaceTable0.add(s);
  }
  if (m->_diagModel == 2) {
    m->_spaceTable0.add(99);
  }
  bool add_last_spacing = false;
  if (add_last_spacing) {
    lastSpacing = 2 * m->_spaceTable0.getLast();
    m->_spaceTable0.add(lastSpacing);
    // DKF 7/27/24 m->_spaceTable0.add(100);
  }

  for (uint dIndex = 0; (dIndex < dTable->getCnt()) && (dIndex < dCnt);
       dIndex++) {
    double r = dTable->get(dIndex);  // layout
    m->_dataTable.add(r);
  }
  if (pTable != nullptr) {
    m->_pTable.resetCnt();
    for (uint pIndex = 0; pIndex < pTable->getCnt(); pIndex++) {
      double p = pTable->get(pIndex);
      m->_pTable.add(p);
    }
  }
}

void extRCModel::allocOverTable(extMeasure* measure)
{
  for (uint ii = 0; ii < measure->_dataTable.getCnt(); ii++) {
    if (!ii) {
      _modelTable[ii]->allocOverTable(measure->_met, &measure->_widthTable0);
    } else {
      _modelTable[ii]->allocOverTable(measure->_met, &measure->_widthTable);
    }
  }
}

void extRCModel::allocDiagUnderTable(extMeasure* measure)
{
  for (uint ii = 0; ii < measure->_dataTable.getCnt(); ii++) {
    if (!ii) {
      if (_diagModel == 2) {
        _modelTable[ii]->allocDiagUnderTable(
            measure->_met,
            &measure->_widthTable0,
            measure->_diagWidthTable0.getCnt(),
            measure->_diagSpaceTable0.getCnt());
      } else if (_diagModel == 1) {
        _modelTable[ii]->allocDiagUnderTable(measure->_met,
                                             &measure->_widthTable0);
      }
    } else {
      if (_diagModel == 2) {
        _modelTable[ii]->allocDiagUnderTable(
            measure->_met,
            &measure->_widthTable,
            measure->_diagWidthTable0.getCnt(),
            measure->_diagSpaceTable0.getCnt());
      } else if (_diagModel == 1) {
        _modelTable[ii]->allocDiagUnderTable(measure->_met,
                                             &measure->_widthTable);
      }
    }
  }
}

void extRCModel::setDiagUnderTables(extMeasure* measure)
{
  for (uint ii = 0; ii < measure->_dataTable.getCnt(); ii++) {
    _modelTable[ii]->setDiagUnderTables(measure->_met,
                                        measure->_overMet,
                                        &measure->_diagWidthTable0,
                                        &measure->_diagSpaceTable0);
  }
}

void extRCModel::allocUnderTable(extMeasure* measure)
{
  for (uint ii = 0; ii < measure->_dataTable.getCnt(); ii++) {
    if (!ii) {
      _modelTable[ii]->allocUnderTable(measure->_met, &measure->_widthTable0);
    } else {
      _modelTable[ii]->allocUnderTable(measure->_met, &measure->_widthTable);
    }
  }
}

void extRCModel::allocOverUnderTable(extMeasure* measure)
{
  for (uint ii = 0; ii < measure->_dataTable.getCnt(); ii++) {
    if (!ii) {
      _modelTable[ii]->allocOverUnderTable(measure->_met,
                                           &measure->_widthTable0);
    } else {
      _modelTable[ii]->allocOverUnderTable(measure->_met,
                                           &measure->_widthTable);
    }
  }
}

uint extMain::writeRules(const char* name, const char* rulesFile)
{
  GenExtRules(rulesFile);
  return 0;
}

uint extRCModel::findBiggestDatarateIndex(double d)
{
  return _dataRateTable->findNextBiggestIndex(d, 1);
}

int extRCModel::findDatarateIndex(double d)
{
  for (uint ii = 0; ii < _modelCnt; ii++) {
    if (d == _dataRateTable->get(ii)) {
      return ii;
    }
    if (d > _dataRateTable->get(ii)) {
      return ii - 1;
    }
  }
  return -1;
}

extDistWidthRCTable* extRCModel::getWidthDistRCtable(uint met,
                                                     int mUnder,
                                                     int mOver,
                                                     int& n,
                                                     double dRate)
{
  int rIndex = 0;
  if (dRate > 0) {
    rIndex = findDatarateIndex(dRate);
    if (rIndex < 0) {
      return nullptr;
    }
  }
  if ((mUnder > 0) && (mOver > 0)) {
    n = getMetIndexOverUnder(met,
                             mUnder,
                             mOver,
                             _layerCnt,
                             _modelTable[rIndex]->_capOverUnder[met]->_metCnt);
    assert(n >= 0);
    return _modelTable[rIndex]->_capOverUnder[met];
  }
  if (mOver) {
    n = mUnder;
    return _modelTable[rIndex]->_capOver[met];
  }
  n = mOver - met - 1;
  return _modelTable[rIndex]->_capUnder[met];
}

}  // namespace rcx
