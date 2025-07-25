// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "db/infra/frPoint.h"
#include "db/obj/frBlockObject.h"
#include "frShape.h"

namespace drt {
class frViaDef;
class frPinAccess;
class frMTerm;
class frAccessPoint : public frBlockObject
{
 public:
  // constructors
  frAccessPoint(const Point& point, frLayerNum layerNum)
      : point_(point), layerNum_(layerNum)
  {
  }
  frAccessPoint() = default;
  frAccessPoint(const frAccessPoint& rhs)
      : frBlockObject(rhs),
        point_(rhs.point_),
        layerNum_(rhs.layerNum_),
        accesses_(rhs.accesses_),
        viaDefs_(rhs.viaDefs_),
        typeL_(rhs.typeL_),
        typeH_(rhs.typeH_),
        pathSegs_(rhs.pathSegs_)
  {
  }
  frAccessPoint& operator=(const frAccessPoint&) = delete;
  // getters
  const Point& getPoint() const { return point_; }
  frLayerNum getLayerNum() const { return layerNum_; }
  bool hasAccess() const
  {
    return (accesses_[0] || accesses_[1] || accesses_[2] || accesses_[3]
            || accesses_[4] || accesses_[5]);
  }
  bool hasAccess(const frDirEnum& dir) const
  {
    switch (dir) {
      case (frDirEnum::E):
        return accesses_[0];
        break;
      case (frDirEnum::S):
        return accesses_[1];
        break;
      case (frDirEnum::W):
        return accesses_[2];
        break;
      case (frDirEnum::N):
        return accesses_[3];
        break;
      case (frDirEnum::U):
        return accesses_[4];
        break;
      case (frDirEnum::D):
        return accesses_[5];
        break;
      default:
        return false;
    }
  }
  bool hasHorzAccess() const { return accesses_[0] || accesses_[2]; }
  bool hasVertAccess() const { return accesses_[1] || accesses_[3]; }
  bool hasViaAccess() const { return accesses_[4] || accesses_[5]; }
  bool hasPlanarAccess() const { return hasVertAccess() || hasHorzAccess(); }
  const std::vector<bool>& getAccess() const { return accesses_; }
  bool hasViaDef(int numCut = 1, int idx = 0) const
  {
    // first check numCuts
    int numCutIdx = numCut - 1;
    if (numCutIdx >= 0 && numCutIdx < (int) viaDefs_.size()) {
      ;
    } else {
      return false;
    }
    // then check idx
    if (idx >= 0 && idx < (int) (viaDefs_[numCutIdx].size())) {
      return true;
    }
    return false;
  }
  // e.g., getViaDefs()     --> get all one-cut viadefs
  // e.g., getViaDefs(1)    --> get all one-cut viadefs
  // e.g., getViaDefs(2)    --> get all two-cut viadefs
  const std::vector<const frViaDef*>& getViaDefs(int numCut = 1) const
  {
    return viaDefs_[numCut - 1];
  }
  std::vector<const frViaDef*>& getViaDefs(int numCut = 1)
  {
    return viaDefs_[numCut - 1];
  }
  const std::vector<std::vector<const frViaDef*>>& getAllViaDefs() const
  {
    return viaDefs_;
  }
  void sortViaDefs(const std::map<const frViaDef*, int> cost_map)
  {
    auto cmp = [&](const frViaDef* a, const frViaDef* b) {
      return cost_map.at(a) < cost_map.at(b);
    };

    for (auto& viaDefsLayer : viaDefs_) {
      std::sort(viaDefsLayer.begin(), viaDefsLayer.end(), cmp);
    }
  }
  // e.g., getViaDef()     --> get best one-cut viadef
  // e.g., getViaDef(1)    --> get best one-cut viadef
  // e.g., getViaDef(2)    --> get best two-cut viadef
  // e.g., getViaDef(1, 1) --> get 1st alternative one-cut viadef
  const frViaDef* getViaDef(int numCut = 1, int idx = 0) const
  {
    return viaDefs_[numCut - 1][idx];
  }
  frPinAccess* getPinAccess() const { return aps_; }
  frCost getCost() const { return (int) typeL_ + 4 * (int) typeH_; }
  frAccessPointEnum getType(bool isL) const
  {
    if (isL) {
      return typeL_;
    }
    return typeH_;
  }
  bool isViaAllowed() const { return allow_via_; }
  // setters
  void setPoint(const Point& in) { point_ = in; }
  void setLayer(const frLayerNum& layerNum) { layerNum_ = layerNum; }
  void setAccess(const frDirEnum& dir, bool isValid = true)
  {
    switch (dir) {
      case (frDirEnum::E):
        accesses_[0] = isValid;
        break;
      case (frDirEnum::S):
        accesses_[1] = isValid;
        break;
      case (frDirEnum::W):
        accesses_[2] = isValid;
        break;
      case (frDirEnum::N):
        accesses_[3] = isValid;
        break;
      case (frDirEnum::U):
        accesses_[4] = isValid;
        break;
      case (frDirEnum::D):
        accesses_[5] = isValid;
        break;
      default:
        std::cout << "Error: unexpected direction in setValidAccess\n";
    }
  }
  template <std::size_t N>
  void setMultipleAccesses(const frDirEnum (&dirArray)[N], bool isValid = true)
  {
    for (std::size_t i = 0; i < N; ++i) {
      setAccess(dirArray[i], isValid);
    }
  }
  void addViaDef(const frViaDef* in);
  void addToPinAccess(frPinAccess* in) { aps_ = in; }
  void setType(frAccessPointEnum in, bool isL = true)
  {
    if (isL) {
      typeL_ = in;
    } else {
      typeH_ = in;
    }
  }
  void setAllowVia(bool in) { allow_via_ = in; }
  // others
  frBlockObjectEnum typeId() const override { return frcAccessPoint; }
  frCoord x() const { return point_.x(); }
  frCoord y() const { return point_.y(); }

  void addPathSeg(const frPathSeg& ps) { pathSegs_.emplace_back(ps); }
  std::vector<frPathSeg>& getPathSegs() { return pathSegs_; }

 private:
  Point point_;
  frLayerNum layerNum_{0};
  // 0 = E, 1 = S, 2 = W, 3 = N, 4 = U, 5 = D
  std::vector<bool> accesses_ = std::vector<bool>(6, false);
  // cut number -> up-via access map
  std::vector<std::vector<const frViaDef*>> viaDefs_;
  frAccessPointEnum typeL_{frAccessPointEnum::OnGrid};
  frAccessPointEnum typeH_{frAccessPointEnum::OnGrid};
  frPinAccess* aps_{nullptr};
  std::vector<frPathSeg> pathSegs_;
  bool allow_via_{false};
  template <class Archive>
  void serialize(Archive& ar, unsigned int version);
  friend class boost::serialization::access;
};

class frPinAccess : public frBlockObject
{
 public:
  frPinAccess() = default;
  frPinAccess(const frPinAccess& rhs) : frBlockObject(rhs), pin_(rhs.pin_)
  {
    aps_.clear();
    for (const auto& ap : rhs.aps_) {
      aps_.push_back(std::make_unique<frAccessPoint>(*ap));
      aps_.back()->addToPinAccess(this);
    }
  }

  // getters
  const std::vector<std::unique_ptr<frAccessPoint>>& getAccessPoints() const
  {
    return aps_;
  }
  frAccessPoint* getAccessPoint(int idx) const { return aps_[idx].get(); }
  int getNumAccessPoints() const { return aps_.size(); }
  frPin* getPin() const { return pin_; }
  // setters
  void addAccessPoint(std::unique_ptr<frAccessPoint> in)
  {
    in->setId(aps_.size());
    in->addToPinAccess(this);
    aps_.push_back(std::move(in));
  }
  void setPin(frPin* in) { pin_ = in; }
  // others
  frBlockObjectEnum typeId() const override { return frcPinAccess; }

 private:
  std::vector<std::unique_ptr<frAccessPoint>> aps_;
  frPin* pin_{nullptr};
  template <class Archive>
  void serialize(Archive& ar, unsigned int version);
  friend class boost::serialization::access;
};
}  // namespace drt
