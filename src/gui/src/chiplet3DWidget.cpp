// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025, The OpenROAD Authors

#include "chiplet3DWidget.h"

#include <QMouseEvent>
#include <QPaintEvent>
#include <QPainter>
#include <algorithm>
#include <any>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "gui/gui.h"
#include "odb/db.h"
#include "odb/dbObject.h"
#include "odb/dbTransform.h"
#include "odb/geom.h"
#include "odb/unfoldedModel.h"
#include "utl/Logger.h"

namespace {
constexpr float kInitialDistanceFactor = 3.0f;
constexpr float kLayerGapFactor = 2.0f;
constexpr float kMinDistanceCheck = 100.0f;
constexpr float kDefaultDistance = 1000.0f;
constexpr float kDefaultSafeSize = 1000.0f;
constexpr float kNearFarFactor = 2.0f;
constexpr float kMinZNear = 10.0f;
constexpr float kMaxZFarOffset = 10000.0f;
constexpr float kGridSizeFactor = 1.5f;
constexpr float kGridSteps = 5.0f;
constexpr float kGridZOffsetFactor = 0.05f;
constexpr int kGridLineCount = 5;
constexpr float kRotationSensitivity = 2.0f;
constexpr float kPanSensitivity = 0.002f;
constexpr float kZoomInFactor = 0.9f;
constexpr float kZoomOutFactor = 1.1f;

static const std::array<QVector3D, 7> kColorPalette = {{
    {0.0f, 1.0f, 0.0f},  // Green
    {0.0f, 0.5f, 0.5f},  // Teal (Replaced Yellow)
    {0.0f, 1.0f, 1.0f},  // Cyan
    {1.0f, 0.0f, 1.0f},  // Magenta
    {1.0f, 0.5f, 0.0f},  // Orange
    {0.5f, 0.5f, 1.0f},  // Blue-ish
    {1.0f, 0.0f, 0.0f}   // Red
}};
}  // namespace

namespace gui {

Chiplet3DWidget::Chiplet3DWidget(QWidget* parent) : QWidget(parent)
{
  // Optimize for painting speed
  setAttribute(Qt::WA_OpaquePaintEvent);
}

void Chiplet3DWidget::setChip(odb::dbChip* chip)
{
  chip_ = chip;
  buildGeometries();
  update();
}

void Chiplet3DWidget::setLogger(utl::Logger* logger)
{
  logger_ = logger;
}

void Chiplet3DWidget::setSelection(const SelectionSet& selection)
{
  selection_ = selection;
  focus_ = Selected();
  update();
}

void Chiplet3DWidget::selectionFocus(const Selected& focus)
{
  focus_ = focus;
  update();
}

void Chiplet3DWidget::buildGeometries()
{
  if (!chip_) {
    return;
  }
  const odb::Cuboid global_cuboid = chip_->getCuboid();
  const odb::UnfoldedModel model(logger_, chip_);
  const odb::dbTransform center_transform
      = odb::dbTransform(odb::Point3D(-global_cuboid.xCenter(),
                                      -global_cuboid.yCenter(),
                                      -global_cuboid.zCenter()));

  vertices_.clear();
  indices_faces_.clear();
  chip_objects_.clear();
  face_to_chip_index_.clear();

  // Center and Camera calculations
  const float cx = (global_cuboid.xMin() + global_cuboid.xMax()) / 2.0f;
  const float cy = (global_cuboid.yMin() + global_cuboid.yMax()) / 2.0f;
  center_ = QVector3D(cx, cy, 0.0f);

  const float dx = global_cuboid.dx();
  const float dy = global_cuboid.dy();
  const float dz = global_cuboid.dz() * kLayerGapFactor;
  bounding_radius_ = std::sqrt(dx * dx + dy * dy + dz * dz) / 2.0f;

  distance_ = bounding_radius_ * kInitialDistanceFactor;
  if (distance_ < kMinDistanceCheck) {
    distance_ = kDefaultDistance;
  }

  int index = 0;
  for (const auto& chip : model.getChips()) {
    const odb::dbObject* obj = nullptr;
    if (!chip.chip_inst_path.empty()) {
      obj = chip.chip_inst_path.back();
    } else {
      obj = chip_;
    }
    chip_objects_.push_back(obj);

    odb::Cuboid draw_cuboid = chip.cuboid;
    center_transform.apply(draw_cuboid);
    // Color by Depth (proportional to Z)
    const QVector3D color = kColorPalette[index++ % kColorPalette.size()];

    const uint32_t base = vertices_.size();
    for (const auto& p : draw_cuboid.getPoints()) {
      vertices_.push_back({QVector3D(p.x(), p.y(), p.z()), color});
    }

    // Add face indices (6 faces, 4 vertices each)
    // Winding order: CCW for outward normal
    const uint32_t faces[24] = {0, 3, 2, 1,   // Bottom
                                4, 5, 6, 7,   // Top
                                0, 1, 5, 4,   // Side 1
                                1, 2, 6, 5,   // Side 2
                                2, 3, 7, 6,   // Side 3
                                3, 0, 4, 7};  // Side 4

    const int current_chip_idx = chip_objects_.size() - 1;
    // 6 faces
    for (int f = 0; f < 6; ++f) {
      face_to_chip_index_.push_back(current_chip_idx);
    }

    for (const uint32_t i : faces) {
      indices_faces_.push_back(base + i);
    }
  }
}

void Chiplet3DWidget::paintEvent(QPaintEvent* event)
{
  QPainter painter(this);
  painter.setRenderHint(QPainter::Antialiasing);

  // 1. Clear Background
  painter.fillRect(rect(), QColor(26, 26, 26));  // Approx 0.1f grey

  // 2. Setup Camera Matrices

  // Use bounding_radius_ which is rotation-invariant
  float safe_size = bounding_radius_;
  if (safe_size < 1.0f) {
    safe_size = kDefaultSafeSize;
  }

  const float zNear
      = std::max(distance_ - safe_size * kNearFarFactor, kMinZNear);
  float zFar = distance_ + safe_size * kNearFarFactor;

  if (zFar < zNear + kMinDistanceCheck) {
    zFar = zNear + kMaxZFarOffset;
  }

  // Projection Matrix
  const qreal aspect = qreal(width()) / qreal(height() ? height() : 1);
  QMatrix4x4 projection;
  projection.perspective(45.0f, aspect, zNear, zFar);

  // ModelView Matrix
  QMatrix4x4 modelView;
  modelView.translate(pan_x_, pan_y_, -distance_);
  modelView.rotate(rotation_);

  // 3. Draw Grid
  painter.setPen(QPen(QColor(76, 76, 76), 1));  // 0.3f grey
  const float grid_size
      = safe_size * kGridSizeFactor;  // Adjust grid to scene size
  const float step = grid_size / kGridSteps;
  const float grid_z = -center_.z() - (distance_ * kGridZOffsetFactor * 0.5f);

  for (int i = -kGridLineCount; i <= kGridLineCount; ++i) {
    const float pos = i * step;
    // Lines parallel to Z (or Y in local, depending on orientation)
    // The original code drew grid on Z plane.
    drawLine3D(painter,
               QVector3D(pos, -grid_size, grid_z),
               QVector3D(pos, grid_size, grid_z),
               QColor(76, 76, 76),
               modelView,
               projection,
               rect());
    drawLine3D(painter,
               QVector3D(-grid_size, pos, grid_z),
               QVector3D(grid_size, pos, grid_z),
               QColor(76, 76, 76),
               modelView,
               projection,
               rect());
  }

  // 4. Draw Chiplets
  if (vertices_.empty()) {
    return;
  }

  // Sort Faces by Depth (Painter's Algorithm)
  // We use the farthest vertex (minimum Z in view space) as the sorting key.
  // This helps when a large polygon (like a base chiplet) spans a large Z
  // range. By using the farthest point, we ensure it's drawn earlier
  // (background) rather than later based on a centroid that might be closer
  // than other objects.
  struct FaceDepth
  {
    size_t index_offset;
    float depth;
  };
  std::vector<FaceDepth> sorted_faces;
  sorted_faces.reserve(indices_faces_.size() / 4);

  for (size_t i = 0; i < indices_faces_.size(); i += 4) {
    float min_z = std::numeric_limits<float>::max();
    for (int j = 0; j < 4; ++j) {
      const uint32_t idx = indices_faces_[i + j];
      if (idx < vertices_.size()) {
        const QVector3D view_pos = modelView * vertices_[idx].position;
        min_z = std::min(view_pos.z(), min_z);
      }
    }
    sorted_faces.push_back({i, min_z});
  }

  // Sort back-to-front (smaller Z in view space is further away if looking down
  // -Z, but standard OpenGL view is looking down -Z. View space Z usually is
  // negative in front of camera. So "further away" means more negative. Wait,
  // let's check coordinate system.
  // Standard: Camera at origin looking at -Z.
  // Objects at Z=-10 are further than Z=-5.
  // We want to draw -10 first. So sort ascending Z (more negative first).
  std::ranges::stable_sort(
      sorted_faces,
      [](const FaceDepth& a, const FaceDepth& b) { return a.depth < b.depth; });

  // Draw Faces
  for (const auto& face : sorted_faces) {
    std::vector<QVector3D> face_verts;
    face_verts.reserve(4);
    QColor color;

    // Determine if highlighted
    const size_t face_idx = face.index_offset / 4;
    bool is_highlighted = false;
    if (face_idx < face_to_chip_index_.size()) {
      const int chip_idx = face_to_chip_index_[face_idx];
      if (chip_idx >= 0 && chip_idx < chip_objects_.size()) {
        const odb::dbObject* obj = chip_objects_[chip_idx];
        if (obj) {
          odb::dbBlock* chiplet_block = nullptr;
          if (obj->getObjectType() == odb::dbChipObj) {
            chiplet_block = ((odb::dbChip*) obj)->getBlock();
          } else if (obj->getObjectType() == odb::dbChipInstObj) {
            chiplet_block
                = ((odb::dbChipInst*) obj)->getMasterChip()->getBlock();
          }

          // Check selection set
          for (const auto& sel : selection_) {
            const std::any& any_obj = sel.getObject();
            if (any_obj.type() == typeid(odb::dbChipInst*)) {
              if (std::any_cast<odb::dbChipInst*>(any_obj) == obj) {
                is_highlighted = true;
                break;
              }
            } else if (any_obj.type() == typeid(odb::dbChip*)) {
              if (std::any_cast<odb::dbChip*>(any_obj) == obj) {
                is_highlighted = true;
                break;
              }
            } else if (any_obj.type() == typeid(odb::dbMarker*)) {
              auto* marker = std::any_cast<odb::dbMarker*>(any_obj);
              // Check ownership
              if (chiplet_block) {
                auto* category = marker->getCategory();
                if (category) {
                  auto* parent = category->getParent();
                  if (parent == chiplet_block) {
                    is_highlighted = true;
                    break;
                  }
                }
              }
              // Check sources
              for (odb::dbObject* src : marker->getSources()) {
                if (src == obj) {
                  is_highlighted = true;
                  break;
                }
                // Check name match for ChipInst vs Inst
                if (obj->getObjectType() == odb::dbChipInstObj
                    && src->getObjectType() == odb::dbInstObj) {
                  const auto* chip_inst = (odb::dbChipInst*) obj;
                  const auto* inst = (odb::dbInst*) src;
                  if (chip_inst && inst
                      && chip_inst->getName() == inst->getName()) {
                    is_highlighted = true;
                    break;
                  }
                }
              }
              if (is_highlighted) {
                break;
              }
            }
          }

          // Check focus object
          if (!is_highlighted && focus_) {
            const std::any& any_obj = focus_.getObject();
            if (any_obj.type() == typeid(odb::dbChipInst*)) {
              if (std::any_cast<odb::dbChipInst*>(any_obj) == obj) {
                is_highlighted = true;
              }
            } else if (any_obj.type() == typeid(odb::dbChip*)) {
              if (std::any_cast<odb::dbChip*>(any_obj) == obj) {
                is_highlighted = true;
              }
            } else if (any_obj.type() == typeid(odb::dbMarker*)) {
              auto* marker = std::any_cast<odb::dbMarker*>(any_obj);
              // Check ownership
              if (chiplet_block) {
                auto* category = marker->getCategory();
                if (category) {
                  auto* parent = category->getParent();
                  if (parent == chiplet_block) {
                    is_highlighted = true;
                  }
                }
              }
              // Check sources
              for (odb::dbObject* src : marker->getSources()) {
                if (src == obj) {
                  is_highlighted = true;
                  break;
                }
                // Check name match for ChipInst vs Inst
                if (obj->getObjectType() == odb::dbChipInstObj
                    && src->getObjectType() == odb::dbInstObj) {
                  const auto* chip_inst = (odb::dbChipInst*) obj;
                  const auto* inst = (odb::dbInst*) src;
                  if (chip_inst && inst
                      && chip_inst->getName() == inst->getName()) {
                    is_highlighted = true;
                    break;
                  }
                }
              }
            }
          }
        }
      }
    }

    if (is_highlighted) {
      // Yellow highlight
      color.setRgbF(1.0f, 1.0f, 0.0f);
    } else {
      // Normal color retrieval
      // We need to retrieve it from vertices since we don't store it elsewhere
      // efficiently here
      const uint32_t first_v_idx = indices_faces_[face.index_offset];
      if (first_v_idx < vertices_.size()) {
        const VertexData& v = vertices_[first_v_idx];
        color.setRgbF(v.color.x(), v.color.y(), v.color.z());
      }
    }

    for (int j = 0; j < 4; ++j) {
      const uint32_t idx = indices_faces_[face.index_offset + j];
      if (idx < vertices_.size()) {
        const VertexData& v = vertices_[idx];
        face_verts.push_back(v.position);
      }
    }

    // Use solid color to prevent confusing transparency overlaps
    color.setAlpha(255);

    drawFace3D(painter, face_verts, color, modelView, projection, rect());
  }
}

// Helper: Projects 3D points to 2D, Handles Z-Clipping
void Chiplet3DWidget::drawFace3D(QPainter& painter,
                                 const std::vector<QVector3D>& face_verts_world,
                                 const QColor& color,
                                 const QMatrix4x4& modelView,
                                 const QMatrix4x4& projection,
                                 const QRect& viewport)
{
  if (face_verts_world.empty()) {
    return;
  }

  QPolygonF polygon;
  bool all_behind = true;
  const float kClipZ = -0.1f;

  // Simple approach: project all points. If any is behind camera, we have to
  // clip the polygon. Clipping a polygon against a plane is non-trivial
  // (Sutherland-Hodgman). For simplicity in this widget, if any point is
  // behind, we might skip or clamp. Given this is a simple visualizer, let's
  // try projecting and if Z > kClipZ, we don't draw or we clamp.
  //
  // Better approach for a simple visualizer: Just project. The QMatrix4x4 map
  // handles w-division. If w <= 0 or close, result is garbage.

  // We will do a simplified clipping: if the centroid is behind, skip.
  // Or check if ALL points are behind.
  // Correct way: Sutherland-Hodgman.
  //
  // Let's implement a very basic check: if any point is behind near plane, skip
  // face. This causes popping but is safe.
  // To improve: implement Sutherland-Hodgman clipping against Z = kClipZ.

  std::vector<QVector3D> current_poly = face_verts_world;

  // Transform to View Space
  for (auto& p : current_poly) {
    p = modelView * p;
    if (p.z() < kClipZ) {
      all_behind = false;
    }
  }

  if (all_behind) {
    return;  // All vertices behind camera (Z > kClipZ? No, Z is negative in
             // front).
    // Wait, Camera looks down -Z. Objects in front have Negative Z.
    // Behind camera have Positive Z (or > -near).
    // kClipZ = -0.1.
    // If p.z() > kClipZ, it is behind the near plane (or very close).
    // So if ALL p.z() > kClipZ, we skip.
  }

  // Check if all are behind (Z > -0.1)
  bool all_clipped = true;
  for (const auto& p : current_poly) {
    if (p.z() <= kClipZ) {
      all_clipped = false;
      break;
    }
  }
  if (all_clipped) {
    return;
  }

  // Back-Face Culling & Flat Shading
  // Compute normal of the face in View Space.
  // We need at least 3 points. Since it's a quad or triangle, we use the
  // first 3. Note: Points are already in View Space in 'current_poly'.
  QColor shaded_color = color;
  if (current_poly.size() >= 3) {
    const QVector3D v1 = current_poly[1] - current_poly[0];
    const QVector3D v2 = current_poly[2] - current_poly[0];
    const QVector3D normal = QVector3D::crossProduct(v1, v2);

    // If normal.z > 0, the face is pointing towards the camera (assuming
    // standard OpenGL view looking down -Z). Wait, in OpenGL View Space, Camera
    // is at (0,0,0) looking at -Z. A face visible to the camera should have a
    // normal pointing roughly towards +Z (back at the camera). So if Normal.z >
    // 0, it is visible. If Normal.z <= 0, it is back-facing.
    if (normal.z() <= 0) {
      return;
    }

    // Flat Shading
    // 1. Fixed Light Direction (in View Space)
    //    We assume light comes from top-right-front relative to camera.
    //    Normalized vector: (1, 1, 1) -> (0.577, 0.577, 0.577)
    const QVector3D light_dir = QVector3D(0.5f, 0.5f, 1.0f).normalized();

    // 2. Normal Vector (already computed 'normal' is not normalized)
    const QVector3D n = normal.normalized();

    // 3. Dot Product (Lambertian)
    //    Ideally, we want dot(n, light_dir).
    //    Since both n and light_dir are in View Space (and light is defined
    //    relative to camera), this works.
    float brightness = QVector3D::dotProduct(n, light_dir);

    // Clamp brightness to [0, 1] range (handle back-facing light slightly
    // gracefully or just clamp) Ambient light: 0.3 Diffuse factor: 0.7
    brightness = std::max(0.0f, brightness);
    const float intensity = 0.3f + 0.7f * brightness;

    // 4. Modulate Color
    shaded_color.setRedF(std::min(1.0, color.redF() * intensity));
    shaded_color.setGreenF(std::min(1.0, color.greenF() * intensity));
    shaded_color.setBlueF(std::min(1.0, color.blueF() * intensity));
    shaded_color.setAlpha(color.alpha());
  }

  // Sutherland-Hodgman Clip against Z = kClipZ plane
  std::vector<QVector3D> clipped_poly;
  for (size_t i = 0; i < current_poly.size(); ++i) {
    const QVector3D& p1
        = current_poly[(i + current_poly.size() - 1) % current_poly.size()];
    const QVector3D& p2 = current_poly[i];

    bool inside1 = (p1.z() <= kClipZ);
    bool inside2 = (p2.z() <= kClipZ);

    if (inside1 != inside2) {
      float t = (kClipZ - p1.z()) / (p2.z() - p1.z());
      clipped_poly.push_back(p1 + (p2 - p1) * t);
    }
    if (inside2) {
      clipped_poly.push_back(p2);
    }
  }

  if (clipped_poly.size() < 3) {
    return;
  }

  // Project to Screen
  for (const auto& p : clipped_poly) {
    QVector3D ndc = projection.map(p);

    const float w = viewport.width();
    const float h = viewport.height();

    const float x = (ndc.x() + 1.0f) * 0.5f * w;
    const float y = (1.0f - ndc.y()) * 0.5f * h;

    polygon << QPointF(x, y);
  }

  painter.setBrush(shaded_color);
  painter.setPen(QPen(shaded_color, 1));
  painter.drawPolygon(polygon);
}

// Helper: Projects 3D points to 2D, Handles Z-Clipping
void Chiplet3DWidget::drawLine3D(QPainter& painter,
                                 const QVector3D& p1_world,
                                 const QVector3D& p2_world,
                                 const QColor& color,
                                 const QMatrix4x4& modelView,
                                 const QMatrix4x4& projection,
                                 const QRect& viewport)
{
  // 1. Transform to View Space
  QVector3D p1_view = modelView * p1_world;
  QVector3D p2_view = modelView * p2_world;

  // 2. Clip against Near Plane
  // In OpenGL view space, camera looks down -Z.
  // Near plane is at z = -nearVal (e.g. -10.0). Points with z > -nearVal are
  // behind the near plane. However, QMatrix4x4::perspective sets up a standard
  // frustum. We can extract the near plane distance from our setup, but a small
  // fixed epsilon usually works for preventing divide-by-zero or "behind head"
  // artifacts.

  // Simple clipping: Clip against Z = -0.1 (very close to camera)
  const float kClipZ = -0.1f;

  const bool p1_visible = p1_view.z() < kClipZ;
  const bool p2_visible = p2_view.z() < kClipZ;

  if (!p1_visible && !p2_visible) {
    return;  // Both behind camera
  }

  if (p1_visible != p2_visible) {
    // Line spans the clip plane. Calculate intersection.
    // t = (kClipZ - z1) / (z2 - z1)
    const float t = (kClipZ - p1_view.z()) / (p2_view.z() - p1_view.z());
    const QVector3D intersection = p1_view + (p2_view - p1_view) * t;

    if (!p1_visible) {
      p1_view = intersection;
    } else {
      p2_view = intersection;
    }
  }

  // NDC or Normalized Device Coordinates
  const QVector3D p1_ndc = projection.map(p1_view);
  const QVector3D p2_ndc = projection.map(p2_view);

  // Map NDC (-1 to 1) to Viewport (0 to width/height)
  const float w = viewport.width();
  const float h = viewport.height();

  // Note: NDC Y is up, Screen Y is down.
  const QPointF s1((p1_ndc.x() + 1.0f) * 0.5f * w,
                   (1.0f - p1_ndc.y()) * 0.5f * h);
  const QPointF s2((p2_ndc.x() + 1.0f) * 0.5f * w,
                   (1.0f - p2_ndc.y()) * 0.5f * h);

  // 4. Draw
  painter.setPen(QPen(color, 2));
  painter.drawLine(s1, s2);
}

void Chiplet3DWidget::mousePressEvent(QMouseEvent* e)
{
  mouse_press_position_ = QVector2D(e->localPos());
}

void Chiplet3DWidget::mouseReleaseEvent(QMouseEvent* e)
{
}

void Chiplet3DWidget::mouseMoveEvent(QMouseEvent* e)
{
  const QVector2D diff = QVector2D(e->localPos()) - mouse_press_position_;
  mouse_press_position_ = QVector2D(e->localPos());

  if (e->buttons() & Qt::LeftButton) {
    // Rotation
    const QVector3D n = QVector3D(diff.y(), diff.x(), 0.0).normalized();
    const float angle = diff.length() / kRotationSensitivity;
    rotation_ = QQuaternion::fromAxisAndAngle(n, angle) * rotation_;
  } else if (e->buttons() & Qt::RightButton) {
    // Pan
    const float scale = distance_ * kPanSensitivity;
    pan_x_ += diff.x() * scale;
    pan_y_ -= diff.y() * scale;
  }
  update();
}

void Chiplet3DWidget::wheelEvent(QWheelEvent* e)
{
  if (e->angleDelta().y() > 0) {
    distance_ *= kZoomInFactor;
  } else {
    distance_ *= kZoomOutFactor;
  }
  update();
}

}  // namespace gui
