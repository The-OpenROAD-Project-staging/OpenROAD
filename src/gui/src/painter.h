// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2023-2025, The OpenROAD Authors

#pragma once

#include <QPainter>
#include <cmath>
#include <string>
#include <vector>

#include "gui/gui.h"
#include "odb/geom.h"
#include "options.h"
#include "ruler.h"

namespace gui {

// This class wraps the QPainter in the abstract Painter API for
// Renderer instances to use.
class GuiPainter : public Painter
{
 public:
  GuiPainter(QPainter* painter,
             Options* options,
             const odb::Rect& bounds,
             qreal pixels_per_dbu,
             int dbu_per_micron)
      : Painter(options, bounds, pixels_per_dbu),
        painter_(painter),
        dbu_per_micron_(dbu_per_micron)
  {
  }

  Color getPenColor() override
  {
    QColor color = painter_->pen().color();
    return Color(color.red(), color.green(), color.blue(), color.alpha());
  }

  void setPen(odb::dbTechLayer* layer, bool cosmetic = false) override
  {
    QPen pen(getOptions()->color(layer));
    pen.setCosmetic(cosmetic);
    painter_->setPen(pen);
  }

  void setPen(const Color& color, bool cosmetic = false, int width = 1) override
  {
    QPen pen(QColor(color.r, color.g, color.b, color.a));
    pen.setCosmetic(cosmetic);
    pen.setWidth(width);
    painter_->setPen(pen);
  }

  void setPenWidth(int width) override
  {
    QPen pen(painter_->pen().color());
    pen.setCosmetic(painter_->pen().isCosmetic());
    pen.setWidth(width);
    painter_->setPen(pen);
  }

  void setBrush(odb::dbTechLayer* layer, int alpha = -1) override
  {
    QColor color = getOptions()->color(layer);
    Qt::BrushStyle brush_pattern = getOptions()->pattern(layer);
    if (alpha >= 0) {
      color.setAlpha(alpha);
    }
    painter_->setBrush(QBrush(color, brush_pattern));
  }

  void setBrush(const Color& color, const Brush& style = Brush::kSolid) override
  {
    const QColor qcolor(color.r, color.g, color.b, color.a);

    Qt::BrushStyle brush_pattern;
    if (color == Painter::kTransparent) {
      // if color is transparent, make it no brush
      brush_pattern = Qt::NoBrush;
    } else {
      switch (style) {
        case kNone:
          brush_pattern = Qt::NoBrush;
          break;
        case kDiagonal:
          brush_pattern = Qt::DiagCrossPattern;
          break;
        case kCross:
          brush_pattern = Qt::CrossPattern;
          break;
        case kDots:
          brush_pattern = Qt::Dense6Pattern;
          break;
        case kSolid:
        default:
          brush_pattern = Qt::SolidPattern;
          break;
      }
    }

    painter_->setBrush(QBrush(qcolor, brush_pattern));
  }

  void setFont(const Font& font) override;

  void saveState() override { painter_->save(); }

  void restoreState() override { painter_->restore(); }

  void drawOctagon(const odb::Oct& oct) override
  {
    std::vector<odb::Point> points = oct.getPoints();
    drawPolygon(points);
  }
  void drawRect(const odb::Rect& rect, int round_x, int round_y) override
  {
    if (round_x > 0 || round_y > 0) {
      painter_->drawRoundedRect(
          QRect(rect.xMin(), rect.yMin(), rect.dx(), rect.dy()),
          round_x,
          round_y,
          Qt::RelativeSize);
    } else {
      painter_->drawRect(QRect(rect.xMin(), rect.yMin(), rect.dx(), rect.dy()));
    }
  }
  void drawPolygon(const odb::Polygon& polygon) override
  {
    drawPolygon(polygon.getPoints());
  }
  void drawPolygon(const std::vector<odb::Point>& points) override
  {
    QPolygon poly;
    for (const auto& pt : points) {
      poly.append(QPoint(pt.x(), pt.y()));
    }
    painter_->drawPolygon(poly);
  }
  void drawLine(const odb::Point& p1, const odb::Point& p2) override
  {
    painter_->drawLine(p1.x(), p1.y(), p2.x(), p2.y());
  }
  using Painter::drawLine;

  void drawCircle(int x, int y, int r) override
  {
    painter_->drawEllipse(QPoint(x, y), r, r);
  }

  void drawX(int x, int y, int size) override
  {
    const int o = size / 2;
    painter_->drawLine(x - o, y - o, x + o, y + o);
    painter_->drawLine(x - o, y + o, x + o, y - o);
  }

  odb::Point determineStringOrigin(int x,
                                   int y,
                                   Anchor anchor,
                                   const QString& text,
                                   bool rotate_90 = false);

  // NOTE: The constant height text s drawn with this function, hence
  //       the transformation is mapped to the base transformation and
  //       the world co-ordinates are mapped to the window co-ordinates
  //       before drawing.
  void drawString(int x,
                  int y,
                  Anchor anchor,
                  const std::string& s,
                  bool rotate_90 = false) override;

  odb::Rect stringBoundaries(int x,
                             int y,
                             Anchor anchor,
                             const std::string& s) override
  {
    const QString text = QString::fromStdString(s);
    const odb::Point origin = determineStringOrigin(x, y, anchor, text);
    const qreal scale_adjust = 1.0 / getPixelsPerDBU();

    const QRect text_bbox = painter_->fontMetrics().boundingRect(text);
    const int x_min = origin.x() - text_bbox.left() * scale_adjust;
    const int y_min = origin.y() - text_bbox.bottom() * scale_adjust;
    const int x_max = x_min + text_bbox.width() * scale_adjust;
    const int y_max = y_min + text_bbox.height() * scale_adjust;
    return {x_min, y_min, x_max, y_max};
  }

  void drawRuler(int x0,
                 int y0,
                 int x1,
                 int y1,
                 bool euclidian = true,
                 const std::string& label = "") override
  {
    if (euclidian) {
      drawRuler(x0, y0, x1, y1, label);
    } else {
      const int x_dist = std::abs(x0 - x1);
      const int y_dist = std::abs(y0 - y1);
      std::string x_label = label;
      std::string y_label;
      if (y_dist > x_dist) {
        std::swap(x_label, y_label);
      }
      const odb::Point mid_pt = Ruler::getManhattanJoinPt({x0, y0}, {x1, y1});
      drawRuler(x0, y0, mid_pt.x(), mid_pt.y(), x_label);
      drawRuler(mid_pt.x(), mid_pt.y(), x1, y1, y_label);
    }
  }

  QPainter* getPainter() { return painter_; }

 private:
  QPainter* painter_;
  int dbu_per_micron_;

  void drawRuler(int x0, int y0, int x1, int y1, const std::string& label);
};

}  // namespace gui
