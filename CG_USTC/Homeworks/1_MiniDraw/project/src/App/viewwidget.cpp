#include "viewwidget.h"
#include "Line.h"
#include "Rect.h"
#include "Circle.h"
#include "Polygon.h"
#include "CEllipse.h"

ViewWidget::ViewWidget(QWidget* parent)
	: QWidget(parent)
{
	ui.setupUi(this);
	draw_status_ = false;
	shape_ = NULL;
	type_ = Shape::kDefault;
}

ViewWidget::~ViewWidget()
{
}

void ViewWidget::setLine()
{
	type_ = Shape::kLine;
}

void ViewWidget::setRect()
{
	type_ = Shape::kRect;
}

void ViewWidget::setCirc()
{
	type_ = Shape::kCirc;
}

void ViewWidget::setPoly()
{
	type_ = Shape::kPoly;
}

void ViewWidget::setEllipse()
{
	type_ = Shape::kEllipse;
}

void ViewWidget::undo()
{
	// now this function is equivalent to 'Delete'
	type_ = Shape::kDefault;
	shape_ = NULL;
	if(!shape_list_.empty())
		shape_list_.pop_back();
}

void ViewWidget::clearAll()
{
	type_ = Shape::kDefault;
	shape_ = NULL;
	shape_list_.clear();
}

void ViewWidget::mousePressEvent(QMouseEvent* event)
{
	if (Qt::LeftButton == event->button())
	{
		if (!draw_status_)
			switch (type_)
			{
			case Shape::kLine:
				shape_ = new Line();
				break;
			case Shape::kDefault:
				break;
			case Shape::kRect:
				shape_ = new Rect();
				break;
			case Shape::kCirc:
				shape_ = new Circle();
				break;
			case Shape::kPoly:
				shape_ =  new CPolygon();
				shape_->set_start(event->pos());
				break;
			case Shape::kEllipse:
				shape_ = new CEllipse();
				shape_->set_start(event->pos());
			}
		if (shape_ != NULL)
		{
			draw_status_ = true;
			start_point_ = end_point_ = event->pos();
			shape_->set_start(start_point_);
			shape_->set_end(end_point_);
		}
	}
	if (Qt::RightButton == event->button())
	{
		// 在绘制多边形时用于终止作图
		if (draw_status_ && type_ == Shape::kPoly)
		{
			draw_status_ = false;
			shape_list_.push_back(shape_);
			// type_ = Shape::kDefault;
		}
	}

	update();
}

void ViewWidget::mouseMoveEvent(QMouseEvent* event)
{
	
	if (draw_status_ && shape_ != NULL)
	{
		end_point_ = event->pos();
		shape_->set_end(end_point_);
	}
}

void ViewWidget::mouseReleaseEvent(QMouseEvent* event)
{
	if (shape_ != NULL && type_ != Shape::kPoly)
	{
		draw_status_ = false;
		shape_list_.push_back(shape_);
		shape_ = NULL;
		// type_ = Shape::kDefault;
	}
}

void ViewWidget::mouseDoubleClickEvent(QMouseEvent* event)
{
	/// 在绘制多边形时用于终止作图 // 已废除，请用单击右键的方式
	/*
	if (draw_status_ && type_ == Shape::kPoly)
	{
		draw_status_ = false;
		shape_list_.push_back(shape_);
		shape_ = NULL;
	}
	*/
}

void ViewWidget::paintEvent(QPaintEvent*)
{
	QPainter painter(this);

	for (auto &shape_tmp : shape_list_ )
	{
		shape_tmp->Draw(painter);
	}

	if (shape_ != NULL) {
		shape_->Draw(painter);
	}

	update();
}