#pragma once

#include <QtGui>

class Shape
{
public:
	Shape();
	virtual ~Shape();
	virtual void Draw(QPainter& paint) = 0;
	//virtual void add_vertex(QPoint p); // used only drawing polygons and freehand shapes
	virtual void set_start(QPoint s);
	virtual void set_end(QPoint e);

public:
	enum Type
	{
		kDefault = 0,
		kLine = 1,
		kRect = 2,
		kCirc = 3,
		kEllipse = 4,
		kPoly = 5,
	};

protected:
	// for lines, circles and rectangles, there are only 2 control points
	//std::vector<QPoint*> controlPoints;
	QPoint start;
	QPoint end;
};

