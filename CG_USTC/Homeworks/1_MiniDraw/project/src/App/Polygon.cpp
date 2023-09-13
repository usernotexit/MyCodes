#include "Polygon.h"

CPolygon::CPolygon()
{
	vertices.clear();
}

CPolygon::CPolygon(QPoint start)
{
	vertices.clear();
	vertices.push_back(start);
	n = vertices.size();
}

CPolygon::~CPolygon()
{
}

void CPolygon::Draw(QPainter& painter)
{
	const auto& points = vertices.data();
	painter.drawPolygon(points, n);
}

void CPolygon::set_start(QPoint p)
{
	//printf("CPolygon set_start\n");
	vertices.push_back(p); n = vertices.size();
}

void CPolygon::set_end(QPoint p)
{
	vertices.pop_back();
	//printf("CPolygon set_end\n");
	vertices.push_back(p); n = vertices.size();
}