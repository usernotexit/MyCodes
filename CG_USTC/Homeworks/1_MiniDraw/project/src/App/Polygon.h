#pragma once
#include "shape.h"

//多边形的绘制与矩形、圆等不同
class CPolygon :
    public Shape
{
public:
    CPolygon();
    CPolygon(QPoint start); // 保证初始点数量
    ~CPolygon();

    void Draw(QPainter &painter);
    void set_start(QPoint p);// 重载shape的set函数
    void set_end(QPoint p);

protected:
    int n;
    std::vector<QPoint> vertices;
};