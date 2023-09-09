#include "Circle.h"

auto get_min = [](int a, int b) { return a > b ? b : a;};
auto sqrt_ave = [](int a, int b) { return sqrt(a * a + b * b);};

Circle::Circle()
{
}

Circle::~Circle()
{
}

void Circle::Draw2(QPainter &painter)
{
	// 确定半径（直径）
	int r2 = get_min(abs(end.x() - start.x()), abs(end.y() - start.y()));
	// 确定位置（左上）
	int pos_x = start.x() < end.x() ? start.x() : start.x() - r2;
	int pos_y = start.y() < end.y() ? start.y() : start.y() - r2;
	painter.drawArc(pos_x, pos_y, r2, r2, 0, 360 * 16);
}

void Circle::Draw(QPainter& painter)
{
	// 更适合数学人体质的画圆方法，起点为圆心，终点在圆上
	painter.setPen(QPen(Qt::blue, 4));

	int r = round(sqrt_ave(abs(end.x() - start.x()), abs(end.y() - start.y())));
	painter.drawArc(start.x() - r, start.y() - r, 2 * r, 2 * r, 0, 360 * 16);

	painter.setPen(QPen());
}