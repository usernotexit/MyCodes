#include "CEllipse.h"

CEllipse::CEllipse()
{
}

CEllipse::~CEllipse()
{
}

void CEllipse::Draw(QPainter& painter)
{
	// 在起点和终点确定的矩形中内嵌的椭圆
	QRect rect(start, end);
	painter.drawEllipse(rect);
}