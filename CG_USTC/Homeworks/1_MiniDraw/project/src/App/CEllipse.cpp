#include "CEllipse.h"

CEllipse::CEllipse()
{
}

CEllipse::~CEllipse()
{
}

void CEllipse::Draw(QPainter& painter)
{
	// �������յ�ȷ���ľ�������Ƕ����Բ
	QRect rect(start, end);
	painter.drawEllipse(rect);
}