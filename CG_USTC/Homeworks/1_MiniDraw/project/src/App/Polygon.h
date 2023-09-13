#pragma once
#include "shape.h"

//����εĻ�������Ρ�Բ�Ȳ�ͬ
class CPolygon :
    public Shape
{
public:
    CPolygon();
    CPolygon(QPoint start); // ��֤��ʼ������
    ~CPolygon();

    void Draw(QPainter &painter);
    void set_start(QPoint p);// ����shape��set����
    void set_end(QPoint p);

protected:
    int n;
    std::vector<QPoint> vertices;
};