#pragma once
#include "shape.h"
class CEllipse :
    public Shape
{
public:
    CEllipse();
    ~CEllipse();

    void Draw(QPainter& painter);
};

