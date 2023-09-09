#pragma once
#include "shape.h"

class Circle :
    public Shape
{
public:
    Circle();
    ~Circle();

    void Draw(QPainter& painter);
    void Draw2(QPainter& painter);
};

