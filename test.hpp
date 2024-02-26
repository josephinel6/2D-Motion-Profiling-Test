#include "position.hpp"

class CubicBezier
{
public:
    CubicBezier(Position p0, Position p1, Position p2, Position p3);
    Position getPoint(double t);
    double getVelocity(double t);
    double getAcceleration(double t);

private:
    Position p0;
    Position p1;
    Position p2;
    Position p3;
};
