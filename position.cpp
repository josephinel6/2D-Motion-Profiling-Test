#include "position.hpp"
#include <cmath>

Position::Position()
{
    this->x = 0;
    this->y = 0;
    this->heading = 0;
}

Position::Position(double x, double y, double heading)
{
    this->x = x;
    this->y = y;
    this->heading = heading;
}

void Position::setPosition(double x, double y, double heading)
{
    this->x = x;
    this->y = y;
    this->heading = heading;
}

double Position::angleTo(Position secondPosition)
{
    return atan2((secondPosition.y - this->y), (secondPosition.x - this->x));
}

double Position::distanceFrom(Position secondPosition)
{
    return sqrt(pow(secondPosition.x - this->x, 2) + pow(secondPosition.y - this->y, 2));
}
