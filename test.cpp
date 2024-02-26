#define NUMBER 0.5

#include <iostream>
#include <math.h>
#include <vector>
// #include "test.hpp"
// #include "position.hpp"

// using namespace std;

class Position
{
public:
    double x;
    double y;
    double heading;
    Position();
    Position(double x, double y, double heading);
    void setPosition(double x, double y, double heading);
    Position operator+(Position secondPosition);
    Position operator-(Position secondPosition);
    double angleTo(Position secondPosition);
    double distanceFrom(Position secondPosition);
};

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

class CubicBezier
{
public:
    CubicBezier();
    CubicBezier(Position p0, Position p1, Position p2, Position p3, double numDivisions = 100);
    Position getPoint(double t);
    double getVelocity(double t);
    double getAcceleration(double t);
    double getCurvature(double t);
    double getLength();
    double getMappedT(double t);

private:
    Position p0;
    Position p1;
    Position p2;
    Position p3;

    double length;
    double numDivisions;
    std::vector<double> arcLengths;
};

CubicBezier::CubicBezier(Position p0, Position p1, Position p2, Position p3, double numDivisions)
{
    this->p0 = p0;
    this->p1 = p1;
    this->p2 = p2;
    this->p3 = p3;

    this->length = getLength();
    this->numDivisions = numDivisions;
}

Position CubicBezier::getPoint(double t)
{
    double x = pow(1 - t, 3) * p0.x + 3 * pow(1 - t, 2) * t * p1.x + 3 * (1 - t) * pow(t, 2) * p2.x + pow(t, 3) * p3.x;
    double y = pow(1 - t, 3) * p0.y + 3 * pow(1 - t, 2) * t * p1.y + 3 * (1 - t) * pow(t, 2) * p2.y + pow(t, 3) * p3.y;
    return Position(x, y, 0);
}

double CubicBezier::getLength()
{
    std::vector<Position> fakePoints;
    fakePoints.push_back(getPoint(0));
    arcLengths.push_back(0);

    double length = 0;

    for (int i = 1; i < numDivisions + 1; i++)
    {
        std::cout << "For loop " << std::endl;
        float t = (double)i / numDivisions;
        Position thisPoint = getPoint(t);
        // std::cout << "(" << thisPoint.x << ", " << thisPoint.y << ")" << std::endl;
        // std::cout << "t: " << t << std::endl;
        // std::cout << "dx = " << getPoint(t).x << " - " << fakePoints[i - 1].x << std::endl;
        // std::cout << "dy = " << getPoint(t).y << " - " << fakePoints[i - 1].y << std::endl;
        double dx = getPoint(t).x - fakePoints[i - 1].x;
        double dy = getPoint(t).y - fakePoints[i - 1].y;
        fakePoints.push_back(thisPoint);
        length += sqrt(pow(dx, 2) + pow(dy, 2));
        arcLengths.push_back(length);
        std::cout << i << "length " << length;
        std::cout << arcLengths[i] << std::endl;
        // std::cout << length << sqrt(pow(dx, 2) + pow(dy, 2)) << std::endl;
    }

    // std::cout << "length is " << length << std::endl;

    return length;
}

// double CubicBezier::getCurvature(double t)
// {
// }

double CubicBezier::getVelocity(double t)
{
    // double velX = 3 * pow(1 - t, 2) * (p1.x - p0.x) + 6 * (1 - t) * t * (p2.x - p1.x) + 3 * pow(t, 2) * (p3.x - p2.x);
    // double velY = 3 * pow(1 - t, 2) * (p1.y - p0.y) + 6 * (1 - t) * t * (p2.y - p1.y) + 3 * pow(t, 2) * (p3.y - p2.y);
    double dxdt = (-3) * pow(1 - t, 2) * p0.x + 3 * pow(1 - t, 2) * p1.x - 6 * t * (1 - t) * p1.x - 3 * pow(t, 2) * p2.x + 6 * t * (1 - t) * p2.x + 3 * pow(t, 2) * p3.x;
    double dydt = (-3) * pow(1 - t, 2) * p0.y + 3 * pow(1 - t, 2) * p1.y - 6 * t * (1 - t) * p1.y - 3 * pow(t, 2) * p2.y + 6 * t * (1 - t) * p2.y + 3 * pow(t, 2) * p3.y;
    // return sqrt(pow(velX, 2) + pow(velY, 2));
    return sqrt(pow(dxdt, 2) + pow(dydt, 2));
}

double CubicBezier::getMappedT(double t)
{
    //* binary search arc lengths to find one that matches the ratio of arc length : total length & t value
    double low = 0;
    double high = numDivisions + 1;
    int pointer = 0;
    double targetLength = length * t;
    while (low < high)
    {
        pointer = low + (high - low) / 2;
        if (arcLengths[pointer] < length * t)
        {
            low = pointer + 1;
        }
        else
        {
            high = pointer;
        }
        if (arcLengths[pointer] > length * t)
        {
            pointer--;
        }
    }
    double lengthBefore = arcLengths[pointer];
    if (lengthBefore == length * t)
    {
        return pointer / (numDivisions + 1);
    }
    else
    {
        std::cout << "Desired" << length * t << std::endl
                  << "Before" << lengthBefore << std::endl
                  << "After" << arcLengths[pointer + 1] << std::endl;
        std::cout << length * t - lengthBefore << std::endl;
        std::cout << arcLengths[pointer + 1] - lengthBefore << std::endl;
        std::cout << ((length * t) - lengthBefore) / (arcLengths[pointer + 1] - lengthBefore) << std::endl;
        return (pointer + (length * t - lengthBefore) / (arcLengths[pointer + 1] - lengthBefore)) / (numDivisions + 1);
    }
}

double CubicBezier::getAcceleration(double t)
{
    double accelX = 6 * (1 - t) * (p2.x - 2 * p1.x + p0.x) + 6 * t * (p3.x - 2 * p2.x + p1.x);
    double accelY = 6 * (1 - t) * (p2.y - 2 * p1.y + p0.y) + 6 * t * (p3.y - 2 * p2.y + p1.y);
    return sqrt(pow(accelX, 2) + pow(accelY, 2));
}

CubicBezier::CubicBezier()
{
}
class TwoDimensionalProfile
{
private:
    CubicBezier cubicBezier;
    float deltaDistance;
    std::vector<Position> calculatePoints();
    void calculateVelocities();

public:
    std::vector<Position> points;
    TwoDimensionalProfile(CubicBezier cubicBezier, float deltaDistance);
};

TwoDimensionalProfile::TwoDimensionalProfile(CubicBezier cubicBezier, float deltaDistance)
{
    this->cubicBezier = cubicBezier;
    this->deltaDistance = deltaDistance;

    this->points = calculatePoints();
}

std::vector<Position> TwoDimensionalProfile::calculatePoints()
{
    int length = cubicBezier.getLength();

    double t = 0;
    while (t < 1)
    {
        double trueT = cubicBezier.getMappedT(t);
        std::cout << " true T is " << trueT << std::endl;
        points.push_back(cubicBezier.getPoint(trueT));
        points[points.size() - 1].heading = cubicBezier.getVelocity(trueT);
        t += this->deltaDistance / length;
        std::cout << " added " << this->deltaDistance / length << std::endl;
    }
    return points;
}

int main()
{
    // CubicBezier cubicBezier(Position(1.5, 1.5, 0), Position(5, 5, 0), Position(5, 1, 0), Position(7, 6, 0));
    // CubicBezier cubicBezier(Position(0, 0, 0), Position(24, 24, 0), Position(24, -24, 0), Position(50, 0, 0));
    CubicBezier cubicBezier(Position(0, 0, 0), Position(24, 24, 0), Position(24, -24, 0), Position(50, 0, 0));
    // std::cout << cubicBezier.getPoint(NUMBER).x << ", " << cubicBezier.getPoint(NUMBER).y << std::endl;
    // std::cout << cubicBezier.getVelocity(NUMBER) << std::endl;
    TwoDimensionalProfile twoDimensionalProfile(cubicBezier, 1);
    for (int i = 0; i < twoDimensionalProfile.points.size(); i++)
    {
        // std::cout << twoDimensionalProfile.points[i].x << ", " << twoDimensionalProfile.points[i].y << ", " << twoDimensionalProfile.points[i].heading << std::endl;
        std::cout << twoDimensionalProfile.points[i].x << ", " << twoDimensionalProfile.points[i].y << std::endl;
    }
}