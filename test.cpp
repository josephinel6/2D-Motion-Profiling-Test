#define NUMBER 0.5
#define TRACK_WIDTH 12.5

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
    Position getVelocity(double t);
    Position getAcceleration(double t);
    double getCurvature(double t);
    double getLength();
    double getMappedT(double t);
    double length;

private:
    Position p0;
    Position p1;
    Position p2;
    Position p3;

    double numDivisions;
    std::vector<double> arcLengths;
};

CubicBezier::CubicBezier(Position p0, Position p1, Position p2, Position p3, double numDivisions)
{
    this->p0 = p0;
    this->p1 = p1;
    this->p2 = p2;
    this->p3 = p3;

    this->numDivisions = numDivisions;
    this->length = getLength();
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

    for (int i = 1; i < numDivisions + 1; i++)
    {
        // std::cout << "For loop " << std::endl;
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
        // std::cout << i << "length " << length;
        std::cout << arcLengths[i] << std::endl;
        // std::cout << length << sqrt(pow(dx, 2) + pow(dy, 2)) << std::endl;
    }

    std::cout << "length is " << length << std::endl;

    return length;
}

// double CubicBezier::getCurvature(double t)
// {
// }

Position CubicBezier::getVelocity(double t)
{
    // double velX = 3 * pow(1 - t, 2) * (p1.x - p0.x) + 6 * (1 - t) * t * (p2.x - p1.x) + 3 * pow(t, 2) * (p3.x - p2.x);
    // double velY = 3 * pow(1 - t, 2) * (p1.y - p0.y) + 6 * (1 - t) * t * (p2.y - p1.y) + 3 * pow(t, 2) * (p3.y - p2.y);
    double dxdt = (-3) * pow(1 - t, 2) * p0.x + 3 * pow(1 - t, 2) * p1.x - 6 * t * (1 - t) * p1.x - 3 * pow(t, 2) * p2.x + 6 * t * (1 - t) * p2.x + 3 * pow(t, 2) * p3.x;
    double dydt = (-3) * pow(1 - t, 2) * p0.y + 3 * pow(1 - t, 2) * p1.y - 6 * t * (1 - t) * p1.y - 3 * pow(t, 2) * p2.y + 6 * t * (1 - t) * p2.y + 3 * pow(t, 2) * p3.y;
    // return sqrt(pow(velX, 2) + pow(velY, 2));
    return Position(dxdt, dydt, sqrt(pow(dxdt, 2) + pow(dydt, 2)));
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
        // std::cout << "Current pointer: " << pointer << " Length * t: " << length * t << std::endl;
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
    // std::cout << "Pointer" << pointer << "Desired" << length * t << std::endl
    //           << "Before" << lengthBefore << std::endl
    //           << "After" << arcLengths[pointer + 1] << std::endl;
    if (lengthBefore == length * t)
    {
        return pointer / (numDivisions + 1);
    }
    else
    {
        std::cout << length * t - lengthBefore << std::endl;
        std::cout << arcLengths[pointer + 1] - lengthBefore << std::endl;
        std::cout << ((length * t) - lengthBefore) / (arcLengths[pointer + 1] - lengthBefore) << std::endl;
        return (pointer + (length * t - lengthBefore) / (arcLengths[pointer + 1] - lengthBefore)) / (numDivisions + 1);
    }
}

Position CubicBezier::getAcceleration(double t)
{
    double accelX = 6 * (1 - t) * (p2.x - 2 * p1.x + p0.x) + 6 * t * (p3.x - 2 * p2.x + p1.x);
    double accelY = 6 * (1 - t) * (p2.y - 2 * p1.y + p0.y) + 6 * t * (p3.y - 2 * p2.y + p1.y);
    return Position(accelX, accelY, sqrt(pow(accelX, 2) + pow(accelY, 2)));
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
    double maxVel, maxAccel;
    void setVelocityProfile();
    double getLinearVelocity(double distance);

    double accelEnd;
    double cruiseEnd;
    double decelEnd;

public:
    std::vector<Position> points;
    std::vector<std::vector<double>> velocities;
    TwoDimensionalProfile(CubicBezier cubicBezier, float deltaDistance, double maxVelocity, double maxAcceleration);
};

TwoDimensionalProfile::TwoDimensionalProfile(CubicBezier cubicBezier, float deltaDistance, double maxVelocity, double maxAcceleration)
{
    this->cubicBezier = cubicBezier;
    this->deltaDistance = deltaDistance;
    this->maxVel = maxVelocity;
    this->maxAccel = maxAcceleration;

    setVelocityProfile();
    this->points = calculatePoints();
}

//* https://math.stackexchange.com/questions/3276910/cubic-b%C3%A9zier-radius-of-curvature-calculation
double CubicBezier::getCurvature(double t)
{
    Position point = getPoint(t);
    Position firstDerivative = getVelocity(t);
    Position secondDerivative = getAcceleration(t);

    double numerator = secondDerivative.x * firstDerivative.y - secondDerivative.y * firstDerivative.x;
    double denominator = pow(pow(firstDerivative.x, 2) + pow(firstDerivative.y, 2), 3 / 2);

    return numerator / denominator;
}

double TwoDimensionalProfile::getLinearVelocity(double distance)
{
    if (distance <= accelEnd)
    {
        std::cout << "Accell" << sqrt(2 * maxAccel * distance) << std::endl;
        return sqrt(2 * maxAccel * distance);
    }
    else if (distance <= cruiseEnd)
    {
        std::cout << "Max Vel" << maxVel << std::endl;
        return maxVel;
    }
    else
    {
        std::cout << "decel" << sqrt(2 * maxAccel * (distance - cubicBezier.length)) << std::endl;
        return sqrt(2 * maxAccel * (cubicBezier.length - distance));
    }
}

void TwoDimensionalProfile::setVelocityProfile()
{
    double accelDist = pow(maxVel, 2) / (2 * maxAccel);
    // double accelDist = maxVel * (maxVel / maxAccel);
    // if (sqrt(2 * maxAccel * accelDist) > maxVel)
    // {
    //     maxAccel = sqrt(maxVel) / (2 * accelDist);
    // }

    if (accelDist > cubicBezier.length / 2) //* if we're out of distance we obviously need to do something
    {
        std::cout << "bezier length " << cubicBezier.length;
        // accelTime = sqrt((dist / 2) / (maxAccel / 2));
        accelDist = cubicBezier.length / 2;
        this->maxVel = sqrt(2 * maxAccel * accelDist);
    }

    //* now recalculate max velocity based on the new acceleration amt
    // this->maxVel = maxAccel * sqrt((cubicBezier.length / 2) / (maxAccel / 2));

    accelEnd = accelDist;
    decelEnd = cubicBezier.length;
    cruiseEnd = cubicBezier.length - accelDist;

    std::cout << "accel end: " << accelEnd;
    std::cout << "cruise end: " << cruiseEnd;
    std::cout << "decel end: " << decelEnd;
    std::cout << " max vel " << maxVel;
}

int sign(double number) //* returns sign (-1, 0, or 1)
{
    if (number >= 0)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

std::vector<Position> TwoDimensionalProfile::calculatePoints()
{
    int length = cubicBezier.length;

    double t = 0;
    double prevT = 0;
    double linearVelocity = 0;

    while (t < 1)
    {
        double trueT = cubicBezier.getMappedT(t);
        // std::cout << " true T is " << trueT << std::endl;

        Position point = cubicBezier.getPoint(trueT);
        // std::cout << "linear velocity" << linearVelocity << " " << cubicBezier.getCurvature(trueT) << std::endl;
        double angularVelocity = linearVelocity * (cubicBezier.getCurvature(trueT));
        // linearVelocity = sqrt(pow(linearVelocity, 2) + cubicBezier.getAcceleration(trueT).heading * 2 * length * trueT);
        linearVelocity = getLinearVelocity(length * trueT);

        double maxVelocityAroundCurve = (9.8 * fabs(cubicBezier.getCurvature(trueT)));
        std::cout << "Max around curve " << maxVelocityAroundCurve << std::endl;
        angularVelocity = sign(angularVelocity) * std::min(fabs(angularVelocity), maxVelocityAroundCurve);
        // linearVelocity = std::max(maxVelocityAroundCurve, linearVelocity);

        velocities.push_back({linearVelocity, angularVelocity});

        // double linearVelocity = sqrt(2 * cubicBezier.getAcceleration(t) * length * trueT);

        points.push_back(point);
        t += this->deltaDistance / length;
        prevT = trueT;
        // std::cout << " added " << this->deltaDistance / length << std::endl;
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
    TwoDimensionalProfile twoDimensionalProfile(cubicBezier, 1, 600, 30000);
    for (int i = 0; i < twoDimensionalProfile.points.size(); i++)
    {
        // std::cout << twoDimensionalProfile.points[i].x << ", " << twoDimensionalProfile.points[i].y << ", " << twoDimensionalProfile.points[i].heading << std::endl;
        // std::cout << twoDimensionalProfile.points[i].x << ", " << twoDimensionalProfile.points[i].y << std::endl;
        std::cout << i << ", " << twoDimensionalProfile.points[i].heading << std::endl;
    }
}