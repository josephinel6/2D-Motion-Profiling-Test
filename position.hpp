
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
