#pragma once
#include <vector>
#include <random>
#include <algorithm>
#include <map>
#include <string>

class Point
{
private:
    std::vector<double> coords;
public:
    Point(std::vector<double> coords_) : coords(coords_) {}

    std::vector<double> getCoords()
    {
        return coords;
    }

    Point operator+(const Point& other_p)
    {
        if (other_p.coords.size() != coords.size())
            throw "The points demensionality must match!";

        std::vector<double> new_p_coords(coords.size());
        for (int i = 0; i < coords.size(); ++i)
            new_p_coords[i] = other_p.coords[i] + coords[i];

        return Point(new_p_coords);
    }

    Point operator*(const Point& other_p)
    {
        if (other_p.coords.size() != coords.size())
            throw "The points demensionality must match!";

        std::vector<double> new_p_coords(coords.size());
        for (int i = 0; i < coords.size(); ++i)
            new_p_coords[i] = other_p.coords[i] * coords[i];

        return Point(new_p_coords);
    }

    Point operator*(const double n)
    {
        std::vector<double> new_p_coords(coords.size());
        for (int i = 0; i < coords.size(); ++i)
            new_p_coords[i] = coords[i] * n;

        return Point(new_p_coords);
    }

    Point operator/(const double n)
    {
        std::vector<double> new_p_coords(coords.size());
        for (int i = 0; i < coords.size(); ++i)
            new_p_coords[i] = coords[i] / n;

        return Point(new_p_coords);
    }

    Point operator-(const Point& other_p)
    {
        if (other_p.coords.size() != coords.size())
            throw "The points demensionality must match!";

        std::vector<double> new_p_coords(coords.size());
        for (int i = 0; i < coords.size(); ++i)
            new_p_coords[i] = coords[i] - other_p.coords[i];

        return Point(new_p_coords);
    }

    double len()
    {
        double s = 0;
        for (int i = 0; i < coords.size(); ++i)
            s += std::powf(coords[i], 2);
        return std::sqrtf(s);
    }
};

class Simplex
{
private:
    std::vector<Point> simplex_points;
    double alpha = 1, beta = 0.5, gamma = 2;
    double threshold = 0.001;
public:
    Simplex(Point, int, double);
    double func(Point);
    std::map<std::string, double> findExtPoints();
    bool didConv();
    Point optimize();
};