#include <iostream>
#include "simplex.h"
#include "simplex_implementation.h"

int main()
{
    std::vector<double> init_coords = { 1.5, 2 };
    Point p = Point(init_coords);

    Simplex simplex = Simplex(p, 2);
    Point result = simplex.optimize();

    for (int i = 0; i < result.getCoords().size(); ++i)
        std::cout << result.getCoords()[i] << ' ';
    return 0;
}