#pragma once
#include "simplex.h"

Simplex::Simplex(Point initial_point, int n, double thrs = 0.00001)
{
    /* During creating of a simplex object there is need to declare an initial point.
    * Pleas, chose it carefully, depending on it you can sink into local extremum.
    * Other simplex's points are calculated from the inital point
    * 
    * Params:
    *       initial_point - what the fuck didn't you understand?
    *       n - the demension of SIMPLEX (not a function)
    */

    threshold = thrs;

    double k = 0.5; // size of simplex side

    this->simplex_points.push_back(initial_point);

    for (int i = 1; i <= n; ++i)
    {
        std::vector<double> new_coords(n);
        new_coords[i - 1] = 1 * k;

        for (int j = 0; j < n; ++j)
            new_coords[j] += initial_point.getCoords()[j];

        Point new_p = Point(new_coords);
        this->simplex_points.push_back(new_p);
    }
}

double Simplex::func(Point p)
{
    /* This funcion is an implementation of source func which we need to optimize.
    * Parameter:
    *       p - point where the function will be calculated
    */

    return std::powf(p.getCoords()[0] + 1, 2) + std::powf(p.getCoords()[1], 2);
}

std::map<std::string, double> Simplex::findExtPoints()
{
    std::map<std::string, double> result{ {"f_h", -1000000}, {"f_g", -1000000}, {"f_l", 1000000},
        {"x_h_index",  0}, {"x_g_index", 0}, { "x_l_index", 0 }
    };

    for (int i = 0; i < simplex_points.size(); ++i)
    {
        double func_val = func(simplex_points[i]);
        if (func_val > result["f_h"])
        {
            result["f_g"] = result["f_h"];
            result["x_g_index"] = result["x_h_index"];
            result["f_h"] = func_val;
            result["x_h_index"] = i;
        }
        else if (func_val > result["f_g"]) // Updated this condition
        {
            result["f_g"] = func_val;
            result["x_g_index"] = i;
        }
        if (func_val < result["f_l"])
        {
            result["f_l"] = func_val;
            result["x_l_index"] = i;
        }
    }

    return result;
}


bool Simplex::didConv()
{
    double f_mean = 0;
    for (int i = 0; i < simplex_points.size(); ++i)
        f_mean += func(simplex_points[i]);
    f_mean = f_mean / simplex_points.size();

    double sigma = 0;
    for (int i = 0; i < simplex_points.size(); ++i)
        sigma += std::powf(func(simplex_points[i]) - f_mean, 2);
    sigma = std::sqrtf(sigma / simplex_points.size());

    return (sigma < threshold) ? true : false;
}

Point Simplex::optimize()
{
    int iterations = 0;
    while (true)
    {
        iterations++;
        // find x which are: the biggest func, the next biggest func and the lowest func
        std::map<std::string, double> ext_simplex_points = findExtPoints();

        // calculate center of mass
        std::vector<double> x0(simplex_points[0].getCoords().size());
        Point x0_ = Point(x0);
        for (int i = 0; i < simplex_points.size() - 1; ++i)
        {
            if (i != ext_simplex_points["x_h_index"])
                x0_ = x0_ + simplex_points[i];
        }
        x0_ = x0_ / (simplex_points.size() - 1);
        double f_x0 = func(x0_);

        // find a point of reflection
        Point xr = x0_ * (1 + alpha) - simplex_points[ext_simplex_points["x_h_index"]] * alpha;
        alpha = (xr - x0_).len() / (x0_ - simplex_points[ext_simplex_points["x_h_index"]]).len();
        double f_xr = func(xr);

        if (f_xr < ext_simplex_points["f_l"]) // extent the simplex
        {
            Point xe = xr * gamma + x0_ * (1 - gamma);
            gamma = (xe - x0_).len() / (xr - x0_).len();
            double f_xe = func(xe);

            if (f_xe < ext_simplex_points["f_l"])
                simplex_points[ext_simplex_points["x_h_index"]] = xe;
            else
                simplex_points[ext_simplex_points["x_h_index"]] = xr;
            if (didConv())
            {
                std::cout << "Iterations: " << iterations << '\n';
                std::cout << "Func value:" << func(simplex_points[ext_simplex_points["x_l_index"]]) << '\n';
                return simplex_points[ext_simplex_points["x_l_index"]];
            }
            continue;
        }
        else if (f_xr <= ext_simplex_points["f_g"])
        {
            simplex_points[ext_simplex_points["x_h_index"]] = xr;
            if (didConv())
            {
                std::cout << "Iterations: " << iterations << '\n';
                std::cout << "Func value:" << func(simplex_points[ext_simplex_points["x_l_index"]]) << '\n';
                return simplex_points[ext_simplex_points["x_l_index"]];
            }
            continue;
        }
        else if (f_xr > ext_simplex_points["f_g"]) // squeeze simplex
        {
            if (f_xr <= ext_simplex_points["f_h"])
            {
                simplex_points[ext_simplex_points["x_h_index"]] = xr;
                ext_simplex_points["f_h"] = f_xr;
            }

            Point xc = simplex_points[ext_simplex_points["x_h_index"]] * beta + x0_ * (1 - beta);
            beta = (xc - x0_).len() / (xr - x0_).len();
            double f_xc = func(xc);

            if (f_xc <= ext_simplex_points["f_h"])
            {
                simplex_points[ext_simplex_points["x_h_index"]] = xc;
                if (didConv())
                {
                    std::cout << "Iterations: " << iterations << '\n';
                    std::cout << "Func value:" << func(simplex_points[ext_simplex_points["x_l_index"]]) << '\n';
                    return simplex_points[ext_simplex_points["x_l_index"]];
                }
                continue;
            }
            else
            {
                for (int i = 0; i < simplex_points.size(); ++i)
                    if (i != ext_simplex_points["x_l_index"])
                        simplex_points[i] = (simplex_points[i] + simplex_points[ext_simplex_points["x_l_index"]]) * 0.5;
                if (didConv())
                {
                    std::cout << "Iterations: " << iterations << '\n';
                    std::cout << "Func value:" << func(simplex_points[ext_simplex_points["x_l_index"]]) << '\n';
                    return simplex_points[ext_simplex_points["x_l_index"]];
                }
                continue;
            }
        }
    }
}
