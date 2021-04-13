
/**
 * Main module. 
 * This program is for calculating and saving data. For animations please use animate.cpp
 */
#include <string>
#include <sstream>
#include <fstream>

#include "Body.h"

std::string filename(int n);

int main(int argc, char *argv[])
{

    return 0;
}

std::string filename(int n)
{
    std::string name;
    std::stringstream n_s;
    n_s << n;
    name = "./" + n_s.str() + ".csv";
    return name;
}