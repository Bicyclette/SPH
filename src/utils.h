#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <memory>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdint>
#include <array>

inline std::string file2String(std::string const & iFile)
{
	std::ifstream fileStream(iFile.c_str());
	if (fileStream.is_open())
	{
		std::stringstream buffer;
		buffer << fileStream.rdbuf();
		return buffer.str();
	}
	else
	{
		std::cerr << "Error: failed opening file \"" << iFile << "\"" << std::endl;
		std::exit(-1);
	}
}

inline float square(const float a) { return a * a; }

inline float cube(const float a) { return a * a * a; }

inline float clamp(const float v, const float vmin, const float vmax)
{
    if (v < vmin) return vmin;
    if (v > vmax) return vmax;
    return v;
}
