#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include "geometry.hpp"
#include "utils.hpp"

/**
 * @brief Function to compute the distance between a point and the scene
 *
 * @param orig Origin of the ray
 * @return float Distance between the point and the scene
 */
float map(const Vec3f &orig, Vec3f *color);

/**
 * @brief Function to trace a ray and check if it intersects with the scene
 *
 * @param orig Origin of the ray
 * @param dir Direction of the ray
 * @param pos Position of the intersection
 * @param color Color of the intersection
 * @return true If the ray intersects with the scene
 * @return false If the ray does not intersect with the scene
 */
bool sphere_trace(const Vec3f &orig, const Vec3f &dir, Vec3f &pos, Vec3f &color);

/**
 * @brief Function to compute the normal of the distance field at a given position
 *
 * @param pos Position to compute the normal from
 * @return Vec3f Normal of the distance field at the given position
 */
Vec3f distance_field_normal(const Vec3f &pos);

/**
 * @brief Function to render the scene
 *
 * @param orig Origin of the ray
 * @param dir Direction of the ray
 * @return Vec3f Color of the pixel
 */
Vec3f render(const Vec3f &orig, const Vec3f &dir);