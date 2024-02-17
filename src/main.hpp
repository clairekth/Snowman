#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include "geometry.hpp"

/**
 * @brief Function to compute the signed distance between a point and a sphere
 * 
 * @param p Point to compute the distance from
 * @param center Center of the sphere
 * @param radius Radius of the sphere
 * @return float Signed distance between the point and the sphere
 */
float sdf_sphere(const Vec3f &p, const Vec3f &center, float radius);

/**
 * @brief Function to compute the signed distance between a point and a sphere with a displacement
 * 
 * @param p Point to compute the distance from
 * @param center Center of the sphere
 * @param radius Radius of the sphere
 * @param noise_amplitude Amplitude of the noise
 * @return float Signed distance between the point and the sphere
 */
float sdf_sphere_displaced(const Vec3f &p, const Vec3f &center, float radius, float noise_amplitude);

/**
 * @brief Function to compute the minimum of two values with a smoothing factor
 * 
 * @param a First value
 * @param b Second value
 * @param k Smoothing factor
 * @return float Minimum of the two values with a smoothing factor
 */
float smooth_min(float a, float b, float k);

/**
 * @brief Function to compute the distance between a point and the scene
 * 
 * @param orig Origin of the ray
 * @return float Distance between the point and the scene
 */
float map(const Vec3f &orig);

/**
 * @brief Function to trace a ray and check if it intersects with the scene
 * 
 * @param orig Origin of the ray
 * @param dir Direction of the ray
 * @param pos Position of the intersection
 * @return true If the ray intersects with the scene
 * @return false If the ray does not intersect with the scene
 */
bool sphere_trace(const Vec3f &orig, const Vec3f &dir, Vec3f &pos);

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