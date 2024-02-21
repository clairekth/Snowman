#ifndef UTILS_HPP
#define UTILS_HPP

#include "geometry.hpp"
#include "main.hpp"
#include <vector>
#define M_PI 3.14159265358979323846

/******************************************************************
 * DEFINE                                                         *
 ******************************************************************/
#define RAY_EPSILON 0.001f
#define RAY_MAX_DIST 20.f
#define RAY_ITERATIONS 512
#define WIDTH 800
#define HEIGHT 600

/******************************************************************
 * FUNCTIONS                                                      *
 ******************************************************************/
template <typename T>
inline T lerp(const T &v0, const T &v1, float t)
{
    return v0 + (v1 - v0) * std::max(0.f, std::min(1.f, t));
}

/**
 * @brief Function to compute the hash of a float
 *
 * @param n Float to compute the hash from
 * @return float Hash of the float
 */
float hash(const float n)
{
    float x = sin(n) * 43758.5453f;
    return x - floor(x);
}

/**
 * @brief Function to compute the noise of a vector
 *
 * @param x Vector to compute the noise from
 * @return float Noise of the vector
 */
float noise(const Vec3f &x)
{
    Vec3f p(floor(x.x), floor(x.y), floor(x.z));
    Vec3f f(x.x - p.x, x.y - p.y, x.z - p.z);
    f = f * (f * (Vec3f(3.f, 3.f, 3.f) - f * 2.f));
    float n = p * Vec3f(1.f, 57.f, 113.f);
    return lerp(lerp(
                    lerp(hash(n + 0.f), hash(n + 1.f), f.x),
                    lerp(hash(n + 57.f), hash(n + 58.f), f.x), f.y),
                lerp(
                    lerp(hash(n + 113.f), hash(n + 114.f), f.x),
                    lerp(hash(n + 170.f), hash(n + 171.f), f.x), f.y),
                f.z);
}

/**
 * @brief Function to rotate a vector
 *
 * @param v Vector to rotate
 * @return Vec3f Rotated vector
 */
Vec3f rotate(const Vec3f &v)
{
    return Vec3f(Vec3f(0.00, 0.80, 0.60) * v, Vec3f(-0.80, 0.36, -0.48) * v, Vec3f(-0.60, -0.48, 0.64) * v);
}

/**
 * @brief Function to compute the fractal brownian motion of a vector
 *
 * @param x Vector to compute the fractal brownian motion from
 * @return float Fractal brownian motion of the vector
 */
float fractal_brownian_motion(const Vec3f &x)
{
    Vec3f p = rotate(x);
    float f = 0;
    f += 0.5000 * noise(p);
    p = p * 2.32;
    f += 0.2500 * noise(p);
    p = p * 3.03;
    f += 0.1250 * noise(p);
    p = p * 2.61;
    f += 0.0625 * noise(p);
    return f / 0.9375;
}

void show_progress(float progress, int width, int height)
{
    int percent = int(ceil(float(progress / (width * height)) * 100.f));
    std::cout << "\r\033[K"; // Erase the entire current line.
    std::cout << "progress : [";
    if (percent % 5 == 0)
    {
        for (int i = 0; i < percent / 5; i++)
        {
            std::cout << "#";
        }
        for (int i = 0; i < 20 - percent / 5; i++)
        {
            std::cout << " ";
        }
        std::cout << "] " << percent << "%" << std::flush;
    }
}

/******************************************************************
 * STRUCTS                                                        *
 ******************************************************************/
struct Camera
{
    Vec3f pos;
    Vec3f dir;

    Camera(const Vec3f &pos, const Vec3f &dir) : pos(pos), dir(dir) {}
};

/******************************************************************
 * CLASSES                                                        *
 ******************************************************************/
class Entity
{
public:
    virtual float sdf(const Vec3f &pos) const = 0;
};

class Sphere : public Entity
{
public:
    Vec3f center;
    float radius;

    Sphere(const Vec3f &c, const float r) : center(c), radius(r) {}

    float sdf(const Vec3f &pos) const
    {
        return (pos - center).norm() - radius;
    }
};

class Shape
{
public:
    const std::vector<Entity *> entities;
    const Vec3f color;
    const float noise_amplitude;
    const float smin_factor;

    Shape(const std::vector<Entity *> e, const Vec3f &c, const float n, const float s) : entities(e), color(c), noise_amplitude(n), smin_factor(s) {}

    float sdf(const Vec3f &pos) const
    {
        float displacement = 0.0;
        if (noise_amplitude - 0.0001 >= 0)
            displacement = -fractal_brownian_motion(pos * 3.4) * noise_amplitude;
        float d = (entities.size() > 0) ? entities[0]->sdf(pos) : 1000.; // TODO: Change to defined value
        for (size_t i = 1; i < entities.size(); i++)
        {
            d = smin(d, entities[i]->sdf(pos), smin_factor);
        }
        return d + displacement;
    }

    float smin(float a, float b, float k) const
    {
        if (k - 0.0001 < 0)
            return std::min(a, b);
        float tmp = k - std::abs(a - b);
        float h = std::max(tmp, 0.f) / k;
        return std::min(a, b) - h * h * h * k * (1. / 6.);
    }
};

#endif // UTILS_HPP