#include <iostream>
#include "vec3.hpp"
#include "color.hpp"

/**
 * @brief Width of the output image
 */
#define WIDTH 1280

/**
 * @brief Value of PI
 */
#define PI 3.14159265358979323846

/**
 * @brief Height of the output image
 */
#define HEIGHT 720

/**
 * @brief Maximum number of steps to ray march
 */
#define MAX_STEPS 1000

/**
 * @brief Maximum distance to the closest object
 */
#define MAX_DIST 1000.

/**
 * @brief Minimum distance to the closest object
 */
#define MIN_DIST 0.01

// TODO: Move to a separate file
float map(const point3 &p);

float deg_to_rad(float deg)
{
    return deg * (PI / 180);
}

/**
 * @brief Signed Distance Function for a sphere
 *
 * @param p Point to compute the distance
 * @param center Center of the sphere
 * @param radius Radius of the sphere
 * @return float Signed Distance to the sphere
 */
float sdf_sphere(const point3 &p, const point3 &center, float radius)
{
    return (p - center).length() - radius;
}

/**
 * @brief Function to compute the normal of a point in the scene
 *
 * @param p Point to compute the normal
 * @return vec3 Normal of the point
 */
vec3 normal(const point3 &p)
{
    float eps = 0.01;
    float x = map(p + vec3(eps, 0., 0.)) - map(p - vec3(eps, 0., 0.));
    float y = map(p + vec3(0., eps, 0.)) - map(p - vec3(0., eps, 0.));
    float z = map(p + vec3(0., 0., eps)) - map(p - vec3(0., 0., eps));

    return unit_vector(vec3(x, y, z));
}

/**
 * @brief Function to compute the minimum of two values with a smooth transition
 *
 * @param a First value
 * @param b Second value
 * @param k Smoothing factor
 * @return float Minimum of the two values with a smooth transition
 */
float smooth_min(float a, float b, float k)
{
    float tmp = k - std::abs(a - b);
    float h = std::max(tmp, float(0.)) / k;
    return std::min(a, b) - h * h * h * k * (1. / 6.);
}

/**
 * @brief Function to map a point to the distance to the closest object.
 *
 * This funtion contains the scene definition => All the objects and their positions
 *
 * @param p Point to map
 * @return float Signed Distance to the closest object
 */
float map(const point3 &p)
{
    float radius = 1.;
    point3 center1 = point3(.0, .5, 0.);
    point3 center2 = point3(1., .5, .1);

    float d = sdf_sphere(p, center1, radius);
    d = smooth_min(d, sdf_sphere(p, center2, radius), .1);
    d = std::min(d, float(p.y() + 1.));

    return d;
}

/**
 * @brief Function to ray march the scene
 *
 * @param origin Origin of the ray
 * @param direction Direction of the ray
 * @return float Distance to the closest object
 */
float ray_march(const point3 &origin, const vec3 &direction)
{
    float dist = 0.;
    for (int i = 0; i < MAX_STEPS; i++)
    {
        point3 p = origin + dist * direction;
        float dist_to_sdf = map(p);

        if (dist_to_sdf < MIN_DIST)
        {
            break;
        }

        dist += dist_to_sdf;

        if (dist > MAX_DIST)
            break;
    }
    return dist;
}

/**
 * @brief Function to render a pixel
 *
 * @param uv UV coordinates of the pixel
 * @return color Color of the pixel
 */
vec3 render(const vec3 &uv)
{
    color pixel_color = color(0., 0., 0.);

    // ray origin
    point3 origin = point3(.0, .0, -2.);
    // ray direction
    vec3 direction = vec3(uv.x(), uv.y(), 1.);

    // camera
    // TODO: A revoir sûrement
    float aspect_ratio = WIDTH / HEIGHT;
    float fov = 90.;
    float scale = tan(deg_to_rad(fov * .5));
    direction = unit_vector(vec3(direction.x(), direction.y(), direction.z() / aspect_ratio));
    direction = unit_vector(vec3(direction.x() * scale, direction.y() * scale, direction.z()));
    float dist = ray_march(origin, direction);

    if (dist <= MAX_DIST)
    {
        pixel_color = color(1., 1., 1.);

        // color the normal
        point3 p = origin + dist * direction;
        vec3 n = normal(p);
        vec3 light = vec3(2., 3., -5.);
        float intensity = std::min(std::max(dot(n, unit_vector(light - p)), 0.0), 1.0);

        intensity *= 5. / dot(light - p, light - p);

        // shadows
        // TODO: A revoir sûrement
        float shadow_dist = ray_march(p + n * .1, light - p);
        if (shadow_dist < (light - p).length() && shadow_dist < MAX_DIST)
        {
            intensity *= .5;
        }

        pixel_color = vec3(pow(intensity, .4545), pow(intensity, .4545), pow(intensity, .4545));
    }

    return pixel_color;
}

int main(int argc, char **argv)
{
    float img_factor = 1.;
    if (argc > 1)
        img_factor = std::stof(argv[1]);

    // Image
    int image_width = WIDTH * img_factor;
    int image_height = HEIGHT * img_factor;

    // Render

    std::cout << "P3\n"
              << image_width << ' ' << image_height << "\n255\n";
#pragma omp parallel for
    for (int j = 0; j < image_height; ++j)
    {
        std::clog << "\rScanlines remaining: " << (image_height - j) << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i)
        {
            vec3 uv = vec3(i, j, 0) / vec3(image_width, image_height, 1);
            uv = uv - vec3(0.5, 0.5, 0.);
            uv *= 2.;
            // set the (-1, 1) coordinate system to the bottom left
            uv[1] *= -1.; // TODO: Check to avoid this

            uv = vec3(uv.x() * image_width / image_height, uv.y(), uv.z());

            color pixel_color = render(uv);
            write_color(std::cout, pixel_color);
        }
    }
    std::clog << "\rDone.                 \n";
}