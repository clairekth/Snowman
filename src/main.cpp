#include "main.hpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define RAY_EPSILON 0.001f
#define RAY_MAX_DIST 20.f
#define RAY_ITERATIONS 512

// #define RAY_EPSILON 0.0000001f
// #define RAY_MAX_DIST 10000.f
// #define RAY_ITERATIONS 2048

#define WIDTH 800
#define HEIGHT 600


float lighting(const Vec3f &p)
{
    Vec3f light_pos = Vec3f(-5, 5, 2);             // Position of the light source
    Vec3f light_dir = (light_pos - p).normalize(); // Direction of the light source from the hit point

    Vec3f n = distance_field_normal(p); // Normal at the hit point

    return std::max(0.f, n.dot(light_dir));
}

float map(const Vec3f &orig, Vec3f *color = nullptr)
{
    Sphere sphere1(Vec3f(0.0, 1.0, -5.0), 1.0);
    Sphere sphere2(Vec3f(0.0, 0.0, -5.0), 1.0);
    Sphere sphere3(Vec3f(0.0, -1.0, -5.0), 1.0);

    Shape shape({&sphere1, &sphere2, &sphere3}, Vec3f(1.0, 0.5, 1.0), 0.0, 0.0);

    Sphere sphere4(Vec3f(-1.0, 1.0, -5.0), 0.7);
    Sphere sphere5(Vec3f(1.0, 1.0, -5.0), 0.7);

    Shape shape2({&sphere4, &sphere5}, Vec3f(0.5, 1.0, 0.5), 0.1, 0.1);

    std::vector<Shape> shapes = {shape, shape2};

    float d = RAY_MAX_DIST;
    for (size_t i = 0; i < shapes.size(); i++)
    {
        float d_tmp = shapes[i].sdf(orig);
        if (d_tmp < d)
        {
            d = d_tmp;
            if (color)
                *color = shapes[i].color;
        }
    }
    return d;
}

bool sphere_trace(const Vec3f &orig, const Vec3f &dir, Vec3f &pos, Vec3f &color)
{
    pos = orig;
    for (int i = 0; i < RAY_ITERATIONS; i++)
    {
        float d = map(pos, &color);
        if (d < RAY_EPSILON)
            return true;
        pos = pos + dir * std::max(d * 0.1f, .01f);
        if (d > RAY_MAX_DIST)
            break;
    }
    return false;
}

Vec3f distance_field_normal(const Vec3f &pos)
{
    const float eps = 0.1;
    float d = map(pos);
    float nx = map(pos + Vec3f(eps, 0, 0)) - d;
    float ny = map(pos + Vec3f(0, eps, 0)) - d;
    float nz = map(pos + Vec3f(0, 0, eps)) - d;
    return Vec3f(nx, ny, nz).normalize();
}

Vec3f render(const Vec3f &orig, const Vec3f &dir)
{
    Vec3f color = Vec3f(1, 1, 1); // Default color
    Vec3f hit;
    if (!sphere_trace(orig, dir, hit, color))
        return Vec3f(0.2, 0.7, 0.8);

    color = color * lighting(hit);

    return color;
}

int main(int argc, char **argv)
{
    int factor = 1;
    if (argc > 1)
        factor = std::stoi(argv[1]);
    const int width = WIDTH * factor, height = HEIGHT * factor;
    const float fov = M_PI / 3.;
    float progress = 0;
    std::vector<Vec3f> framebuffers(width * height);
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
    for (int j = 0; j < height; j++)
    { // actual rendering loop
        // TODO: Add a progress bar
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
        for (int i = 0; i < width; i++)
        {
            Vec3f origin(0, 0, 3);
            float dir_x = (i + 0.5) - width / 2.;
            float dir_y = -(j + 0.5) + height / 2.; // this flips the image at the same time
            float dir_z = -height / (2. * tan(fov / 2.));
            Vec3f dir = Vec3f(dir_x, dir_y, dir_z).normalize();
            framebuffers[i + j * width] = render(origin, dir);
            progress++;
        }
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::vector<unsigned char> pixmap(width * height * 3);
    for (int i = 0; i < height * width; ++i)
    {
        Vec3f &c = framebuffers[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max > 1)
            c = c * (1. / max);
        for (int j = 0; j < 3; j++)
        {
            pixmap[i * 3 + j] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffers[i][j])));
        }
    }
    stbi_write_jpg("out.jpg", width, height, 3, pixmap.data(), 100);
    std::chrono::duration<double> time_span = t2 - t1;
    std::cout << "\noutput generated in " << time_span.count() << " seconds." << std::endl;

    return 0;
}
