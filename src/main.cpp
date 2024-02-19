#include "main.hpp"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define RAY_EPSILON 0.000001f
#define RAY_MAX_DIST 1000.f
#define RAY_ITERATIONS 1024

#define WIDTH 640
#define HEIGHT 480

template <typename T>
inline T lerp(const T &v0, const T &v1, float t)
{
    return v0 + (v1 - v0) * std::max(0.f, std::min(1.f, t));
}

float hash(const float n)
{
    float x = sin(n) * 43758.5453f;
    return x - floor(x);
}

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

Vec3f rotate(const Vec3f &v)
{
    return Vec3f(Vec3f(0.00, 0.80, 0.60) * v, Vec3f(-0.80, 0.36, -0.48) * v, Vec3f(-0.60, -0.48, 0.64) * v);
}

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

float sdf_sphere(const Vec3f &p, const Vec3f &center, float radius)
{
    return (p - center).norm() - radius;
}

float sdf_sphere_displaced(const Vec3f &p, const Vec3f &center, float radius, float noise_amplitude)
{
    float displacement = (sin(16 * p.x) * sin(16 * p.y) * sin(16 * p.z) + 1.) / 2.;
    return (p - center).norm() - radius + displacement * noise_amplitude;
}

float sdf_sphere_displaced_noised(const Vec3f &p, const Vec3f &center, float radius, float noise_amplitude)
{
    float displacement = -fractal_brownian_motion(p * 3.4) * noise_amplitude;
    return (p - center).norm() - radius - displacement;
}

float smooth_min(float a, float b, float k)
{
    float tmp = k - std::abs(a - b);
    float h = std::max(tmp, float(0.)) / k;
    return std::min(a, b) - h * h * h * k * (1. / 6.);
}

float map(const Vec3f &orig)
{
    Vec3f sphere1(0, 0, -2);
    float radius1 = 1.5;

    Vec3f sphere2(0, 1., -2);
    float radius2 = 1.5;

    float d = sdf_sphere_displaced_noised(orig, sphere1, radius1, .1);
    d = smooth_min(d, sdf_sphere(orig, sphere2, radius2), 0.2);

    return d;
}

bool sphere_trace(const Vec3f &orig, const Vec3f &dir, Vec3f &pos)
{
    pos = orig;
    for (int i = 0; i < RAY_ITERATIONS; i++)
    {
        float d = map(pos);
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
    Vec3f hit;
    if (!sphere_trace(orig, dir, hit))
        return Vec3f(0.2, 0.7, 0.8);
    Vec3f light_dir = (Vec3f(10, 10, 10) - hit).normalize(); // light is placed at (10,10,10)
    float light_intensity = std::max(0.4f, light_dir * distance_field_normal(hit));

    Vec3f color;
    if (map(hit) == sdf_sphere_displaced_noised(hit, Vec3f(0, 0, -2), 1.5, .1))
    {
        color = Vec3f(1, 0, 0); // couleur de la première sphère
    }
    else
    {
        color = Vec3f(0, 1, 0); // couleur de la deuxième sphère
    }

    return color * light_intensity;
}

int main(int argc, char **argv)
{
    int factor = 1;
    if (argc > 1)
        factor = std::stoi(argv[1]);
    const int width = WIDTH * factor, height = HEIGHT * factor;
    const float fov = M_PI / 3.;
    float progess = 0;
    std::vector<Vec3f> framebuffer(width * height);
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
    for (int j = 0; j < height; j++)
    { // actual rendering loop
        // TODO: Add a progress bar
        int percent = int(ceil(float(progess / (width * height)) * 100.f));
        std::cout << "\r\033[K"; // Erase the entire current line.
        std::cout << "progess : [";
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
            framebuffer[i + j * width] = render(origin, dir);
            progess++;
        }
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::vector<unsigned char> pixmap(width * height * 3);
    for (int i = 0; i < height * width; ++i)
    {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max > 1)
            c = c * (1. / max);
        for (int j = 0; j < 3; j++)
        {
            pixmap[i * 3 + j] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    stbi_write_jpg("out.jpg", width, height, 3, pixmap.data(), 100);
    std::chrono::duration<double> time_span = t2 - t1;
    std::cout << "\noutput generated in " << time_span.count() << " seconds." << std::endl;

    return 0;
}
