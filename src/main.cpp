#include "main.hpp"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


float sdf_sphere(const Vec3f &p, const Vec3f &center, float radius)
{
    return (p - center).norm() - radius;
}

float sdf_sphere_displaced(const Vec3f &p, const Vec3f &center, float radius, float noise_amplitude)
{
    float displacement = (sin(16 * p.x) * sin(16 * p.y) * sin(16 * p.z) + 1.) / 2.;
    return (p - center).norm() - radius + displacement * noise_amplitude;
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

    float d = sdf_sphere_displaced(orig, sphere1, radius1, 0.1);
    d = smooth_min(d, sdf_sphere(orig, sphere2, radius2), 0.2);

    return d;
}

bool sphere_trace(const Vec3f &orig, const Vec3f &dir, Vec3f &pos)
{
    pos = orig;
    for (size_t i = 0; i < 128; i++)
    {
        float d = map(pos);
        if (d < 0.)
            return true;
        pos = pos + dir * std::max(d * 0.1f, .01f);
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

    return Vec3f(1, 1, 1) * light_intensity;
}

int main(int argc, char **argv)
{
    int factor = 1;
    if (argc > 1)
        factor = std::stoi(argv[1]);
    const int width = 640 * factor, height = 480 * factor;
    const float fov = M_PI / 3.;
    std::vector<Vec3f> framebuffer(width * height);
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
    for (int j = 0; j < height; j++)
    { // actual rendering loop
        // TODO: Add a progress bar
        // FIXME: Is printed in random order
        // std::cout << height - j << " lines remaining\n" << std::flush;
        for (int i = 0; i < width; i++)
        {
            Vec3f origin(0, 0, 3);
            float dir_x = (i + 0.5) - width / 2.;
            float dir_y = -(j + 0.5) + height / 2.; // this flips the image at the same time
            float dir_z = -height / (2. * tan(fov / 2.));
            Vec3f dir = Vec3f(dir_x, dir_y, dir_z).normalize();
            framebuffer[i + j * width] = render(origin, dir);
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
    // TODO: Add output to png with stb_image : stbi_write_png --> need a pixmap array instead of a framebuffer
    std::cout << "output generated in " << time_span.count() << " seconds." << std::endl;

    return 0;
}
