#include "main.hpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

float lighting(const Vec3f &p)
{
    Vec3f light_pos = Vec3f(-5, 5, 2);             // Position of the light source
    Vec3f light_dir = (light_pos - p).normalize(); // Direction of the light source from the hit point

    Vec3f n = distance_field_normal(p); // Normal at the hit point

    return std::max(0.f, n.dot(light_dir));
}

float map(const Vec3f &orig, Vec3f *color = nullptr)
{
    Sphere sphere_top(Vec3f(0.0, 1.55, -5.0), .65);
    Sphere sphere_middle(Vec3f(0.0, 0.3, -5.0), 1.0);
    Sphere sphere_bottom(Vec3f(0.0, -1.5, -5.0), 1.5);
    Shape snow({&sphere_top, &sphere_middle, &sphere_bottom}, Vec3f(1, 1., 1), 0.1, 0.1);

    Sphere right_eye(Vec3f(0.25, 1.8, -4.5), 0.12);
    Sphere left_eye(Vec3f(-0.25, 1.8, -4.5), 0.12);
    Shape eyes({&right_eye, &left_eye}, Vec3f(0, 0, 0), 0., 0.);

    Sphere button1(Vec3f(0.0, 0.5, -3.0), 0.1);
    Sphere button2(Vec3f(0.0, 0, -3.0), 0.1);
    Sphere button3(Vec3f(0.0, -0.4, -0.5), 0.08);
    Sphere button4(Vec3f(0.0, -.8, -.5), 0.08);
    Shape buttons({&button1, &button2, &button3, &button4}, Vec3f(0.05, 0.05, 0.05), 0., 0.);

    Cylinder hat_base(Vec3f(0.0, 2.1, -5.0), 1, 0.01);
    Cylinder hat_top(Vec3f(0.0, 2.6, -5.0), 0.5, 1.);
    Shape hat({&hat_base, &hat_top}, Vec3f(0.5, 0.5, 0.5), 0., 0.);

    Cylinder ribbon(Vec3f(0.0, 2.6, -5.0), 0.52, 0.2);
    Shape ribbons({&ribbon}, Vec3f(0.8, 0.5, 0.5), 0., 0.);

    Capsule left_arm(Vec3f(-0.7, 0.3, -4.5), Vec3f(-2.5, -0.3, -4.2), 0.05);
    Capsule left_finger_1(Vec3f(-1.7, -0.05, -4.3), Vec3f(-2.2, 0.29, -4.25), 0.05);
    Capsule left_finger_2(Vec3f(-2.1, -0.20, -4.25), Vec3f(-2.35, -0.55, -4.18), 0.05);

    Capsule right_arm(Vec3f(.7, 0.3, -4.5), Vec3f(2.5, -0.3, -4.2), 0.05);
    Capsule right_finger_1(Vec3f(1.7, -0.05, -4.3), Vec3f(2.2, 0.29, -4.25), 0.05);
    Capsule right_finger_2(Vec3f(2.1, -0.20, -4.26), Vec3f(2.35, -0.55, -4.18), 0.05);

    Shape arms({&left_arm, &right_arm, &left_finger_1, &left_finger_2, &right_finger_1, &right_finger_2}, Vec3f(0.467, 0.275, 0.204), 0.1, 0.001);

    Sphere nose_base(Vec3f(0.0, 1.50, -4.3), 0.15);
    Sphere nose_top(Vec3f(0.0, 1.50, -3.9), 0.02);
    Shape nose({&nose_base, &nose_top}, Vec3f(1, 0.5, 0.31), 0.0001, 0.85);

    std::vector<Shape> shapes = {nose, snow, eyes, arms, buttons, hat, ribbons};

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

bool sphere_trace(Camera &camera, Vec3f &pos, Vec3f &color)
{
    pos = camera.pos;
    for (int i = 0; i < RAY_ITERATIONS; i++)
    {
        float d = map(pos, &color);
        if (d < RAY_EPSILON)
            return true;
        pos = pos + camera.dir * d;
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

Vec3f render(Camera &camera, std::vector<Vec3f> *env_map, const bool is_loaded, const int env_map_width, const int env_map_height)
{
    Vec3f color = Vec3f(1, 1, 1); // Default color
    Vec3f hit;
    if (!sphere_trace(camera, hit, color) && is_loaded)
    {
        // FIXME: This is not working
        // float theta = acosf(camera.dir.y) / M_PI;
        // float phi = (atan2f(camera.dir.z, camera.dir.x) + M_PI) / (2 * M_PI);
        // int x = std::min(int(phi * env_map_width), env_map_width - 1);
        // int y = std::min(int(theta * env_map_height), env_map_height - 1);
        int x = ((atan2(camera.dir.z, camera.dir.x) + 5) / (2 * M_PI)) * env_map_width;
        int y = ((acos(camera.dir.y) / M_PI)) * env_map_height;

        color = (*env_map)[x + y * env_map_width];
    }
    else
    {
        color = color * lighting(hit);
    }

    return color;
}

bool load_env_map(std::vector<Vec3f> *env_map, int *env_map_width, int *env_map_height)
{
    int width, height, channels;
    unsigned char *data = stbi_load("../envmap.jpg", &width, &height, &channels, 0);
    if (!data)
    {
        std::cerr << "Error: failed to load texture" << std::endl;
        return false;
    }
    env_map->resize(width * height);
    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            int idx = (i + j * width) * channels;
            (*env_map)[i + j * width] = Vec3f(data[idx] / 255.f, data[idx + 1] / 255.f, data[idx + 2] / 255.f);
        }
    }
    stbi_image_free(data);
    *env_map_width = width;
    *env_map_height = height;
    return true;
}

int main(int argc, char **argv)
{
    int factor = 1;
    if (argc > 1)
        factor = std::stoi(argv[1]);
    const int width = WIDTH * factor, height = HEIGHT * factor;
    float progress = 0;
    std::vector<Vec3f> framebuffers(width * height);
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    std::vector<Vec3f> env_map;
    int env_map_width, env_map_height;
    bool is_loaded = load_env_map(&env_map, &env_map_width, &env_map_height);
#pragma omp parallel for
    for (int j = 0; j < height; j++)
    { // actual rendering loop
        show_progress(progress, width, height);
        for (int i = 0; i < width; i++)
        {
            // Normalize point to [-1, 1]
            float x = i / (float)width * 2 - 1;
            float y = j / (float)height * 2 - 1;
            y *= -height / (float)width; // Swap y to match the orientation of the image and correct the aspect ratio

            const float fov = 1.;
            Vec3f origin(0., 1.5, 1.);
            Vec3f target(0.0, 0.0, -5.0);

            Camera camera(origin, target, fov, x, y);
            framebuffers[i + j * width] = render(camera, &env_map, is_loaded, env_map_width, env_map_height);
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
