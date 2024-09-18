#include "geometry.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#define M_PI 3.14159265358979323846

struct Light {
    Light(const Vec3f &position, const float &intensity) : position(position), intensity(intensity) {}
    Vec3f position;
    float intensity;
};

struct Material {
    Material(const Vec3f &color) : diffuse_color(color) {}
    Material() : diffuse_color() {}
    Vec3f diffuse_color;
};

struct Sphere {
    Vec3f center;
    float radius;

    Material material;

    Sphere(
        const Vec3f &center,
        const float &radius,
        const Material &material) : center(center), radius(radius), material(material) {}

    bool ray_intersect(
        const Vec3f &origin,
        const Vec3f &direction,
        float &t0) const {

        Vec3f L = center - origin;
        float tca = L * direction;
        float d2 = L * L - tca * tca;

        if (d2 > radius * radius) {
            return false;
        }

        float thc = sqrtf(radius * radius - d2);
        t0 = tca - thc;
        float t1 = tca + thc;

        if (t0 < 0) {
            t0 = t1;
        }

        if (t0 < 0) {
            return false;
        }
        return true;
    }
};

bool scene_intersect(
    const Vec3f &origin,
    const Vec3f &direction,
    const std::vector<Sphere> &spheres,
    Vec3f &hit,
    Vec3f &N,
    Material &material) {

    float spheres_dist = std::numeric_limits<float>::max();

    for (size_t i = 0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].ray_intersect(origin, direction, dist_i) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = origin + direction * dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }

    return spheres_dist < 1000;
}

Vec3f cast_ray(
    const Vec3f &origin,
    const Vec3f &direction,
    const std::vector<Sphere> &spheres,
    const std::vector<Light> lights) {

    Vec3f point, N;
    Material material;

    if (!scene_intersect(origin, direction, spheres, point, N, material)) {
        return Vec3f(0.2, 0.7, 0.8);
    }

    float diffuse_light_intensity = 0;
    for (size_t i = 0; i < lights.size(); i++) {
        Vec3f lightDirection = (lights[i].position - point).normalize();
    }

    return material.diffuse_color;
}

void render(const std::vector<Sphere> &spheres) {
    const int width = 1024;
    const int height = 728;
    const float fov = M_PI / 2.;
    std::vector<Vec3f> framebuffer(width * height);

    for (size_t j = 0; j < height; j++) {
        for (size_t i = 0; i < width; i++) {
            float x = (2 * (i + 0.5) / (float)width - 1) * tan(fov / 2.) * width / (float)height;
            float y = -(2 * (j + 0.5) / (float)height - 1) * tan(fov / 2.);
            Vec3f direction = Vec3f(x, y, -1).normalize();

            framebuffer[i + j * width] = cast_ray(Vec3f(0, 0, 0), direction, spheres);
        }
    }

    std::ofstream ofs;
    ofs.open("./out.ppm", std::ios::binary);
    ofs << "P6\n"
        << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height * width; ++i) {
        for (size_t j = 0; j < 3; j++) {
            ofs << static_cast<char>(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

int main() {
    Material ivory(Vec3f(0.4, 0.4, 0.3));
    Material red_rubber(Vec3f(0.3, 0.1, 0.1));

    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(Vec3f(-3, 0, -16), 2, ivory));
    spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2, red_rubber));
    spheres.push_back(Sphere(Vec3f(1.5, -0.5, -18), 3, red_rubber));
    spheres.push_back(Sphere(Vec3f(7, 5, -18), 4, ivory));
    render(spheres);

    return 0;
}