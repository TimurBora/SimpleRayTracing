#include "geometry.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#define M_PI 3.14159265358979323846

enum LightType {
    Ambient,
    Point,
    Directional,
};

struct AmbientLight {
    const float &intensity;
};

struct PointLight {
    const Vec3f &position;
    const float &intensity;
};

struct DirectionalLight {
    const Vec3f &direction;
    const float &intensity;
};

// using LightEnum = std::variant<AmbientLight, PointLight, DirectionalLight>;

struct Light {
    Light(LightType lightType, const float &intensity)
        : lightType(lightType), intensity(intensity), direction(), position() {
        if (lightType != Ambient) {
            throw std::invalid_argument("This constructor is only for Ambient lights");
        }
    }

    Light(LightType lightType, const Vec3f &position, const float &intensity)
        : lightType(lightType), position(position), intensity(intensity), direction() {
        if (lightType != Point) {
            throw std::invalid_argument("This constructor is only for Point lights");
        }
    }

    // Light(LightType lightType, const Vec3f &direction, const float &intensity)
    //     : lightType(lightType), direction(direction), intensity(intensity), position() {
    //     if (lightType != Directional) {
    //         throw std::invalid_argument("This constructor is only for Directional lights");
    //     }
    // }

    LightType lightType;
    Vec3f position;
    Vec3f direction;
    float intensity;
};

struct Material {
    Material(
        const Vec4f &a,
        const Vec3f &color,
        const Vec3f &ambient,
        const float &spec,
        const float &rf) : albedo(a), diffuseColor(color), specularExponent(spec), refractiveIndex(rf), ambientColor(ambient) {}
    Material() : albedo(1, 0, 0, 0), diffuseColor(), specularExponent(), refractiveIndex(1.0), ambientColor(0.1, 0.1, 0.1) {}
    Vec4f albedo;
    Vec3f diffuseColor;
    Vec3f ambientColor;
    float specularExponent;
    float refractiveIndex;
};

struct Sphere {
    Vec3f center;
    float radius;

    Material material;

    Sphere(
        const Vec3f &center,
        const float &radius,
        const Material &material) : center(center), radius(radius), material(material) {}

    bool RayIntersect(
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

Vec3f Reflect(const Vec3f &direction, const Vec3f &N) {
    return direction - N * 2.f * (direction * N);
}

Vec3f Refract(const Vec3f &direction, const Vec3f &N, const float &refractiveIndex) {
    float cosi = -std::max(-1.f, std::min(1.f, direction * N));
    float etai = 1, etat = refractiveIndex;
    Vec3f n = N;
    if (cosi < 0) {
        cosi *= -1;
        std::swap(etai, etat);
        n = -N;
    }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? Vec3f(0, 0, 0) : direction * eta + n * (eta * cosi - sqrtf(k));
}

bool SceneIntersect(
    const Vec3f &origin,
    const Vec3f &direction,
    const std::vector<Sphere> &spheres,
    Vec3f &hit,
    Vec3f &N,
    Material &material) {

    float spheresDist = std::numeric_limits<float>::max();

    for (size_t i = 0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].RayIntersect(origin, direction, dist_i) && dist_i < spheresDist) {
            spheresDist = dist_i;
            hit = origin + direction * dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }

    return spheresDist < 1000;
}

Vec3f CastRay(
    const Vec3f &origin,
    const Vec3f &direction,
    const std::vector<Sphere> &spheres,
    const std::vector<Light> lights,
    size_t depth = 0) {

    Vec3f point, N;
    Material material;

    if (depth > 4 || !SceneIntersect(origin, direction, spheres, point, N, material)) {
        return Vec3f(0.2, 0.7, 0.8);
    }

    Vec3f reflectDirection = Reflect(direction, N).normalize();
    Vec3f reflectOrigin = reflectDirection * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
    Vec3f reflectColor = CastRay(reflectOrigin, reflectDirection, spheres, lights, depth + 1);

    Vec3f refractDirection = Refract(direction, N, material.refractiveIndex);
    Vec3f refractOrigin = refractDirection * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
    Vec3f refractColor = CastRay(refractOrigin, refractDirection, spheres, lights, depth + 1);

    float diffuseLightIntensity = 0;
    float specularLightIntensity = 0;
    float ambientLightIntensity = 0;
    for (size_t i = 0; i < lights.size(); i++) {
        if (lights[i].lightType == Ambient) {
            ambientLightIntensity += lights[i].intensity;
            continue;
        }
        Vec3f lightDirection = (lights[i].position - point).normalize();
        float reflect = Reflect(lightDirection, N) * direction;
        float lightDistance = (lights[i].position - point).norm();

        Vec3f shadowOrigin = lightDirection * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
        Vec3f shadowPoint, shadowN;

        Material tempMaterial;
        if (SceneIntersect(shadowOrigin, lightDirection, spheres, shadowPoint, shadowN, tempMaterial) &&
            (shadowPoint - shadowOrigin).norm() < lightDistance) {
            continue;
        }

        diffuseLightIntensity += std::max(0.0f, lightDirection * N) * lights[i].intensity;
        specularLightIntensity += powf(std::max(0.0f, reflect), material.specularExponent) * lights[i].intensity;
    }

    return material.ambientColor * ambientLightIntensity +
           material.diffuseColor * diffuseLightIntensity * material.albedo[0] +
           Vec3f(1.0, 1.0, 1.0) * specularLightIntensity * material.albedo[1] +
           reflectColor * material.albedo[2] +
           refractColor * material.albedo[3];
}

void Render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
    const int width = 1024;
    const int height = 728;
    const float fov = M_PI / 2.0;
    std::vector<Vec3f> framebuffer(width * height);

    for (size_t j = 0; j < height; j++) {
        for (size_t i = 0; i < width; i++) {
            float x = (2 * (i + 0.5) / (float)width - 1) * tan(fov / 2.) * width / (float)height;
            float y = -(2 * (j + 0.5) / (float)height - 1) * tan(fov / 2.);
            Vec3f direction = Vec3f(x, y, -1).normalize();

            framebuffer[i + j * width] = CastRay(Vec3f(0, 0, 0), direction, spheres, lights, 2);
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
    Vec3f ambientColor = Vec3f(0.2, 0.2, 0.2);
    Material ivory(Vec4f(0.6, 0.3, 0.1, 0.0), Vec3f(0.4, 0.4, 0.3), ambientColor, 50.0, 1.0);
    Material glass(Vec4f(0.0, 0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8), ambientColor, 125.0, 1.5);
    Material redRubber(Vec4f(0.9, 0.1, 0.0, 0.0), Vec3f(0.3, 0.1, 0.1), ambientColor, 10.0, 1.0);
    Material mirror(Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0), ambientColor, 1425.0, 1.0);

    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(Vec3f(-3, 0, -16), 2, ivory));
    spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2, glass));
    spheres.push_back(Sphere(Vec3f(1.5, -0.5, -18), 3, redRubber));
    spheres.push_back(Sphere(Vec3f(7, 5, -18), 4, mirror));
    spheres.push_back(Sphere(Vec3f(-5, -9000, -30), 8995, mirror));

    std::vector<Light> lights;
    lights.push_back(Light(Point, Vec3f(-20, 20, 20), 1.5));
    lights.push_back(Light(Point, Vec3f(30, 50, -25), 1.8));
    lights.push_back(Light(Point, Vec3f(30, 20, 30), 1.7));
    lights.push_back(Light(Ambient, 0.1));

    Render(spheres, lights);

    return 0;
}