
#include "geometry.h"
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <vector>

#define M_PI 3.14159265358979323846

const int MAX_DEPTH = 4;
const double EPSILON = 1e-3;
const Vec3f BACKGROUND_COLOR = Vec3f(0.2, 0.7, 0.8);

struct Material;
class Shape;
class Box;
class Sphere;
class InfinityPlane;
class Light;
class PointLight;
class DirectionalLight;
class SphereLight;
class AmbientClass;

Vec3f Reflect(const Vec3f &direction, const Vec3f &N);
Vec3f Refract(const Vec3f &direction, const Vec3f &N, const double &refractiveIndex);
Vec3f AdjustRayOrigin(Vec3f direction, Vec3f point, Vec3f N);
void Render(const std::vector<std::shared_ptr<Shape>> &shapes, const std::vector<std::shared_ptr<Light>> &lights);

bool SceneIntersect(
    const Vec3f &origin,
    const Vec3f &direction,
    const std::vector<std::shared_ptr<Shape>> &shapes,
    Vec3f &hit,
    Vec3f &N,
    Material &material);

Vec3f CastRay(
    const Vec3f &origin,
    const Vec3f &direction,
    const std::vector<std::shared_ptr<Shape>> &shapes,
    const std::vector<std::shared_ptr<Light>> &lights,
    size_t depth);

Vec3f CalculateFinalColor(
    const Material &material,
    const double &ambientLightIntensity,
    const double &diffuseLightIntensity,
    const double &specularLightIntensity,
    const Vec3f &reflectColor,
    const Vec3f &refractColor);

bool IsInShadow(
    const Vec3f &N,
    const Vec3f &point,
    const Vec3f &lightDirection,
    const double &lightDistance,
    const std::vector<std::shared_ptr<Shape>> &shapes);

struct Material {
    Material(
        const Vec4f &a,
        const Vec3f &color,
        const Vec3f &ambient,
        const double &spec,
        const double &rf) : albedo(a),
                            diffuseColor(color),
                            specularExponent(spec),
                            refractiveIndex(rf),
                            ambientColor(ambient) {}

    Material() : albedo(1, 0, 0, 0), diffuseColor(), specularExponent(), refractiveIndex(1.0), ambientColor(0.1, 0.1, 0.1) {}

    Vec4f albedo;
    Vec3f diffuseColor;
    Vec3f ambientColor;
    double specularExponent;
    double refractiveIndex;
};

class Shape {
protected:
    Material material;

public:
    virtual ~Shape() = default;

    Shape(const Material &material) : material(material) {}

    const Material &GetMaterial() const {
        return material;
    }

    virtual const Vec3f GetNormal(const Vec3f &hitPoint) const = 0;

    virtual bool RayIntersect(
        const Vec3f &origin,
        const Vec3f &direction,
        double &t0) const = 0;
};

class Sphere : public Shape {
private:
    Vec3f center;
    double radius;

    Material material;

public:
    Sphere(
        const Vec3f &center,
        const double &radius,
        const Material &material) : center(center), radius(radius), Shape(material) {}

    const Vec3f &GetCenter() const {
        return center;
    }

    const double &GetRadius() const {
        return radius;
    }

    const Vec3f GetNormal(const Vec3f &hitPoint) const override {
        return (hitPoint - center).normalize();
    }

    bool RayIntersect(
        const Vec3f &origin,
        const Vec3f &direction,
        double &t0) const override {

        Vec3f L = center - origin;
        double tca = L * direction;
        double d2 = L * L - tca * tca;

        if (d2 > radius * radius) {
            return false;
        }

        double thc = sqrtf(radius * radius - d2);
        t0 = tca - thc;
        double t1 = tca + thc;

        if (t0 < 0) {
            t0 = t1;
        }

        if (t0 < 0) {
            return false;
        }

        return true;
    }
};

class Box : public Shape {
private:
    Vec3f boxMax;
    Vec3f boxMin;
    Material material;

public:
    Box(
        Vec3f boxMax,
        Vec3f boxMin,
        Material material) : boxMax(boxMax), boxMin(boxMin), Shape(material) {}

    const Vec3f GetNormal(const Vec3f &hitPoint) const override {

        Vec3f normal;
        if (std::abs(hitPoint.x - boxMin.x) < EPSILON)
            normal = Vec3f(-1, 0, 0); // левая грань
        else if (std::abs(hitPoint.x - boxMax.x) < EPSILON)
            normal = Vec3f(1, 0, 0); // правая грань
        else if (std::abs(hitPoint.y - boxMin.y) < EPSILON)
            normal = Vec3f(0, -1, 0); // нижняя грань
        else if (std::abs(hitPoint.y - boxMax.y) < EPSILON)
            normal = Vec3f(0, 1, 0); // верхняя грань
        else if (std::abs(hitPoint.x - boxMin.z) < EPSILON)
            normal = Vec3f(0, 0, -1); // задняя грань
        else if (std::abs(hitPoint.z - boxMax.z) < EPSILON)
            normal = Vec3f(0, 0, 1); // передняя грань

        return normal;
    }

    bool RayIntersect(
        const Vec3f &origin,
        const Vec3f &direction,
        double &t0) const override {

        Vec3f inversDirection = Vec3f(1.0 / direction.x, 1.0 / direction.y, 1.0 / direction.z);

        double t1 = (boxMin.x - origin.x) * inversDirection.x;
        double t2 = (boxMax.x - origin.x) * inversDirection.x;
        double t3 = (boxMin.y - origin.y) * inversDirection.y;
        double t4 = (boxMax.y - origin.y) * inversDirection.y;
        double t5 = (boxMin.z - origin.z) * inversDirection.z;
        double t6 = (boxMax.z - origin.z) * inversDirection.z;

        double tmin = std::max(std::max(std::min(t1, t2), std::min(t3, t4)), std::min(t5, t6));
        double tmax = std::min(std::min(std::max(t1, t2), std::max(t3, t4)), std::max(t5, t6));

        if (tmax < 0) {
            return false;
        }

        if (tmin > tmax) {
            return false;
        }

        t0 = tmin;
        if (tmin < 0.0) {
            t0 = tmax;
        }

        return true;
    }
};

class InfinityPlane : public Shape {
public:
    InfinityPlane(
        const Vec3f &position,
        const Vec3f &N,
        const Material &material) : position(position), N((position - N).normalize()), Shape(material) {}

    const Vec3f GetNormal(const Vec3f & /*hitPoint*/) const override {
        return this->N;
    }

    bool RayIntersect(
        const Vec3f &origin,
        const Vec3f &direction,
        double &t0) const override {

        double rayPoint = direction * N;
        if (rayPoint < EPSILON) {
            return false;
        }

        double s = (this->N * (this->position - origin)) / rayPoint;

        t0 = (origin + direction * s).norm();

        return true;
    }

private:
    Vec3f position;
    Vec3f N;
};

class Light {
protected:
    double intensity;

public:
    Light(double intensity) : intensity(intensity) {}
    virtual ~Light() = default;

    virtual Vec3f getDirection(const Vec3f &point) const = 0;
    virtual double getDistance(const Vec3f &point) const = 0;
    virtual bool isAmbient() const { return false; }
    double getIntensity() const { return intensity; }
};

class AmbientLight : public Light {
public:
    AmbientLight(double intensity) : Light{intensity} {}

    Vec3f getDirection(const Vec3f & /*point*/) const override {
        return Vec3f(0, 0, 0);
    }

    double getDistance(const Vec3f & /*point*/) const override {
        return 0.0f;
    }

    bool isAmbient() const override { return true; }
};

class PointLight : public Light {
private:
    Vec3f position;

public:
    PointLight(double intensity, const Vec3f &position) : Light{intensity}, position(position) {}

    Vec3f getDirection(const Vec3f &point) const override {
        return (position - point).normalize();
    }

    double getDistance(const Vec3f &point) const override {
        return (position - point).norm();
    }
};

class DirectionalLight : public Light {
private:
    Vec3f direction;

public:
    DirectionalLight(double intensity, const Vec3f &direction) : Light{intensity}, direction(direction) {}

    Vec3f getDirection(const Vec3f &point) const override {
        return direction;
    }

    double getDistance(const Vec3f &point) const override {
        return std::numeric_limits<double>::infinity();
    }
};

// class SphereLight : public Light {
// public:
//     SphereLight(double intensity, const Sphere &sphere) : Light{intensity}, sphere(sphere) {}

//     Vec3f getDirection(const Vec3f &point) const override {
//         Vec3f randomPoint = sampleRandomPointOnSphere();
//         return (randomPoint - point).normalize();
//     }

//     double getDistance(const Vec3f &point) const override {
//         Vec3f randomPoint = sampleRandomPointOnSphere();
//         return (randomPoint - point).norm();
//     }

// private:
//     Vec3f sampleRandomPointOnSphere() const {
//         double u = static_cast<double>(rand()) / RAND_MAX;
//         double v = static_cast<double>(rand()) / RAND_MAX;

//         double theta = 2 * M_PI * u;
//         double phi = acos(2 * v - 1);

//         double x = sphere.GetCenter().x + sphere.GetRadius() * sin(phi) * cos(theta);
//         double y = sphere.GetCenter().y + sphere.GetRadius() * sin(phi) * sin(theta);
//         double z = sphere.GetCenter().z + sphere.GetRadius() * cos(phi);

//         return Vec3f(x, y, z);
//     }
//     const Sphere &sphere;
// };

Vec3f Reflect(const Vec3f &direction, const Vec3f &N) {
    return direction - N * 2.0 * (direction * N);
}

Vec3f Refract(const Vec3f &direction, const Vec3f &N, const double &refractiveIndex) {
    double cosi = -std::max(-1.0, std::min(1.0, direction * N));
    double etai = 1, etat = refractiveIndex;
    Vec3f n = N;
    if (cosi < 0) {
        cosi *= -1;
        std::swap(etai, etat);
        n = -N;
    }
    double eta = etai / etat;
    double k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? Vec3f(0, 0, 0) : direction * eta + n * (eta * cosi - sqrtf(k));
}

bool SceneIntersect(
    const Vec3f &origin,
    const Vec3f &direction,
    const std::vector<std::shared_ptr<Shape>> &shapes,
    Vec3f &hit,
    Vec3f &N,
    Material &material) {

    double spheresDist = std::numeric_limits<double>::max();

    for (size_t i = 0; i < shapes.size(); i++) {
        double dist_i;
        if (shapes[i]->RayIntersect(origin, direction, dist_i) && dist_i < spheresDist) {
            spheresDist = dist_i;
            hit = origin + direction * dist_i;
            N = shapes[i]->GetNormal(hit);
            material = shapes[i]->GetMaterial();
        }
    }

    return spheresDist < 1000;
}

Vec3f CastRay(
    const Vec3f &origin,
    const Vec3f &direction,
    const std::vector<std::shared_ptr<Shape>> &shapes,
    const std::vector<std::shared_ptr<Light>> &lights,
    size_t depth = 0) {

    Vec3f point, N;
    Material material;

    if (depth > MAX_DEPTH || !SceneIntersect(origin, direction, shapes, point, N, material)) {
        return BACKGROUND_COLOR;
    }

    Vec3f reflectDirection = Reflect(direction, N).normalize();
    Vec3f reflectOrigin = AdjustRayOrigin(reflectDirection, point, N);
    Vec3f reflectColor = CastRay(
        reflectOrigin,
        reflectDirection,
        shapes, lights, depth + 1);

    Vec3f refractDirection = Refract(direction, N, material.refractiveIndex);
    Vec3f refractOrigin = AdjustRayOrigin(refractDirection, point, N);
    Vec3f refractColor = CastRay(
        refractOrigin,
        refractDirection,
        shapes, lights, depth + 1);

    double diffuseLightIntensity = 0;
    double specularLightIntensity = 0;
    double ambientLightIntensity = 0;
    for (const auto &light : lights) {
        if (light->isAmbient()) {
            ambientLightIntensity += light->getIntensity();
            continue;
        }

        Vec3f lightDirection = light->getDirection(point);
        double lightDistance = light->getDistance(point);
        double reflect = Reflect(lightDirection, N) * direction;

        if (IsInShadow(N, point, lightDirection, lightDistance, shapes)) {
            continue;
        }

        diffuseLightIntensity += std::max(0.0, lightDirection * N) * light->getIntensity();
        specularLightIntensity += powf(std::max(0.0, reflect), material.specularExponent) * light->getIntensity();
    }

    return CalculateFinalColor(
        material, ambientLightIntensity,
        diffuseLightIntensity, specularLightIntensity,
        reflectColor, refractColor);
}

bool IsInShadow(
    const Vec3f &N,
    const Vec3f &point,
    const Vec3f &lightDirection,
    const double &lightDistance,
    const std::vector<std::shared_ptr<Shape>> &shapes) {

    Vec3f shadowOrigin = AdjustRayOrigin(lightDirection, point, N);
    Vec3f shadowPoint, shadowN;

    Material tempMaterial;

    return SceneIntersect(shadowOrigin, lightDirection, shapes, shadowPoint, shadowN, tempMaterial) &&
           (shadowPoint - shadowOrigin).norm() < lightDistance;
}

Vec3f AdjustRayOrigin(Vec3f direction, Vec3f point, Vec3f N) {
    return direction * N < 0 ? point - N * EPSILON : point + N * EPSILON;
}

Vec3f CalculateFinalColor(
    const Material &material,
    const double &ambientLightIntensity,
    const double &diffuseLightIntensity,
    const double &specularLightIntensity,
    const Vec3f &reflectColor,
    const Vec3f &refractColor) {

    return material.ambientColor * ambientLightIntensity +
           material.diffuseColor * diffuseLightIntensity * material.albedo[0] +
           Vec3f(1.0, 1.0, 1.0) * specularLightIntensity * material.albedo[1] +
           reflectColor * material.albedo[2] +
           refractColor * material.albedo[3];
}

void Render(const std::vector<std::shared_ptr<Shape>> &shapes, const std::vector<std::shared_ptr<Light>> &lights) {
    const int width = 1024;
    const int height = 720;
    const double fov = M_PI / 2.0;
    std::vector<Vec3f> framebuffer(width * height);

    for (size_t j = 0; j < height; j++) {
        for (size_t i = 0; i < width; i++) {
            double x = (2 * (i + 0.5) / (double)width - 1) * tan(fov / 2.) * width / (double)height;
            double y = -(2 * (j + 0.5) / (double)height - 1) * tan(fov / 2.);
            Vec3f direction = Vec3f(x, y, -1).normalize();

            framebuffer[i + j * width] = CastRay(Vec3f(0, 0, 0), direction, shapes, lights);
        }
    }

    std::ofstream ofs;
    ofs.open("./out.ppm", std::ios::binary);
    ofs << "P6\n"
        << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height * width; ++i) {
        for (size_t j = 0; j < 3; j++) {
            ofs << static_cast<char>(255 * std::max(0.0, std::min(1.0, framebuffer[i][j])));
        }
    }
    ofs.close();
}

int main() {
    Vec3f ambientColor = Vec3f(0.2, 0.2, 0.2);
    Material ivory(Vec4f(0.6, 0.3, 0.3, 0.0), Vec3f(0.4, 0.4, 0.3), ambientColor, 1000.0, 1.0);
    Material glass(Vec4f(0.0, 0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8), ambientColor, 125.0, 1.5);
    Material redRubber(Vec4f(0.9, 0.1, 0.0, 0.0), Vec3f(0.3, 0.1, 0.1), ambientColor, 10.0, 1.0);
    Material mirror(Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0), ambientColor, 1425.0, 1.0);

    // std::vector<Sphere> spheres;
    // spheres.push_back(Sphere(Vec3f(-3, 0, -16), 2, ivory));
    // spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2, glass));
    // spheres.push_back(Sphere(Vec3f(1.5, -0.5, -18), 3, redRubber));
    // spheres.push_back(Sphere(Vec3f(7, 5, -18), 4, mirror));
    // spheres.push_back(Sphere(Vec3f(-5, -900, -30), 895, ivory));

    std::vector<std::shared_ptr<Shape>> shapes;
    // shapes.push_back(std::make_shared<Sphere>(Vec3f(-3, 0, -16), 2, ivory));
    // shapes.push_back(std::make_shared<Sphere>(Vec3f(-1.0, -1.5, -12), 2, glass));
    // shapes.push_back(std::make_shared<Sphere>(Vec3f(1.5, -0.5, -18), 3, redRubber));
    shapes.push_back(std::make_shared<Sphere>(Vec3f(7, 5, -18), 4, mirror));
    shapes.push_back(std::make_shared<Sphere>(Vec3f(-5, -9000, -30), 8995, mirror));
    
    
    shapes.push_back(std::make_shared<Box>(Vec3f(-3, 2, -5), Vec3f(-1.0, -1.5, -12), ivory));
    shapes.push_back(std::make_shared<Box>(Vec3f(-3, 0, -16), Vec3f(5, 5, -20), redRubber));

    //shapes.push_back(std::make_shared<InfinityPlane>(Vec3f(0, -100, 0), Vec3f(24, -25, -3), ivory));

    std::vector<std::shared_ptr<Light>> lights;
    lights.push_back(std::make_shared<PointLight>(PointLight(1.5, Vec3f(-20, 20, 20))));
    lights.push_back(std::make_shared<PointLight>(PointLight(1.8, Vec3f(30, 50, -25))));
    lights.push_back(std::make_shared<PointLight>(PointLight(1.7, Vec3f(30, 20, 30))));
    // lights.push_back(std::make_shared<SphereLight>(SphereLight(5, Sphere(Vec3f(1.5, -0.5, -18), 3, redRubber))));
    lights.push_back(std::make_shared<AmbientLight>(AmbientLight(0.1)));
    // lights.push_back(std::make_shared<DirectionalLight>(DirectionalLight(0.1, Vec3f(-6, 5, 0.78))));

    Render(shapes, lights);

    return 0;
}
