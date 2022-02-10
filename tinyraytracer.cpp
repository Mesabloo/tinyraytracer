#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "csg.h"
#include "geometry.h"
#include "material.h"
#include "model.h"

#define WITH_DUCK 0

int envmap_width, envmap_height;
std::vector<Vec3f> envmap;

struct Light {
    Light(const Vec3f &p, const float i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

static inline Vec3f reflect(const Vec3f &I, const Vec3f &N) { return I - N * 2.f * (I * N); }

static inline Vec3f refract(const Vec3f &I, const Vec3f &N, const float eta_t,
                            const float eta_i = 1.f) { // Snell's law
    const float cosi = -std::clamp(I * N, -1.f, 1.f);  // std::max(-1.f, std::min(1.f, I * N));
    if (cosi < 0)
        return refract(I, -N, eta_i, eta_t); // if the ray comes from the inside the
                                             // object, swap the air and the media
    const float eta = eta_i / eta_t;
    const float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? Vec3f(1, 0, 0) : I * eta + N * (eta * cosi - sqrtf(k)); // k<0 = total reflection, no ray to
                                                                           // refract. I refract it anyways,
                                                                           // this has no physical meaning
}

bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const std::shared_ptr<CSGTree<Shape>> &csg_tree, Vec3f &hit,
                     Vec3f &N, Material &material) {
    const std::set<Interval> intersections = csg_tree->ray_intersect(orig, dir);

    if (intersections.empty())
        return false;

    const Interval &i = *std::cbegin(intersections);

    std::clog << i << std::endl;

    hit = i.from;
    N = i.compute_normal(hit);
    material = i.material;

    return (orig - hit).norm() < 1000;
}

#define MAX_DEPTH 4

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::shared_ptr<CSGTree<Shape>> &csg_tree,
               const std::vector<Light> &lights, size_t depth = 0) {
    Vec3f point, N;
    Material material;

    if (depth > MAX_DEPTH || !scene_intersect(orig, dir, csg_tree, point, N, material)) {
        int a =
            std::clamp(static_cast<int>((atan2(dir.z, dir.x) / (2 * M_PI) + .5) * envmap_width), 0, envmap_width - 1);
        int b = std::clamp(static_cast<int>(acos(dir.y) / M_PI * envmap_height), 0, envmap_height - 1);
        return envmap[a + b * envmap_width]; // background color
    }

    const Vec3f reflect_dir = reflect(dir, N).normalize();
    const Vec3f refract_dir = refract(dir, N, material.refractive_index).normalize();
    const Vec3f reflect_orig =
        reflect_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3; // offset the original point to avoid occlusion by
                                                                   // the object itself
    const Vec3f refract_orig = refract_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
    const Vec3f reflect_color = cast_ray(reflect_orig, reflect_dir, csg_tree, lights, depth + 1);
    const Vec3f refract_color = cast_ray(refract_orig, refract_dir, csg_tree, lights, depth + 1);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
//#pragma omp parallel for
    for (const Light &light : lights) {
        Vec3f vec = light.position - point;
        const Vec3f light_dir = vec.normalize();
        const float light_distance = vec.norm();

        Vec3f shadow_orig = light_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3; // checking if the point lies in
                                                                                     // the shadow of the light
        Vec3f shadow_pt, shadow_N;
        Material tmpmaterial;
        if (scene_intersect(shadow_orig, light_dir, csg_tree, shadow_pt, shadow_N, tmpmaterial) &&
            (shadow_pt - shadow_orig).norm() < light_distance)
            continue;

        diffuse_light_intensity += light.intensity * std::max(0.f, light_dir * N);
        specular_light_intensity +=
            powf(std::max(0.f, -reflect(-light_dir, N) * dir), material.specular_exponent) * light.intensity;
    }
    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] +
           Vec3f(1., 1., 1.) * specular_light_intensity * material.albedo[1] + reflect_color * material.albedo[2] +
           refract_color * material.albedo[3];
}

static inline void render_normal(const std::shared_ptr<CSGTree<Shape>> &csg_tree, const std::vector<Light> &lights,
                                 const Vec3f &camera, const int i) {
    const int width = 1024;
    const int height = 768;
    const float fov = M_PI / 3.;
    std::vector<Vec3f> framebuffer(width * height);

#pragma omp parallel for
    for (size_t j = 0; j < height; j++) { // actual rendering loop
        for (size_t i = 0; i < width; i++) {
            const float dir_x = (i + 0.5) - width / 2.;
            const float dir_y = -(j + 0.5) + height / 2.; // this flips the image at the same time
            const float dir_z = -height / (2. * tan(fov / 2.));

            framebuffer[i + j * width] = cast_ray(camera, Vec3f{dir_x, dir_y, dir_z}.normalize(), csg_tree, lights);
        }
    }

    std::vector<unsigned char> pixmap(width * height * 3);
#pragma omp parallel for
    for (size_t i = 0; i < height * width; ++i) {
        Vec3f &c = framebuffer[i];
        const float max = std::max({c[0], c[1], c[2]});
        if (max > 1)
            c = c * (1. / max);
        for (size_t j = 0; j < 3; j++) {
            pixmap[i * 3 + j] = (unsigned char)(255 * std::clamp(framebuffer[i][j], 0.f, 1.f));
            // std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    stbi_write_jpg(("out_normal_" + std::to_string(i) + ".jpg").c_str(), width, height, 3, pixmap.data(), 100);
}

static inline void render_parallax(const std::shared_ptr<CSGTree<Shape>> &csg_tree, const std::vector<Light> &lights,
                                   const Vec3f &camera, const int i) {
    const int delta = 60; // focal distance 3
    const int width = 1024 + delta;
    const int height = 768;
    const float eyesep = 0.2;
    const float fov = M_PI / 3.;
    std::vector<Vec3f> framebuffer1(width * height);
    std::vector<Vec3f> framebuffer2(width * height);

#pragma omp parallel for
    for (size_t j = 0; j < height; j++) { // actual rendering loop
        for (size_t i = 0; i < width; i++) {
            float dir_x = (i + 0.5) - width / 2.;
            float dir_y = -(j + 0.5) + height / 2.; // this flips the image at the same time
            float dir_z = -height / (2. * tan(fov / 2.));
            // framebuffer[i+j*width] = cast_ray(Vec3f(0,0,0), Vec3f(dir_x, dir_y,
            // dir_z).normalize(), spheres, lights);
            framebuffer1[i + j * width] = cast_ray(Vec3f(camera.x - eyesep / 2, camera.y, camera.z),
                                                   Vec3f(dir_x, dir_y, dir_z).normalize(), csg_tree, lights);
            framebuffer2[i + j * width] = cast_ray(Vec3f(camera.x + eyesep / 2, camera.y, camera.z),
                                                   Vec3f(dir_x, dir_y, dir_z).normalize(), csg_tree, lights);
        }
    }

    std::vector<unsigned char> pixmap((width - delta) * height * 3);
    for (size_t j = 0; j < height; j++) {
        for (size_t i = 0; i < width - delta; i++) {
            Vec3f c1 = framebuffer1[i + delta + j * width];
            Vec3f c2 = framebuffer2[i + j * width];

            float max1 = std::max({c1[0], c1[1], c1[2]});
            if (max1 > 1)
                c1 = c1 * (1. / max1);
            float max2 = std::max({c2[0], c2[1], c2[2]});
            if (max2 > 1)
                c2 = c2 * (1. / max2);
            float avg1 = (c1.x + c1.y + c1.z) / 3.;
            float avg2 = (c2.x + c2.y + c2.z) / 3.;

            pixmap[(j * (width - delta) + i) * 3] = 255 * avg1;
            pixmap[(j * (width - delta) + i) * 3 + 1] = 0;
            pixmap[(j * (width - delta) + i) * 3 + 2] = 255 * avg2;
        }
    }
    stbi_write_jpg(("out_parallax_" + std::to_string(i) + ".jpg").c_str(), width - delta, height, 3, pixmap.data(),
                   100);
}

static inline void render_stereoscope(const std::shared_ptr<CSGTree<Shape>> &csg_tree, const std::vector<Light> &lights,
                                      const Vec3f &camera, const int i) {
    const int delta = 60; // focal distance 3
    const int width = 1024 + delta;
    const int height = 768;
    const float eyesep = 0.2;
    const float fov = M_PI / 3.;
    std::vector<Vec3f> framebuffer1(width * height);
    std::vector<Vec3f> framebuffer2(width * height);

#pragma omp parallel for
    for (size_t j = 0; j < height; j++) { // actual rendering loop
        for (size_t i = 0; i < width; i++) {
            float dir_x = (i + 0.5) - width / 2.;
            float dir_y = -(j + 0.5) + height / 2.; // this flips the image at the same time
            float dir_z = -height / (2. * tan(fov / 2.));
            // framebuffer[i+j*width] = cast_ray(Vec3f(0,0,0), Vec3f(dir_x, dir_y,
            // dir_z).normalize(), spheres, lights);
            framebuffer1[i + j * width] = cast_ray(Vec3f(camera.x - eyesep / 2, camera.y, camera.z),
                                                   Vec3f(dir_x, dir_y, dir_z).normalize(), csg_tree, lights);
            framebuffer2[i + j * width] = cast_ray(Vec3f(camera.x + eyesep / 2, camera.y, camera.z),
                                                   Vec3f(dir_x, dir_y, dir_z).normalize(), csg_tree, lights);
        }
    }

    const float k1 = 0.12;
    const float k2 = 0.10;

    const float xc = (width - delta) / 2.;
    const float yc = height / 2.;
    const float R = std::min(width - delta, height) / 2.f;

    std::vector<unsigned char> pixmap((width - delta) * height * 3 * 2);
    for (size_t j = 0; j < height; j++) {
        for (size_t i = 0; i < width - delta; i++) {
            float xd = i;
            float yd = j;
            float r = std::sqrt(pow(xd - xc, 2) + pow(yd - yc, 2)) / R;

            int xu = xc + (xd - xc) * (1 + k1 * pow(r, 2) + k2 * pow(r, 4));
            int yu = yc + (yd - yc) * (1 + k1 * pow(r, 2) + k2 * pow(r, 4));

            Vec3f c1(0, 0, 0), c2(0, 0, 0);
            if (xu >= 0 && xu < width - delta && yu >= 0 && yu < height) {
                c1 = framebuffer2[xu + yu * width + delta];
                c2 = framebuffer2[xu + yu * width];
            }

            float max2 = std::max({c2[0], c2[1], c2[2]});
            if (max2 > 1)
                c2 = c2 * (1. / max2);
            float max1 = std::max({c1[0], c1[1], c1[2]});
            if (max1 > 1)
                c1 = c1 * (1. / max1);

            for (size_t d = 0; d < 3; d++) {
                pixmap[(j * (width - delta) * 2 + i) * 3 + d] = 255 * c1[d];
                pixmap[(j * (width - delta) * 2 + i + width - delta) * 3 + d] = 255 * c2[d];
            }
        }
    }
    stbi_write_jpg(("out_stereo_" + std::to_string(i) + ".jpg").c_str(), (width - delta) * 2, height, 3, pixmap.data(),
                   100);
}

#define NB_IMAGES 1

static inline void render_video_rebond(void (*renderer)(const std::shared_ptr<CSGTree<Shape>> &,
                                                        const std::vector<Light> &, const Vec3f &, const int),
                                       std::vector<Sphere> &spheres, const std::shared_ptr<CSGTree<Shape>> &csg_tree,
                                       const std::vector<Light> &lights) {
    int first = 0;
    for (int angle = 1; angle <= NB_IMAGES + 1; ++angle) {

        std::clog << "\033[0G\033[2K[" << angle << "/" << NB_IMAGES << "] Generating video...";

        renderer(csg_tree, lights, Vec3f{0, 0, 0}, angle);

        switch (first) {
        case 0:
            spheres[2].center.x -= 1;
            spheres[2].center.y += 1;
            break;
        case 1:
            spheres[2].center.x -= 1;
            spheres[2].center.y -= 1;
            break;
        case 2:
            spheres[2].center.z += 1;
            spheres[2].center.y += 1;
            break;
        case 3:
            spheres[2].center.z += 1;
            spheres[2].center.y -= 1;
            break;
        case 4:
            spheres[2].center.x += 1;
            spheres[2].center.y += 1;
            break;
        case 5:
            spheres[2].center.x += 1;
            spheres[2].center.y -= 1;
            break;
        case 6:
            spheres[2].center.z -= 1;
            spheres[2].center.y += 1;
            break;
        case 7:
            spheres[2].center.z -= 1;
            spheres[2].center.y -= 1;
            break;
        }

        // changement tous les 5 images
        if (angle % 5 == 0) {
            first = (first + 1) % 8;
        }
    }
}

int main() {
    int n = -1;
    unsigned char *pixmap = stbi_load("../envmap.jpg", &envmap_width, &envmap_height, &n, 0);
    if (!pixmap || 3 != n) {
        std::cerr << "Error: can not load the environment map" << std::endl;
        return -1;
    }
    envmap = std::vector<Vec3f>(envmap_width * envmap_height);
#pragma omp parallel for
    for (int j = envmap_height - 1; j >= 0; j--) {
        for (int i = 0; i < envmap_width; i++) {
            const int dup_comp = i + j * envmap_width;
            const int offset = dup_comp * 3;
            envmap[dup_comp] = Vec3f(pixmap[offset + 0], pixmap[offset + 1], pixmap[offset + 2]) * (1 / 255.);
        }
    }
    stbi_image_free(pixmap);

    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(Vec3f(-3, 0, -16), 2, ivory));
    spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2, glass));
    spheres.push_back(Sphere(Vec3f(1.5, -0.5, -18), 3, red_rubber));
    spheres.push_back(Sphere(Vec3f(7, 5, -18), 4, mirror));
    spheres.push_back(Sphere(Vec3f(-10, 7, -22), 3, ivory));     // premier cercle
    spheres.push_back(Sphere(Vec3f(-8, 8, -20), 1, red_rubber)); // deuxieme cercle

    Checkerboard checkerboard(4);

    DuckObject duck("../duck.obj");

    const std::shared_ptr<CSGTree<Shape>> death_star = std::make_shared<CSGTreeNode<Shape>>( // death star
        std::make_shared<CSGTreeLeaf<Shape>>(spheres[4]), std::make_shared<CSGTreeLeaf<Shape>>(spheres[5]),
        CSGTreeNodeType::DIFFERENCE);

    const std::shared_ptr<CSGTree<Shape>> csg_tree = std::make_shared<CSGTreeNode<Shape>>(
        std::make_shared<CSGTreeNode<Shape>>(
            death_star,
            std::make_shared<CSGTreeNode<Shape>>(std::make_shared<CSGTreeLeaf<Shape>>(spheres[0]),
                                                 std::make_shared<CSGTreeLeaf<Shape>>(spheres[1]),
                                                 CSGTreeNodeType::UNION),
            CSGTreeNodeType::UNION),
        std::make_shared<CSGTreeNode<Shape>>(
            std::make_shared<CSGTreeNode<Shape>>(std::make_shared<CSGTreeLeaf<Shape>>(spheres[2]),
                                                 std::make_shared<CSGTreeLeaf<Shape>>(spheres[3]),
                                                 CSGTreeNodeType::UNION),
#if WITH_DUCK
            std::make_shared<CSGTreeNode<Shape>>(std::make_shared<CSGTreeLeaf<Shape>>(duck),
                                                 std::make_shared<CSGTreeLeaf<Shape>>(checkerboard),
                                                 CSGTreeNodeType::UNION),
#else
            std::make_shared<CSGTreeLeaf<Shape>>(checkerboard),
#endif
            CSGTreeNodeType::UNION),

        CSGTreeNodeType::UNION);

    const std::shared_ptr<CSGTree<Shape>> csg_tree2 = std::make_shared<CSGTreeNode<Shape>>(
        std::make_shared<CSGTreeLeaf<Shape>>(duck), std::make_shared<CSGTreeLeaf<Shape>>(checkerboard),
        CSGTreeNodeType::UNION);

    std::vector<Light> lights;
    lights.push_back(Light(Vec3f(-20, 20, 20), 1.5));
    lights.push_back(Light(Vec3f(30, 50, -25), 1.8));
    lights.push_back(Light(Vec3f(30, 20, 30), 1.7));

#if 0
    std::clog << "Rendering normal image..." << std::endl;
    render_normal(csg_tree, lights, 0);
    std::clog << "Rendering parallax image..." << std::endl;
    render_parallax(csg_tree, lights, 0);
    std::clog << "Rendering stereoscope image..." << std::endl;
    render_stereoscope(csg_tree, lights, 0);
#endif
    render_video_rebond(&render_normal, spheres, csg_tree2, lights);

    return 0;
}
