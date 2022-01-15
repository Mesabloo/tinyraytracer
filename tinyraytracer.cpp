#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "geometry.h"
#include "model.h"

#define WITH_DUCK 0

int envmap_width, envmap_height;
std::vector<Vec3f> envmap;
//      _        _        _
//   __(.)<   __(.)>   __(.)=
//   \___)    \___)    \___)
Model duck("../duck.obj");

struct Light {
  Light(const Vec3f &p, const float i) : position(p), intensity(i) {}
  Vec3f position;
  float intensity;
};

struct Material {
  Material(const float r, const Vec4f &a, const Vec3f &color, const float spec)
      : refractive_index(r), albedo(a), diffuse_color(color),
        specular_exponent(spec) {}
  Material()
      : refractive_index(1), albedo(1, 0, 0, 0), diffuse_color(),
        specular_exponent() {}
  float refractive_index;
  Vec4f albedo;
  Vec3f diffuse_color;
  float specular_exponent;
};

const Material ivory(1.0, Vec4f(0.6, 0.3, 0.1, 0.0), Vec3f(0.4, 0.4, 0.3), 50.);
const Material glass(1.5, Vec4f(0.0, 0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8),
                     125.);
const Material red_rubber(1.0, Vec4f(0.9, 0.1, 0.0, 0.0), Vec3f(0.3, 0.1, 0.1),
                          10.);
const Material mirror(1.0, Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0),
                      1425.);

struct Sphere {
  Vec3f center;
  float radius;
  Material material;

  Sphere(const Vec3f &c, const float r, const Material &m)
      : center(c), radius(r), material(m) {}

  inline bool ray_intersect(const Vec3f &orig, const Vec3f &dir,
                            float &t0) const {
    const Vec3f L = center - orig;
    const float tca = L * dir;
    const float d2 = L * L - tca * tca;
    if (d2 > radius * radius)
      return false;
    const float thc = sqrtf(radius * radius - d2);
    t0 = tca - thc;
    const float t1 = tca + thc;
    if (t0 < 0)
      t0 = t1;
    return t0 >= 0;
  }
};

static inline Vec3f rotate_camera(const Vec3f &orig, const Vec3f &dir, const Vec3f &target) {
    const Vec3f zAxis = (orig - target).normalize();
    const Vec3f xAxis = cross(Vec3f(0., 1., 0.), zAxis).normalize();
    const Vec3f yAxis = cross(zAxis, xAxis).normalize();
    const Mat4f transform = Mat4f(Vec4f(xAxis.x, xAxis.y, xAxis.z, 0.), Vec4f(yAxis.x, yAxis.y, yAxis.z, 0.), Vec4f(zAxis.x, zAxis.y, zAxis.z, 0.), Vec4f(orig.x, orig.y, orig.z, 1.));
    const Vec4f transformed = transform * Vec4f(dir.x, dir.y, dir.z, 0.);
    return Vec3f(transformed[0], transformed[1], transformed[2]);
}

static inline Vec3f reflect(const Vec3f &I, const Vec3f &N) {
  return I - N * 2.f * (I * N);
}

static inline Vec3f refract(const Vec3f &I, const Vec3f &N, const float eta_t,
                            const float eta_i = 1.f) { // Snell's law
  const float cosi =
      -std::clamp(I * N, -1.f, 1.f); // std::max(-1.f, std::min(1.f, I * N));
  if (cosi < 0)
    return refract(I, -N, eta_i, eta_t); // if the ray comes from the inside the
                                         // object, swap the air and the media
  const float eta = eta_i / eta_t;
  const float k = 1 - eta * eta * (1 - cosi * cosi);
  return k < 0 ? Vec3f(1, 0, 0)
               : I * eta + N * (eta * cosi -
                                sqrtf(k)); // k<0 = total reflection, no ray to
                                           // refract. I refract it anyways,
                                           // this has no physical meaning
}

bool scene_intersect(const Vec3f &orig, const Vec3f &dir,
                     const std::vector<Sphere> &spheres, Vec3f &hit, Vec3f &N,
                     Material &material) {
  float spheres_dist = std::numeric_limits<float>::max();
  //#pragma omp parallel for
  for (const Sphere &sphere : spheres) {
    float dist_i;
    if (sphere.ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist) {
      spheres_dist = dist_i;
      hit = orig + dir * dist_i;
      N = (hit - sphere.center).normalize();
      material = sphere.material;
    }
  }

  float checkerboard_dist = std::numeric_limits<float>::max();
  if (fabs(dir.y) > 1e-3) {
    const float d =
        -(orig.y + 4) / dir.y; // the checkerboard plane has equation y = -4
    const Vec3f pt = orig + dir * d;
    if (d > 0 && fabs(pt.x) < 10 && pt.z < -10 && pt.z > -30 &&
        d < spheres_dist) {
      checkerboard_dist = d;
      hit = pt;
      N = Vec3f(0, 1, 0);
      material.diffuse_color = (int(.5 * hit.x + 1000) + int(.5 * hit.z)) & 1
                                   ? Vec3f(.3, .3, .3)
                                   : Vec3f(.3, .2, .1);
    }
  }

#if WITH_DUCK
  //      _        _        _
  //   __(.)<   __(.)>   __(.)=
  //   \___)    \___)    \___)
  float duck_dist = std::numeric_limits<float>::max();
  const int nfaces = duck.nfaces();
  //#pragma omp parallel for
  for (int i = 0; i < nfaces; ++i) {
    float duck_dist_i;
    // On cherche à afficher la face du canard qui est la plus proche,
    // uniquement si le rayon la traverse
    if (duck.ray_triangle_intersect(i, orig, dir, duck_dist_i) &&
        duck_dist_i < duck_dist) {
      duck_dist = duck_dist_i;
      hit = orig + dir * duck_dist_i;

      // Les 3 points du triangle
      const Vec3f a = duck.point(duck.vert(i, 0));
      const Vec3f b = duck.point(duck.vert(i, 1));
      const Vec3f c = duck.point(duck.vert(i, 2));

      // Le vecteur normal à un triangle est un vecteur orthogonal à 2 de ses
      // côtés
      //
      // ```
      //                                    /
      //                             (a)   /
      //                                  /
      //                          _,.-'"\/
      //                    _,.-'"      /\
      //              _,.-'"           /  \
      //        _,.-'"                /    \
      // (c)   ('-._                 /      \
      //            '-._            .        \
      //                '-._                  \
      //                    '-._               \
      //                        '-._            \
      //                        /   '-._         \
      //                       /        '-:_      \
      //                      /             '-._   \
      //                     /                  '-._)
      //                    /                          (b)
      //                   /
      //
      //                 (N)
      //
      // ```
      //
      // Ici,
      // * `b - a` représente le vecteur `AB`
      // * `c - a` représente le vecteur `AC`
      //
      // Le produit vectoriel des vecteurs `AB` et `AC` correspond au vecteur
      // normal `N`
      N = cross(b - a, c - a).normalize();

      material = glass;
    }
  }

  return std::min({duck_dist, spheres_dist, checkerboard_dist}) < 1000;
#else
  return std::min({spheres_dist, checkerboard_dist}) < 1000;
#endif
}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir,
               const std::vector<Sphere> &spheres,
               const std::vector<Light> &lights, size_t depth = 0) {
  Vec3f point, N;
  Material material;

  if (depth > 4 || !scene_intersect(orig, dir, spheres, point, N, material)) {
    int a =
        std::clamp(static_cast<int>((atan2(dir.z, dir.x) / (2 * M_PI) + .5) *
                                    envmap_width),
                   0, envmap_width - 1);
    int b = std::clamp(static_cast<int>(acos(dir.y) / M_PI * envmap_height), 0,
                       envmap_height - 1);
    return envmap[a + b * envmap_width]; // background color
  }

  const Vec3f reflect_dir = reflect(dir, N).normalize();
  const Vec3f refract_dir =
      refract(dir, N, material.refractive_index).normalize();
  const Vec3f reflect_orig =
      reflect_dir * N < 0
          ? point - N * 1e-3
          : point + N * 1e-3; // offset the original point to avoid occlusion by
                              // the object itself
  const Vec3f refract_orig =
      refract_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
  const Vec3f reflect_color =
      cast_ray(reflect_orig, reflect_dir, spheres, lights, depth + 1);
  const Vec3f refract_color =
      cast_ray(refract_orig, refract_dir, spheres, lights, depth + 1);

  float diffuse_light_intensity = 0, specular_light_intensity = 0;
#pragma omp parallel for
  for (const Light &light : lights) {
    Vec3f vec = light.position - point;
    const Vec3f light_dir = vec.normalize();
    const float light_distance = vec.norm();

    Vec3f shadow_orig = light_dir * N < 0
                            ? point - N * 1e-3
                            : point + N * 1e-3; // checking if the point lies in
                                                // the shadow of the light
    Vec3f shadow_pt, shadow_N;
    Material tmpmaterial;
    if (scene_intersect(shadow_orig, light_dir, spheres, shadow_pt, shadow_N,
                        tmpmaterial) &&
        (shadow_pt - shadow_orig).norm() < light_distance)
      continue;

    diffuse_light_intensity += light.intensity * std::max(0.f, light_dir * N);
    specular_light_intensity +=
        powf(std::max(0.f, -reflect(-light_dir, N) * dir),
             material.specular_exponent) *
        light.intensity;
  }
  return material.diffuse_color * diffuse_light_intensity * material.albedo[0] +
         Vec3f(1., 1., 1.) * specular_light_intensity * material.albedo[1] +
         reflect_color * material.albedo[2] +
         refract_color * material.albedo[3];
}

static inline void render_normal(const std::vector<Sphere> &spheres,
                                 const std::vector<Light> &lights,
                                 const Vec3f &camera, const Vec3f &dir,
                                 const int i) {
  const int width = 1024;
  const int height = 768;
  const float fov = M_PI / 3.;
  std::vector<Vec3f> framebuffer(width * height);

#pragma omp parallel for
  for (size_t j = 0; j < height; j++) { // actual rendering loop
    for (size_t i = 0; i < width; i++) {
      const float dir_x = (i + 0.5) - width / 2.;
      const float dir_y =
          -(j + 0.5) + height / 2.; // this flips the image at the same time
      const float dir_z = -height / (2. * tan(fov / 2.));

      //Vec3f dir_{dir_x, dir_y, dir_z};
      //dir_ = rotate_camera(camera, dir_.normalize(), Vec3f{0, 0, 0}).normalize();
      // const float theta = 45.f * M_PI / 180.f;
      // const float new_dir_x = dir_.x * cos(theta) - dir_.z * sin(theta);
      // const float new_dir_z = -dir_.x * sin(theta) + dir_.z * cos(theta);
      // dir_.x = new_dir_x;
      // dir_.z = new_dir_z;
      //dir_.normalize();

      framebuffer[i + j * width] = cast_ray(camera, Vec3f{dir_x, dir_y, dir_z}.normalize(), spheres, lights);
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
      pixmap[i * 3 + j] =
          (unsigned char)(255 * std::clamp(framebuffer[i][j], 0.f, 1.f));
      // std::max(0.f, std::min(1.f, framebuffer[i][j])));
    }
  }
  stbi_write_jpg(("out_normal_" + std::to_string(i) + ".jpg").c_str(), width,
                 height, 3, pixmap.data(), 100);
}

static inline void render_parallax(const std::vector<Sphere> &spheres,
                                   const std::vector<Light> &lights,
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
      float dir_y =
          -(j + 0.5) + height / 2.; // this flips the image at the same time
      float dir_z = -height / (2. * tan(fov / 2.));
      // framebuffer[i+j*width] = cast_ray(Vec3f(0,0,0), Vec3f(dir_x, dir_y,
      // dir_z).normalize(), spheres, lights);
      framebuffer1[i + j * width] =
          cast_ray(Vec3f(camera.x - eyesep / 2, camera.y, camera.z),
                   Vec3f(dir_x, dir_y, dir_z).normalize(), spheres, lights);
      framebuffer2[i + j * width] =
          cast_ray(Vec3f(camera.x + eyesep / 2, camera.y, camera.z),
                   Vec3f(dir_x, dir_y, dir_z).normalize(), spheres, lights);
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
  stbi_write_jpg(("out_parallax_" + std::to_string(i) + ".jpg").c_str(),
                 width - delta, height, 3, pixmap.data(), 100);
}

static inline void render_stereoscope(const std::vector<Sphere> &spheres,
                                      const std::vector<Light> &lights,
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
      float dir_y =
          -(j + 0.5) + height / 2.; // this flips the image at the same time
      float dir_z = -height / (2. * tan(fov / 2.));
      // framebuffer[i+j*width] = cast_ray(Vec3f(0,0,0), Vec3f(dir_x, dir_y,
      // dir_z).normalize(), spheres, lights);
      framebuffer1[i + j * width] =
          cast_ray(Vec3f(camera.x - eyesep / 2, camera.y, camera.z),
                   Vec3f(dir_x, dir_y, dir_z).normalize(), spheres, lights);
      framebuffer2[i + j * width] =
          cast_ray(Vec3f(camera.x + eyesep / 2, camera.y, camera.z),
                   Vec3f(dir_x, dir_y, dir_z).normalize(), spheres, lights);
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
        pixmap[(j * (width - delta) * 2 + i + width - delta) * 3 + d] =
            255 * c2[d];
      }
    }
  }
  stbi_write_jpg(("out_stereo_" + std::to_string(i) + ".jpg").c_str(),
                 (width - delta) * 2, height, 3, pixmap.data(), 100);
}

#define NB_IMAGES 40

static inline void render_video(
    void (*renderer)(const std::vector<Sphere> &, const std::vector<Light> &,
                     const Vec3f &, const Vec3f &, const int),
    const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
  Vec3f camera{0, 0, 0};


  for (int angle = 0; angle <= NB_IMAGES; ++angle) {
    // TODO: correctly rotate the camera around the center of the image
    
    std::clog << "\033[0G\033[2K[" << angle << "/" << NB_IMAGES
              << "] Generating video...";

    renderer(spheres, lights, camera, Vec3f{0, 0, 0}, angle);

    // On essaie de tourner autour du point {0, 0, -14}
    const Vec3f center{0, 0, -14};

    float theta = 1.f * M_PI / 180.f;
    float x_rot = camera.x * cos(theta) + (camera.z + center.z) * sin(theta);
    float z_rot = -camera.x * sin(theta) + (camera.z + center.z) * cos(theta)
    - center.z;

    /*camera.x -= sin(1./4.);
    camera.y += 0.;
    camera.z += cos(1./4.);*/
    //camera.z += z_rot;
  }
}

static inline void render_video_rebond(
        void (*renderer)(const std::vector<Sphere> &, const std::vector<Light> &,
                         const Vec3f &, const Vec3f &, const int),
        std::vector<Sphere> &spheres, const std::vector<Light> &lights){
    int first = 0;
    for (int angle = 1; angle <= NB_IMAGES+1; ++angle) {

        std::clog << "\033[0G\033[2K[" << angle << "/" << NB_IMAGES
                  << "] Generating video...";

        renderer(spheres, lights, Vec3f {0, 0, 0}, Vec3f{0, 0, 0}, angle);


        switch (first) {
            case 0 :
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

        //changement tous les 5 images
        if(angle % 5 == 0) {
            first = (first + 1) % 8;
        }

    }
}

int main() {
  int n = -1;
  unsigned char *pixmap =
      stbi_load("../envmap.jpg", &envmap_width, &envmap_height, &n, 0);
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
      envmap[dup_comp] =
          Vec3f(pixmap[offset + 0], pixmap[offset + 1], pixmap[offset + 2]) *
          (1 / 255.);
    }
  }
  stbi_image_free(pixmap);

  std::vector<Sphere> spheres;
  spheres.push_back(Sphere(Vec3f(-3, 0, -16), 2, ivory));
  spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2, glass));
  spheres.push_back(Sphere(Vec3f(1.5, -0.5, -18), 3, red_rubber));
  spheres.push_back(Sphere(Vec3f(7, 5, -18), 4, mirror));

  std::vector<Light> lights;
  lights.push_back(Light(Vec3f(-20, 20, 20), 1.5));
  lights.push_back(Light(Vec3f(30, 50, -25), 1.8));
  lights.push_back(Light(Vec3f(30, 20, 30), 1.7));

#if 0
  std::clog << "Rendering normal image..." << std::endl;
  render_normal(spheres, lights, 0);
  std::clog << "Rendering parallax image..." << std::endl;
  render_parallax(spheres, lights, 0);
  std::clog << "Rendering stereoscope image..." << std::endl;
  render_stereoscope(spheres, lights, 0);
#endif
  //render_video(&render_normal, spheres, lights);
  render_video_rebond(&render_normal, spheres, lights);

  return 0;
}
