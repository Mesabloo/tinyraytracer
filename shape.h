#pragma once
#ifndef SHAPE_H
#define SHAPE_H

#include "interval.h"
#include "model.h"
#include <cfloat>
#include <optional>
#include <set>

struct Shape {
    virtual ~Shape() = default;
    virtual std::set<Interval> ray_intersect(const Vec3f &orig, const Vec3f &dir) const = 0;
};

struct Sphere : public Shape {
    Vec3f center;
    float radius;
    Material material;

    Sphere(const Vec3f &c, const float r, const Material &m) : center(c), radius(r), material(m) {}

    virtual std::set<Interval> ray_intersect(const Vec3f &orig, const Vec3f &dir) const override {
        // https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection

        const Vec3f L = center - orig;

        const float tca = L * dir;
        if (tca < 0)
            return {};

        const float radius2 = radius * radius;
        const float d2 = L * L - tca * tca;

        if (d2 > radius2)
            return {};

        const float thc = sqrt(radius2 - d2);

        float t0 = tca - thc;
        float t1 = tca + thc;

        const std::function<Vec3f(const Vec3f &)> normal = [this](const Vec3f &hit) {
            return (hit - center).normalize();
        };

        if (t0 < 0) {
            if (t1 < 0)
                // pas de point d'intersection : t0 et t1 sont derrière le rayon
                return {};

            // un seul point d'intersection
            const Vec3f hit = orig + dir * t1;
            return {Interval{hit, hit, normal, material, orig}};
        }

        // on veut que le point le plus proche soit `t0`
        if (t0 > t1)
            std::swap(t0, t1);

        const Vec3f from = orig + dir * t0;
        const Vec3f to = orig + dir * t1;
        return {Interval{from, to, normal, material, orig}};
    }
};

struct Checkerboard : public Shape {
    float yAxis;

    Checkerboard(const float y) : yAxis(y) {}

    virtual std::set<Interval> ray_intersect(const Vec3f &orig, const Vec3f &dir) const override {
        if (fabs(dir.y) > 1e-3) {
            const float d = -(orig.y + yAxis) / dir.y;
            const Vec3f pt = orig + dir * d;

            if (d > 0 && fabs(pt.x) < 10 && pt.z < -10 && pt.z > -30) {
                Material material;
                material.diffuse_color =
                    (int(.5 * pt.x + 1000) + int(.5 * pt.z)) & 1 ? Vec3f(.3, .3, .3) : Vec3f(.3, .2, .1);

                return {Interval{pt, pt + 1e-3f, [](const Vec3f &) { return Vec3f(0, 1, 0); }, material, orig}};
            }
        }

        return {};
    }
};

struct DuckObject : public Shape {
    //      _        _        _
    //   __(.)<   __(.)>   __(.)=
    //   \___)    \___)    \___)
    Model duck;

    DuckObject() : duck("../duck.obj") {}
    DuckObject(const char *path) : duck(path) {}

    virtual std::set<Interval> ray_intersect(const Vec3f &orig, const Vec3f &dir) const override {
        std::set<Interval> inter{};

        //      _        _        _
        //   __(.)<   __(.)>   __(.)=
        //   \___)    \___)    \___)
        const int nfaces = duck.nfaces();

#pragma omp parallel for
        for (int i = 0; i < nfaces; ++i) {
            float dist;
            if (duck.ray_triangle_intersect(i, orig, dir, dist)) {
                const Vec3f hit = orig + dir * dist;

                inter.insert(Interval{hit, hit + 1e-3f,
                                      [i, this](const Vec3f &) {
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

                                          return cross(b - a, c - a).normalize();
                                      },
                                      glass, orig});
            }
        }

        return inter;
    }
};

#endif
