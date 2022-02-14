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

static inline bool solve_quadratic(const float &a, const float &b, const float &c, float &x0, float &x1) 
{ 
    float discr = b * b - 4 * a * c; 
    if (discr < 0) return false; 
    else if (discr == 0) x0 = x1 = - 0.5 * b / a; 
    else { 
        float q = (b > 0) ? 
            -0.5 * (b + sqrt(discr)) : 
            -0.5 * (b - sqrt(discr)); 
        x0 = q / a; 
        x1 = c / q; 
    } 
    if (x0 > x1) std::swap(x0, x1); 
 
    return true; 
} 

struct Sphere : public Shape {
    Vec3f center;
    float radius, radius2;
    Material material;

    Sphere(const Vec3f &c, const float r, const Material &m) : center(c), radius(r), radius2(r * r), material(m) {}

    virtual std::set<Interval> ray_intersect(const Vec3f &orig, const Vec3f &dir) const override {
        // https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection

        const Vec3f L = center - orig;

#if 0
        const float a = dir * dir;
        const float b = 2 * (L * dir);
        const float c = L * L - radius2;

        float t0, t1;
        if (!::solve_quadratic(a, b, c, t0, t1))
            return {};
        
#else
        const float tca = L * dir;
        if (tca < 0)
            return {};

        const float d2 = L * L - tca * tca;

        if (d2 > radius2)
            return {};

        const float thc = sqrt(radius2 - d2);

        float t0 = tca - thc;
        float t1 = tca + thc;

        // on veut que le point le plus proche soit `t0`
        if (t0 > t1)
            std::swap(t0, t1);
#endif

        const std::function<Vec3f(const Vec3f &)> normal = [this](const Vec3f &hit) {
            return (hit - center).normalize();
        };

        if (t0 < 0) {
            t0 = t1;
            if (t0 < 0)
                // pas de point d'intersection : t0 et t1 sont derrière le rayon
                return {};

            // un seul point d'intersection
            const Vec3f from = orig + dir * t0;
            const Vec3f to = orig + dir * t1;
            return {Interval{from, to, normal, material, orig}};
        } else {
            const Vec3f from = orig + dir * t0;
            const Vec3f to = orig + dir * t1;
            return {Interval{from, to, normal, material, orig}};
        }
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

                return {Interval{pt, pt, [](const Vec3f &) { return Vec3f(0, 1, 0); }, material, orig}};
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
                const Vec3f normal = cross(b - a, c - a).normalize();

                inter.insert(Interval{hit, hit, [normal](const Vec3f &) { return normal; }, glass, orig});
            }
        }

        return inter;
    }
};

#endif
