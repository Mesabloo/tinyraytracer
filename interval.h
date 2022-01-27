#pragma once
#ifndef INTERVAL_H
#define INTERVAL_H

#include "geometry.h"
#include "material.h"
#include <cfloat>
#include <functional>

//! Vérifie si le point `p3` se situe entre les points `p1` et `p2`.
static inline constexpr bool is_point_between(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3) {
    // voir https://math.stackexchange.com/a/3210793
    //
    // Soit `t = (x3 - x1)/(x2 - x1)`, on vérifie que `t = (y3 - y1)/(y2 - y1)` et `t = (z3 - z1)/(z2 - z1)`
    // et que `t ∈ [0, 1]`.

    const float tx = (p3.x - p1.x) / (p2.x - p1.x);
    const float ty = (p3.y - p1.y) / (p2.y - p1.y);
    const float tz = (p3.z - p1.z) / (p2.z - p1.z);

    return tx - ty < FLT_EPSILON && tx - tz < FLT_EPSILON && tx >= 0 && tx <= 1;
}

struct Interval {
    Vec3f from, to;
    std::function<Vec3f(const Vec3f &)> compute_normal;
    Material material;

    Interval(const Vec3f &from, const Vec3f &to, const std::function<Vec3f(const Vec3f &)> &normal,
             const Material &material)
        : from(from), to(to), compute_normal(normal), material(material) {}

    //! Est-ce que les deux intervalles se chevauchent ? (en supposant que les deux sont sur le même segment)
    bool overlaps(const Interval &i) const {
        // On vérifie si au moins une des deux extrémités de `i` est incluse dans l'intervalle
        return ::is_point_between(from, to, i.from) || ::is_point_between(from, to, i.to);
    }

    //! Est-ce que l'intervalle donné est complètement inclus dans l'intervalle courant ?
    // (en supposant que les deux sont sur le même segment)
    bool includes(const Interval &i) const {
        // On vérifie si les deux points de `i` sont entre les deux points de l'intervalle
        // - si oui, `i` est inclus dans l'intervalle
        // - si non, `i` n'est pas inclus dans l'intervalle

        return ::is_point_between(from, to, i.from) && ::is_point_between(from, to, i.to);
    }
};

//! Vérifie si les deux intervalles commencent et finissent aux mêmes points
inline bool operator==(const Interval &i1, const Interval &i2) { return i1.from == i2.from && i1.to == i2.to; }

#endif
