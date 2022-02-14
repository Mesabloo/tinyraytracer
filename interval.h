#pragma once
#ifndef INTERVAL_H
#define INTERVAL_H

#include "geometry.h"
#include "material.h"
#include <cfloat>
#include <functional>

//! Vérifie si le point `p3` se situe entre les points `p1` et `p2`.
static inline bool is_point_between(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3) {
    // https://stackoverflow.com/a/33156422/6718698
    const Vec3f ba = p2 - p3;
    const Vec3f bc = p2 - p1;
    const float t = (ba * bc) / (bc * bc);
    return 0 < t && t < 1;
}

struct Interval {
    Vec3f from, to;
    std::function<Vec3f(const Vec3f &)> compute_normal;
    Material material;
    const Vec3f &orig;

    Interval(const Vec3f &from, const Vec3f &to, const std::function<Vec3f(const Vec3f &)> &normal,
             const Material &material, const Vec3f &orig)
        : from(from), to(to), compute_normal(normal), material(material), orig(orig) {}

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

inline bool operator<(const Interval &i1, const Interval &i2) {
    const float dist1 = (i1.orig - i1.from).norm();
    const float dist2 = (i2.orig - i2.from).norm();

    return dist1 < dist2;
}

inline std::ostream &operator<<(std::ostream &os, const Interval &i) {
    return os << "[" << i.from << ", " << i.to << "]";
}

#endif
