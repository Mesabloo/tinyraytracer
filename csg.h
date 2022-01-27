#pragma once
#ifndef CSG_H
#define CSG_H

#include <memory>
#include <tuple>
#include <vector>

#include "geometry.h"
#include "interval.h"
#include "material.h"
#include "shape.h"

namespace std {
template <class InputIt1, class InputIt2, class OutputIt>
inline void cross(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, OutputIt out) {
    for (auto f1it = first1; f1it != last1; f1it++) {
        size_t i = 0;

        while (first2 + i != last2) {
            *out++ = std::make_pair(*f1it, *(first2 + i++));
        }
    }
}
} // namespace std

static inline void interval_order_insert(std::vector<Interval> &vec, const Vec3f &orig, const Interval &i) {
    auto cit = std::begin(vec);
    const auto cend = std::end(vec);

    while (cit != cend) {
        const Interval &i2 = *cit;
        const Vec3f &from = i.from;
        const Vec3f &from2 = i2.from;

        // si la distance `from <-> orig` < la distance `from2 <-> orig` alors on insère `i` avant `i2`
        const float norm_i_orig = (from - orig).norm();
        const float norm_i2_orig = (from2 - orig).norm();
        if (norm_i_orig < norm_i2_orig) {
            vec.insert(cit, i);
            return;
        }

        cit++;
    }
    // si on ne l'a pas inséré au milieu, on le met à la fin
    vec.push_back(i);
}

enum class CSGTreeNodeType {
    //! Union de deux primitives biaisée par la gauche (`P1 ∪ P2` favorise `P1`)
    UNION,
    //! Différence entre deux primitives (`P1 - P2`)
    DIFFERENCE,
    //! Intersection entre deux primitives (`P1 ∩ P2`)
    INTERSECTION
};

template <class T, typename = std::enable_if_t<std::is_base_of_v<Shape, T>>> struct CSGTree {
    virtual ~CSGTree() = default;
    virtual std::vector<Interval> ray_intersect(const Vec3f &orig, const Vec3f &dir) const = 0;
};

template <class T, typename = std::enable_if_t<std::is_base_of_v<Shape, T>>> struct CSGTreeNode : public CSGTree<T> {
    const std::shared_ptr<CSGTree<T>> left, right;
    const CSGTreeNodeType node_type;

    CSGTreeNode(const std::shared_ptr<CSGTree<T>> &left, const std::shared_ptr<CSGTree<T>> &right,
                const CSGTreeNodeType node_type)
        : left(left), right(right), node_type(node_type) {}

    virtual std::vector<Interval> ray_intersect(const Vec3f &orig, const Vec3f &dir) const override {
        const std::vector<Interval> left_intersect = left->ray_intersect(orig, dir);
        const std::vector<Interval> right_intersect = right->ray_intersect(orig, dir);

        // Si l'arbre de CSG gauche est vide, et qu'on calcule autre chose qu'une
        // union, alors on peut simplement retourner un ensemble d'intervalles vide.
        // Si on calcule une union, alors on retourne la composante de droite.
        if (left_intersect.empty())
            return node_type == CSGTreeNodeType::UNION ? right_intersect : std::vector<Interval>{};
        // Si on fait une union et que les intersections à droite sont vides, alors
        // on retourne les intersections à gauche.
        if (node_type == CSGTreeNodeType::UNION && right_intersect.empty())
            return left_intersect;

        // Fusion de deux intervalles :
        //
        // - UNION(I1, I2) :
        //   - si deux intervalles se chevauchent :
        //     - chevauchement total : les deux intervalles sont les mêmes -> on en prend un au hasard (???)
        //     - chevauchement partiel [I1, I2] : on coupe aléatoirement I1 ou I2 à l'intersection (I1.to ou I2.from)
        //       entre les deux et on réinsère dans le vecteur
        //     - inclusion [I1, I2, I1] : on supprime I2 et on garde I1 (I2 ne sera en principe pas visible)
        //   - sinon :
        //     - on réinsère les deux intervalles dans le vecteur
        // - DIFFERENCE(I1, I2) :
        //   - si deux intervalles se chevauchent :
        //     - chevauchement total : on enlève I1 complètement
        //     - chevauchement partiel [I1, I2] : on garde (I1.from, I2.from)
        //     - inclusion [I1, I2, I1] : on sépare I1 en I3=(I1.from, I2.from) et I4=(I2.to, I1.to) et on réinsère
        //       [I3;I4] dans le vecteur
        //   - sinon :
        //     - on garde I1 entier
        // - INTERSECTION(I1, I2) :
        //   - si deux intervalles se chevauchent :
        //     - chevauchement total : on choisit I1 ou I2 aléatoirement (???)
        //     - chevauchement partiel [I1, I2] : on garde (I2.from, I1.to)
        //     - inclusion [I1, I2, I1] : on garde I2
        //   - sinon :
        //     - on enlève les deux intervalles

        // Fusion de deux vecteurs d'intervalles (ordonnés) :
        //
        // -

        if (node_type == CSGTreeNodeType::INTERSECTION) {
            std::vector<Interval> result{};
            std::vector<std::pair<Interval, Interval>> tmp{}, tmp2{};

            std::cross(std::cbegin(left_intersect), std::cend(left_intersect), std::cbegin(right_intersect),
                       std::cend(right_intersect), std::back_inserter(tmp));
            std::copy_if(std::cbegin(tmp), std::cend(tmp), std::back_inserter(tmp2),
                         [](const auto &pair) { return std::get<0>(pair).overlaps(std::get<1>(pair)); });

            // pas de chevauchement trouvé : on retire tous les intervalles
            // -- ne rien faire

            for (const auto &[i2, i] : tmp2) {
                // - chevauchement total : on choisit `i2`
                if (i == i2) {
                    interval_order_insert(result, orig, i2);
                    continue;
                }

                // - inclusion [i2, i, i2] : on garde `i`
                if (i2.includes(i)) {
                    interval_order_insert(result, orig, i);
                    continue;
                }

                const float dist_i_orig = (i.from - orig).norm();
                const float dist_i2_orig = (i2.from - orig).norm();

                if (dist_i2_orig < dist_i_orig) {
                    // - chevauchement partiel [i2, i] : on garde (i.from, i2.to)
                    interval_order_insert(result, orig, Interval{i.from, i2.to, i.compute_normal, i.material});
                } else {
                    // - chevauchement partiel [i, i2] : on garde (i2.from, i.to)
                    interval_order_insert(result, orig, Interval{i2.from, i.to, i2.compute_normal, i2.material});
                }
            }

            return result;
        } else {
            std::vector<Interval> result = left_intersect;

            std::for_each(
                std::cbegin(right_intersect), std::cend(right_intersect), [&result, &orig, this](const Interval &i) {
                    const auto cit = std::find_if(std::cbegin(result), std::cend(result),
                                                  [&i, &orig](const Interval &i2) { return i2.overlaps(i); });

                    switch (node_type) {
                    case CSGTreeNodeType::UNION: {
                        // pas de chevauchement trouvé : on ajoute `i` en conservant l'ordre
                        if (cit == std::cend(result)) {
                            interval_order_insert(result, orig, i);
                            return;
                        }

                        // chevauchement trouvé : on résout le conflit selon la règle `UNION(i2, i)` donnée précédemment
                        const Interval &i2 = *cit;

                        // - chevauchement total : les deux intervalles sont les mêmes -> on garde `i2`
                        if (i == i2) {
                            return;
                        }
                        // - inclusion [i2, i, i2] : on supprime `i` et on garde `i2` (`i` ne sera en principe pas
                        // visible)
                        if (i2.includes(i)) {
                            return;
                        }

                        const float dist_i_orig = (i.from - orig).norm();
                        const float dist_i2_orig = (i2.from - orig).norm();

                        if (dist_i2_orig < dist_i_orig) {
                            // - chevauchement partiel [i2, i] : on coupe `i` à l'intersection
                            //   `i2.to` et on réinsère dans le vecteur
                            interval_order_insert(result, orig, Interval{i.from, i2.to, i.compute_normal, i.material});
                        } else {
                            // - chevauchement partiel [i, i2] : on coupe `i2` à l'intersection
                            //   `i.to` et on réinsère dans le vecteur
                            result.erase(cit);
                            interval_order_insert(result, orig,
                                                  Interval{i2.from, i.to, i2.compute_normal, i2.material});
                            interval_order_insert(result, orig, i);
                        }

                        break;
                    }
                    case CSGTreeNodeType::DIFFERENCE: {
                        // pas de chevauchement trouvé : on garde les intervalles de base
                        if (cit == std::cend(result))
                            return;

                        // chevauchement trouvé : on résout le conflit selon la règle `DIFFERENCE(i2, i)` donnée
                        // précédemment
                        const Interval i2 = *cit;
                        result.erase(cit);

                        // - chevauchement total : on enlève `i2` complètement
                        if (i == i2) {
                            return;
                        }
                        // - inclusion [i2, i, i2] : on sépare i2 en i3=(i2.from, i.from) et i4=(i.to, i2.to) et on
                        //   réinsère [i3, i4] dans le vecteur
                        if (i2.includes(i)) {
                            interval_order_insert(result, orig,
                                                  Interval{i2.from, i.from, i2.compute_normal, i2.material});
                            interval_order_insert(result, orig,
                                                  Interval{i.to, i2.to,
                                                           [&i](const Vec3f &hit) { return -i.compute_normal(hit); },
                                                           i2.material});
                            return;
                        }

                        const float dist_i_orig = (i.from - orig).norm();
                        const float dist_i2_orig = (i2.from - orig).norm();

                        if (dist_i2_orig < dist_i_orig) {
                            // - chevauchement partiel [i2, i] : on garde (i2.from, i.from)
                            interval_order_insert(result, orig,
                                                  Interval{i2.from, i.from, i2.compute_normal, i2.material});
                        } else {
                            // - chevauchement partiel [i, i2] : on garde (i.to, i2.to)
                            interval_order_insert(result, orig,
                                                  Interval{i.to, i2.to,
                                                           [&i](const Vec3f &hit) { return -i.compute_normal(hit); },
                                                           i2.material});
                        }

                        break;
                    }
                    default: {
                    }
                    }
                });

            return result;
        }

        return {};
    }
};

template <class T, typename = std::enable_if_t<std::is_base_of_v<Shape, T>>> struct CSGTreeLeaf : public CSGTree<T> {
    const T &value;

    CSGTreeLeaf(const T &val) : value(val) {}

    virtual std::vector<Interval> ray_intersect(const Vec3f &orig, const Vec3f &dir) const override {
        return value.ray_intersect(orig, dir);
    }
};

#endif
