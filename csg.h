#pragma once
#ifndef CSG_H
#define CSG_H

#include <memory>
#include <set>
#include <tuple>
#include <vector>

#include "geometry.h"
#include "interval.h"
#include "material.h"
#include "shape.h"

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
    virtual std::set<Interval> ray_intersect(const Vec3f &orig, const Vec3f &dir) const = 0;
};

template <class T, typename = std::enable_if_t<std::is_base_of_v<Shape, T>>> struct CSGTreeNode : public CSGTree<T> {
    const std::shared_ptr<CSGTree<T>> left, right;
    const CSGTreeNodeType node_type;

    CSGTreeNode(const std::shared_ptr<CSGTree<T>> &left, const std::shared_ptr<CSGTree<T>> &right,
                const CSGTreeNodeType node_type)
        : left(left), right(right), node_type(node_type) {}

    virtual std::set<Interval> ray_intersect(const Vec3f &orig, const Vec3f &dir) const override {
        const std::set<Interval> left_intersect = left->ray_intersect(orig, dir);
        const std::set<Interval> right_intersect = right->ray_intersect(orig, dir);

        // Si l'arbre de CSG gauche est vide, et qu'on calcule autre chose qu'une
        // union, alors on peut simplement retourner un ensemble d'intervalles vide.
        // Si on calcule une union, alors on retourne la composante de droite.
        if (left_intersect.empty())
            return node_type == CSGTreeNodeType::UNION ? right_intersect : std::set<Interval>{};
        // Si on fait une union et que les intersections à droite sont vides, alors
        // on retourne les intersections à gauche.
        if (node_type == CSGTreeNodeType::UNION && right_intersect.empty())
            return left_intersect;

        // Fusion de deux intervalles :
        //
        // - UNION(I1, I2) :
        //   - si deux intervalles se chevauchent :
        //     - chevauchement total : les deux intervalles sont les mêmes -> on en prend un au hasard (???)
        //     - chevauchement partiel [I1, I2] : on coupe aléatoirement I1 ou I2 à l'intersection (I1.to ou
        //     I2.from)
        //       entre les deux et on réinsère dans le vecteur
        //     - inclusion [I1, I2, I1] : on garde I1 et I2 (I2 peut être visible à travers du verre)
        //     - inclusion [I2, I1, I2] : on garde I1 et I2 (I1 peut être visible à travers du verre)
        //   - sinon :
        //     - on réinsère les deux intervalles dans le vecteur
        // - DIFFERENCE(I1, I2) :
        //   - si deux intervalles se chevauchent :
        //     - chevauchement total : on enlève I1 complètement
        //     - chevauchement partiel [I1, I2] : on garde (I1.from, I2.from)
        //     - inclusion [I1, I2, I1] : on sépare I1 en I3=(I1.from, I2.from) et I4=(I2.to, I1.to) et on réinsère
        //       [I3;I4] dans le vecteur
        //     - inclusion [I2, I1, I2] : on retire I1
        //   - sinon :
        //     - on garde I1 entier
        // - INTERSECTION(I1, I2) :
        //   - si deux intervalles se chevauchent :
        //     - chevauchement total : on choisit I1 ou I2 aléatoirement (???)
        //     - chevauchement partiel [I1, I2] : on garde (I2.from, I1.to)
        //     - inclusion [I1, I2, I1] : on garde I2
        //   - sinon :
        //     - on enlève les deux intervalles

        std::set<Interval> result{};
        switch (node_type) {
        case CSGTreeNodeType::INTERSECTION: {
            std::set<std::pair<Interval, Interval>> tmp2{};

            for (auto const &i1 : left_intersect) {
                for (auto const &i2 : right_intersect) {
                    if (i1.overlaps(i2))
                        tmp2.insert(std::make_pair(i1, i2));
                }
            }

            for (const auto &[i2, i] : tmp2) {
                // - chevauchement total : on choisit `i2`
                if (i == i2) {
                    result.insert(i2);
                    continue;
                }

                // - inclusion [i2, i, i2] : on garde `i`
                if (i2.includes(i)) {
                    Interval i3 = i;
                    i3.material = i2.material;
                    result.insert(i3);
                    continue;
                }

                const float dist_i_orig = (i.from - orig).norm();
                const float dist_i2_orig = (i2.from - orig).norm();

                if (dist_i2_orig <= dist_i_orig) {
                    // - chevauchement partiel [i2, i] : on garde (i.from, i2.to)
                    result.insert(Interval{i.from, i2.to, i.compute_normal, i2.material, orig});
                } else {
                    // - chevauchement partiel [i, i2] : on garde (i2.from, i.to)
                    result.insert(Interval{i2.from, i.to, i2.compute_normal, i2.material, orig});
                }
            }

            return result;
        }
        case CSGTreeNodeType::UNION: {
            for (const Interval &il : left_intersect) {
                bool any_overlap = false;

                for (const Interval &ir : right_intersect) {
                    // - chevauchement total : les deux intervalles sont les mêmes -> on garde `il`
                    if (il == ir)
                        break;

                    // - inclusion [il, ir, il] : on garde les deux intervalles au cas où `ir` soit
                    //   visible (à travers du verre par exemple)
                    // - inclusion [ir, il, ir] : on garde les deux intervalles au cas où `il` soit
                    //   visible (à travers du verre par exemple)
                    if (il.includes(ir) || ir.includes(il)) {
                        result.insert(ir);
                        break;
                    }

                    if (il.overlaps(ir) || ir.overlaps(il)) {
                        any_overlap = true;

                        const float dist_ir_orig = (ir.from - orig).norm();
                        const float dist_il_orig = (il.from - orig).norm();

                        if (dist_il_orig <= dist_ir_orig) {
                            // - chevauchement partiel [il, ir] : on coupe `ir` à l'intersection
                            //   `il.to` et on insère dans le vecteur
                            result.insert(il);
                            result.insert(Interval{il.to, ir.to, ir.compute_normal, ir.material, orig});
                        } else {
                            // - chevauchement partiel [ir, il] : on coupe `il` à l'intersection
                            //   `ir.to` et on insère dans le vecteur
                            result.insert(il);
                            result.insert(Interval{ir.from, il.from, ir.compute_normal, ir.material, orig});
                        }

                        break;
                    }
                }

                if (!any_overlap)
                    result.insert(il);
            }
            return result;
        }
        case CSGTreeNodeType::DIFFERENCE: {
            for (const Interval &il : left_intersect) {
                bool any_overlap = false;

                for (const Interval &ir : right_intersect) {
                    // - chevauchement total : on enlève `ir` complètement
                    if (il == ir) {
                        any_overlap = true;
                        break;
                    }

                    // - inclusion [il, ir, il] : on sépare il en i3=(il.from, ir.from) et i4=(ir.to, il.to)
                    //   et on réinsère [i3, i4] dans le vecteur
                    if (il.includes(ir)) {
                        any_overlap = true;

                        const Interval i1{il.from, ir.from, il.compute_normal, il.material, orig};
                        const Interval i2{ir.to, il.to, [ir](const Vec3f &hit) { return -ir.compute_normal(hit); },
                                          il.material, orig};

                        result.insert(i1);
                        result.insert(i2);
                        break;
                    }

                    // - inclusion [ir, il, ir] : on supprime `il`
                    if (ir.includes(il)) {
                        any_overlap = true;
                        break;
                    }

                    if (il.overlaps(ir) || ir.overlaps(il)) {
                        any_overlap = true;

                        const float dist_ir_orig = (ir.from - orig).norm();
                        const float dist_il_orig = (il.from - orig).norm();

                        if (dist_il_orig <= dist_ir_orig) {
                            // - chevauchement partiel [il, ir] : on garde (il.from, ir.from)
                            result.insert(Interval{il.from, ir.from, il.compute_normal, il.material, orig});
                        } else {
                            // - chevauchement partiel [ir, il] : on garde (ir.to, il.to)
                            result.insert(Interval{ir.to, il.to,
                                                   [ir](const Vec3f &hit) { return -ir.compute_normal(hit); },
                                                   il.material, orig});
                        }

                        break;
                    }
                }

                if (!any_overlap)
                    result.insert(il);
            }

            return result;
        }
        }

        return {};
    }
};

template <class T, typename = std::enable_if_t<std::is_base_of_v<Shape, T>>> struct CSGTreeLeaf : public CSGTree<T> {
    const T &value;

    CSGTreeLeaf(const T &val) : value(val) {}

    virtual std::set<Interval> ray_intersect(const Vec3f &orig, const Vec3f &dir) const override {
        return value.ray_intersect(orig, dir);
    }
};

#endif
