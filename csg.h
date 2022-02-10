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

#if 0
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

            // pas de chevauchement trouvé : on retire tous les intervalles
            // -- ne rien faire

            for (const auto &[i2, i] : tmp2) {
                // - chevauchement total : on choisit `i2`
                if (i == i2) {
                    result.insert(i2);
                    continue;
                }

                // - inclusion [i2, i, i2] : on garde `i`
                if (i2.includes(i)) {
                    result.insert(i);
                    continue;
                }

                const float dist_i_orig = (i.from - orig).norm();
                const float dist_i2_orig = (i2.from - orig).norm();

                if (dist_i2_orig < dist_i_orig) {
                    // - chevauchement partiel [i2, i] : on garde (i.from, i2.to)
                    result.insert(Interval{i.from, i2.to, i.compute_normal, i.material, orig});
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
                    if (il == ir)
                        break;

                    if (il.includes(ir))
                        break;
                  
                    if (il.overlaps(ir)) {
                        any_overlap = true;

                        const float dist_ir_orig = (ir.from - orig).norm();
                        const float dist_il_orig = (il.from - orig).norm();

                        if (dist_il_orig < dist_ir_orig) {
                            // - chevauchement partiel [il, ir] : on coupe `ir` à l'intersection
                            //   `il.to` et on insère dans le vecteur
                            result.insert(il);
                            result.insert(Interval{il.to, ir.to, ir.compute_normal, ir.material, orig});
                        } else {
                            // - chevauchement partiel [ir, il] : on coupe `il` à l'intersection
                            //   `ir.to` et on insère dans le vecteur
                            result.insert(il);
                            result.insert(Interval{ir.from, il.from, il.compute_normal, il.material, orig});
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
            return result;
        }
        }
#endif
#if 1
        if (node_type == CSGTreeNodeType::INTERSECTION) {
            std::set<Interval> result{};
            std::set<std::pair<Interval, Interval>> tmp2{};

            for (auto const &i1 : left_intersect) {
                for (auto const &i2 : right_intersect) {
                    if (i1.overlaps(i2))
                        tmp2.insert(std::make_pair(i1, i2));
                }
            }

            // pas de chevauchement trouvé : on retire tous les intervalles
            // -- ne rien faire

            for (const auto &[i2, i] : tmp2) {
                // - chevauchement total : on choisit `i2`
                if (i == i2) {
                    result.insert(i2);
                    continue;
                }

                // - inclusion [i2, i, i2] : on garde `i`
                if (i2.includes(i)) {
                    result.insert(i);
                    continue;
                }

                const float dist_i_orig = (i.from - orig).norm();
                const float dist_i2_orig = (i2.from - orig).norm();

                if (dist_i2_orig < dist_i_orig) {
                    // - chevauchement partiel [i2, i] : on garde (i.from, i2.to)
                    result.insert(Interval{i.from, i2.to, i.compute_normal, i.material, orig});
                } else {
                    // - chevauchement partiel [i, i2] : on garde (i2.from, i.to)
                    result.insert(Interval{i2.from, i.to, i2.compute_normal, i2.material, orig});
                }
            }

            return result;
        } else {
            std::set<Interval> result = left_intersect;

            std::for_each(std::cbegin(right_intersect), std::cend(right_intersect),
                          [&result, &orig, this](const Interval &i) {
                              const auto cit = std::find_if(std::cbegin(result), std::cend(result),
                                                            [&i, &orig](const Interval &i2) { return i2.overlaps(i); });

                              switch (node_type) {
                              case CSGTreeNodeType::UNION: {
                                  // pas de chevauchement trouvé : on ajoute `i` en conservant l'ordre
                                  if (cit == std::cend(result)) {
                                      result.insert(i);
                                      return;
                                  }

                                  // chevauchement trouvé : on résout le conflit selon la règle `UNION(i2, i)` donnée
                                  // précédemment
                                  const Interval &i2 = *cit;

                                  // - chevauchement total : les deux intervalles sont les mêmes -> on garde `i2`
                                  if (i == i2) {
                                      return;
                                  }
                                  // - inclusion [i2, i, i2] : on supprime `i` et on garde `i2` (`i` ne sera en principe
                                  // pas visible)
                                  if (i2.includes(i)) {
                                      return;
                                  }

                                  const float dist_i_orig = (i.from - orig).norm();
                                  const float dist_i2_orig = (i2.from - orig).norm();

                                  if (dist_i2_orig < dist_i_orig) {
                                      // - chevauchement partiel [i2, i] : on coupe `i` à l'intersection
                                      //   `i2.to` et on réinsère dans le vecteur
                                      result.insert(Interval{i.from, i2.to, i.compute_normal, i.material, orig});
                                  } else {
                                      // - chevauchement partiel [i, i2] : on coupe `i2` à l'intersection
                                      //   `i.to` et on réinsère dans le vecteur
                                      result.erase(cit);
                                      result.insert(Interval{i2.from, i.to, i2.compute_normal, i2.material, orig});
                                      result.insert(i);
                                  }

                                  break;
                              }
                              case CSGTreeNodeType::DIFFERENCE: {
                                  // pas de chevauchement trouvé : on garde les intervalles de base
                                  if (cit == std::cend(result))
                                      return;

                                  // chevauchement trouvé : on résout le conflit selon la règle `DIFFERENCE(i2, i)`
                                  // donnée précédemment
                                  const Interval i2 = *cit;
                                  result.erase(cit);

                                  // - chevauchement total : on enlève `i2` complètement
                                  if (i == i2) {
                                      return;
                                  }
                                  // - inclusion [i2, i, i2] : on sépare i2 en i3=(i2.from, i.from) et i4=(i.to, i2.to)
                                  // et on
                                  //   réinsère [i3, i4] dans le vecteur
                                  if (i2.includes(i)) {
                                      result.insert(Interval{i2.from, i.from, i2.compute_normal, i2.material, orig});
                                      result.insert(Interval{i.to, i2.to,
                                                             [&i](const Vec3f &hit) { return -i.compute_normal(hit); },
                                                             i2.material, orig});
                                      return;
                                  }

                                  const float dist_i_orig = (i.from - orig).norm();
                                  const float dist_i2_orig = (i2.from - orig).norm();

                                  if (dist_i2_orig < dist_i_orig) {
                                      // - chevauchement partiel [i2, i] : on garde (i2.from, i.from)
                                      result.insert(Interval{i2.from, i.from, i2.compute_normal, i2.material, orig});
                                  } else {
                                      // - chevauchement partiel [i, i2] : on garde (i.to, i2.to)
                                      result.insert(Interval{i.to, i2.to,
                                                             [&i](const Vec3f &hit) { return -i.compute_normal(hit); },
                                                             i2.material, orig});
                                  }

                                  break;
                              }
                              default: {
                              }
                              }
                          });

            return result;
        }
#endif

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
