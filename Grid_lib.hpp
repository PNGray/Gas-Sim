#pragma once
#include<omp.h>
#include<math.h>
#include "Vec_lib.hpp"
#include "Linked_lists.hpp"

const double r0 = 0.15;
const double r02 = r0 * r0;
const double thresh = pow(20, 2) * r02;
const double m = 1;
int getzone(V3d &pos, iV3d &b_side, double gridsize, V3d &origin){
  V3d v = pos - origin;
  int x = (int)(v.x / gridsize);
  int y = (int)(v.y / gridsize);
  int z = (int)(v.z / gridsize);
  return x * b_side.x + y * b_side.y + z;
}

void check_grid(sLink *box, V3d *pos, int *zone, iV3d &b_side, double gridsize, V3d &origin){
  #pragma omp parallel for
  for (int i = 0; i < b_side.z; i++){
    sLink *current = box + i;
    sLink *target;
    while (current->next != nullptr) {
      int tag = current->next->val;
      int correct_zone = getzone(pos[tag], b_side, gridsize, origin);
      if (correct_zone != zone[tag]){
        zone[tag] = correct_zone;
        target = current->remove_next();
        if (correct_zone >=0 && correct_zone <b_side.z){
          #pragma omp critical
          {
            box[correct_zone].add(target);
          }
        }
      }
      current = current->next;
      if (current == nullptr) break;
    }
  }
}

void init_ps_links(sLink *pss, double n){
  for (int i = 0; i < n; i++){
    (pss + i)->val = i;
  }
}

void init_grid(sLink *box, sLink *pss, V3d *pos, int *zone, iV3d &b_side, double gridsize, V3d &origin, int n){
  #pragma omp parallel for
  for (int i = 0; i < n; i++){
    int z = getzone(pos[i], b_side, gridsize, origin);
    zone[i] = z;
    sLink *p = pss + i;
    #pragma omp critical
    {
      box[z].add(p);
    }
  }
}

double V(double r2, double gridsize){
  double g2 = gridsize * gridsize;
  if (r2 > g2) return V(g2 , gridsize);
  r2 /= r02;
  double r6 = r2 * r2 * r2;
  return (4 * r0 * (1 / r6 / r6) - 1 / r6);
}

double temperature(V3d *vs, int n){
  double ke = 0;
  V3d avg_v(0);
  for (int i = 0; i < n; i++){
    avg_v.add(*(vs + i));
  }
  avg_v.div(n);
  for (int i = 0; i < n; i++){
    V3d rel_v = *(vs + i) - avg_v;
    ke += 0.5 * m * rel_v.lensqr();
  }
  return ke;
}

double potential_energy(V3d *ps, int n, double gridsize){
  double pe = 0;
  double r2;
  V3d dist;
  for (int i = 0; i < n; i++){
    #pragma omp parallel for private(r2, dist)
    for (int j = 0; j < n; j++){
      if (i != j) {
        dist = *(ps + i) - *(ps + j);
        r2 = dist.lensqr();
        #pragma omp critical
        {
          pe += V(r2, gridsize);
        }
      }
    }
  }
  return pe;
}
void gas_force(V3d &p1, V3d &p2, V3d &a, int &num_inter, double gridsize){
  V3d r = p1 - p2;
  double rlen2 = r.lensqr();
  if (rlen2 > gridsize * gridsize) return;
  double runit2 = rlen2 / r02;
  double runit6 = runit2 * runit2 * runit2;
  double runit8 = runit2 * runit6;
  double runit14 = runit6 * runit8;
  double force = 24 * ((2 / runit14) - (1 / runit8)) / r0;
  if (runit2 < thresh) num_inter++;
  r.mul(force);
  a.add(r);
}

void apply_grid(sLink *grid, int index, V3d *pos, V3d *acc, int *num_inter, double gridsize){
  sLink *current = grid->next;
  while (current != nullptr){
    if (current->val != index) {
      gas_force(pos[index], pos[current->val], acc[index], num_inter[index], gridsize);
    }
    current = current->next;
  }
}

void apply(sLink *box, V3d *pos, V3d *acc, int *num_inter, int index, int *zone, iV3d &b_side, double gridsize){
  int current;
  num_inter[index] = 0;
  for (int i = -1; i <= 1; i++){
    for (int j = -1; j <=1; j++){
      for (int k = -1; k <=1; k++){
        current = zone[index] + b_side.x * i + b_side.y * j + k;
        if (current >= 0 && current < b_side.z){
          apply_grid(box + current, index, pos, acc, num_inter, gridsize);
        }
      }
    }
  }
}
