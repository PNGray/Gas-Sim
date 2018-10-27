#include<omp.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "Vec_lib.hpp"
#include "Linked_lists.hpp"
#include "Grid_lib.hpp"

const double dt = 0.001;
const double dthalf = dt / 2;
int count_len(FILE *file, double &size){
  double dump;
  double dump2 = 0;
  int scan_result = 0, count = 0;
  while (true){
    scan_result = fscanf(file, "%lf %lf %lf %lf %lf %lf %lf", &dump2, &dump, &dump, &dump, &dump, &dump, &dump);
    if (scan_result == 7){
      dump2 = 0;
      count++;
    }
    else {
      if (dump2 != 0) size = dump2;
      break;
    }
  }
  rewind(file);
  return count;
}

void draw_box(double L, FILE *outfile) {
  // baaaaarf.
  fprintf(outfile,"l3 -%e -%e -%e -%e -%e  %e\n",L,L,L,L,L,L);
  fprintf(outfile,"l3 -%e -%e  %e -%e  %e  %e\n",L,L,L,L,L,L);
  fprintf(outfile,"l3 -%e  %e  %e -%e  %e -%e\n",L,L,L,L,L,L);
  fprintf(outfile,"l3 -%e  %e -%e -%e -%e -%e\n",L,L,L,L,L,L);
  fprintf(outfile,"l3  %e -%e -%e  %e -%e  %e\n",L,L,L,L,L,L);
  fprintf(outfile,"l3  %e -%e  %e  %e  %e  %e\n",L,L,L,L,L,L);
  fprintf(outfile,"l3  %e  %e  %e  %e  %e -%e\n",L,L,L,L,L,L);
  fprintf(outfile,"l3  %e  %e -%e  %e -%e -%e\n",L,L,L,L,L,L);
  fprintf(outfile,"l3  %e -%e -%e -%e -%e -%e\n",L,L,L,L,L,L);
  fprintf(outfile,"l3  %e  %e -%e -%e  %e -%e\n",L,L,L,L,L,L);
  fprintf(outfile,"l3  %e -%e  %e -%e -%e  %e\n",L,L,L,L,L,L);
  fprintf(outfile,"l3  %e  %e  %e -%e  %e  %e\n",L,L,L,L,L,L);
}

int main(int argc, char **argv){
  if (argc < 8) {
    printf("infile outfile savefile pad size numstep numthreads\n");
    return 0;
  }

  char* infilename = argv[1];
  char* outfilename = argv[2];
  char* savefilename = argv[3];
  double pad = atof(argv[4]);
  double size = atof(argv[5]);
  int num_step = atoi(argv[6]);
  omp_set_num_threads(atoi(argv[7]));

  FILE *infile = fopen(infilename, "r");
  FILE *outfile = strncmp(outfilename, "stdout", 10) == 0 ? stdout :
    (strncmp(outfilename, "null", 5) == 0 ? nullptr : fopen(outfilename, "w"));
  FILE *kinetic = fopen("kinetic.txt", "w");
  int n = count_len(infile, size);

  V3d ps[n];
  V3d vs[n];
  V3d as[n];
  int num_inter[n];
  bool color = false;
  int blob = 5;
  V3d origin(-size / 2);
  double gravity = -0.1;
  double ms[n];
  double e = 0;
  int zone[n];

  int side = (int)(size / pad / r0);
  iV3d b_side(side * side, side, side * side * side);
  int blen = side * side * side;
  double gridsize = size / side;
  double bound = (size - gridsize) / 2.0;
  double conduct = 0.05;
  double env_vel = 1.5;

  sLink *box = new sLink[blen];
  // sLink box[blen];
  sLink pss[n];
  init_ps_links(pss, n);
  generate(ps, vs, ms, n, infile);
  init_grid(box, pss, ps, zone, b_side, gridsize, origin, n);

  for (int i = 0; i <= num_step; i++){
    double t = i * dt;

    #pragma omp parallel for
    for (int i = 0; i < n; i++){
      V3d d = dthalf * vs[i];
      ps[i].add(d);
    }

    check_grid(box, ps, zone, b_side, gridsize, origin, blen);

    #pragma omp parallel for
    for (int i = 0; i < n; i++){
      apply(box, ps, as, num_inter, i, zone, b_side, gridsize);
    }

    #pragma omp parallel for
    for (int i = 0; i < n; i++){
      as[i].z += gravity;
      V3d *p = ps + i;
      double x = fabs(p->x);
      double y = fabs(p->y);
      double z = fabs(p->z);
      double f;
      if (x > bound){
        f = 50000 * (bound - x) / x * p->x;
        as[i].x += f;
        // double len = vs[i].len();
        // V3d dv = (env_vel - len) * conduct * vs[i];
        // vs[i].add(dv);
      }
      if (y > bound){
        f = 50000 * (bound - y) / y * p->y;
        as[i].y += f;
        // double len = vs[i].len();
        // V3d dv = (env_vel - len) * conduct * vs[i];
        // vs[i].add(dv);
      }
      if (z > bound){
        f = 50000 * (bound - z) / z * p->z;
        as[i].z += f;
        double len = vs[i].len();
        if (ps[i].z < 0){
          V3d dv = (env_vel - len) * conduct * vs[i];
          vs[i].add(dv);
        }
      }
      // vs[i].mul(1.00005);
    }


    #pragma omp parallel for
    for (int i = 0; i < n; i++){
      as[i].mul(dt);
      vs[i].add(as[i]);
      as[i].reset();
    }

    fprintf(kinetic, "%f %f\n", t, kinetic_energy(vs, n));

    #pragma omp parallel for
    for (int i = 0; i < n; i++){
      V3d d = dthalf * vs[i];
      ps[i].add(d);
    }

    if (i % 100 == 0){
      // if (i % 200 == 0) e = energy(ps, vs, n, gridsize);
      printf("%d\n", i);
      if (outfile == nullptr) continue;
      for (int j = 0; j < n; j++){
        V3d *p = ps + j;
        if (j == 0){
          fprintf(outfile, "C 1 0 0\n");
          fprintf(outfile, "ct3 0 %f %f %f 0.1\n", p->x, p->y, p->z);
        }
        else {
          if (j == 1) {
            fprintf(outfile, "C 1 1 1\n");
            color = false;
          }
          if (num_inter[j] > blob && !color) {
            fprintf(outfile, "C 0 1 1\n");
            color = true;
          }
          else if (num_inter[j] < blob && color) {
            fprintf(outfile, "C 1 1 1\n");
            color = false;
          }
          fprintf(outfile, "c3 %f %f %f 0.1\n", p->x, p->y, p->z);
        }
      }
      draw_box(bound, outfile);
      fprintf(outfile, "T -0.8 0.8\nt = %.2f\te = %f\n", dt * i, e);
      fprintf(outfile, "F\n");
    }
  }
  fclose(infile);
  FILE *savefile = strncmp(savefilename, "stdout", 10) == 0 ? stdout :
    (strncmp(savefilename, "null", 5) == 0 ? nullptr : fopen(savefilename, "w"));
  if (savefile != nullptr) {
    for (int i = 0; i < n; i++){
      V3d *p = ps + i;
      V3d *v = vs + i;
      fprintf(savefile, "%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n", p->x, p->y, p->z, v->x, v->y, v->z, ms[i]);
    }
    fprintf(savefile, "%.15lf\n", size);
  }
}
