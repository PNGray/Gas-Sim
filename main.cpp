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
      if (dump2 != 0 && size < dump2) size = dump2;
      break;
    }
  }
  rewind(file);
  return count;
}

void draw_box(double xy, double z, FILE *outfile) {
  // baaaaarf.
  fprintf(outfile,"l3 -%e -%e  0  -%e -%e  %e\n",xy,xy,xy,xy,z);
  fprintf(outfile,"l3 -%e -%e  %e -%e  %e  %e\n",xy,xy,z,xy,xy,z);
  fprintf(outfile,"l3 -%e  %e  %e -%e  %e  0 \n",xy,xy,z,xy,xy);
  fprintf(outfile,"l3 -%e  %e  0  -%e -%e  0 \n",xy,xy,xy,xy);
  fprintf(outfile,"l3  %e -%e  0   %e -%e  %e\n",xy,xy,xy,xy,z);
  fprintf(outfile,"l3  %e -%e  %e  %e  %e  %e\n",xy,xy,z,xy,xy,z);
  fprintf(outfile,"l3  %e  %e  %e  %e  %e  0 \n",xy,xy,z,xy,xy);
  fprintf(outfile,"l3  %e  %e  0   %e -%e  0 \n",xy,xy,xy,xy);
  fprintf(outfile,"l3  %e -%e  0  -%e -%e  0 \n",xy,xy,xy,xy);
  fprintf(outfile,"l3  %e  %e  0  -%e  %e  0 \n",xy,xy,xy,xy);
  fprintf(outfile,"l3  %e -%e  %e -%e -%e  %e\n",xy,xy,z,xy,xy,z);
  fprintf(outfile,"l3  %e  %e  %e -%e  %e  %e\n",xy,xy,z,xy,xy,z);
}

int main(int argc, char **argv){
  if (argc < 8) {
    printf("infile outfile savefile pad size numstep numthreads\n");
    return 0;
  }
  //inputs variables
  char* infilename = argv[1];
  char* outfilename = argv[2];
  char* savefilename = argv[3];
  double pad = atof(argv[4]);
  double sizexy = atof(argv[5]);
  double sizez = atof(argv[6]);
  int num_step = atoi(argv[7]);
  omp_set_num_threads(atoi(argv[8]));

  //file io
  FILE *infile = fopen(infilename, "r");
  FILE *outfile = strncmp(outfilename, "stdout", 10) == 0 ? stdout :
    (strncmp(outfilename, "null", 5) == 0 ? nullptr : fopen(outfilename, "w"));
  FILE *kinetic = fopen("kinetic.txt", "w");

  //particle variables
  int n = count_len(infile, sizexy);
  V3d ps[n];
  V3d vs[n];
  V3d as[n];
  double ms[n];
  int num_inter[n];
  int zone[n];

  //color variables
  bool color = false;
  int blob = 5;

  //environment variable
  double gravity = -0.1;
  double e = 0;
  double conduct = 0.015;
  double env_vel = 2.5;
  double last_ke = 0;
  double current_ke = 0;
  double pe = 0;
  double energy_added = 0;

  //box variables
  V3d origin(-sizexy / 2, -sizexy / 2, 0);
  int sidexy = (int)(sizexy / pad / r0);
  double gridsize = sizexy / sidexy;
  int sidez = (int)(sizez / gridsize);
  iV3d b_side(sidez * sidexy, sidez, sidexy * sidexy * sidez);
  double boundxy = (sizexy - gridsize) / 2.0;
  double boundz = sizez - gridsize;
  sLink *box = new sLink[b_side.z];
  printf("%f %f %d %d\n", boundxy, boundz, sidexy, sidez);
  // sLink box[blen];
  sLink pss[n];

  //initialization
  init_ps_links(pss, n);
  generate(ps, vs, ms, n, infile);
  init_grid(box, pss, ps, zone, b_side, gridsize, origin, n);

  //main loop
  for (int i = 0; i <= num_step; i++){
    double t = i * dt;

    #pragma omp parallel for
    for (int i = 0; i < n; i++){
      V3d d = dthalf * vs[i];
      ps[i].add(d);
    }

    check_grid(box, ps, zone, b_side, gridsize, origin);

    #pragma omp parallel for
    for (int i = 0; i < n; i++){
      apply(box, ps, as, num_inter, i, zone, b_side, gridsize);
    }

    // if (boundxy > 0.7) boundxy -= 0.0001;
    #pragma omp parallel for
    for (int i = 0; i < n; i++){
      V3d *p = ps + i;
      double x = fabs(p->x);
      double y = fabs(p->y);
      double f;

      /* if (num_inter[i] > blob)  */as[i].z += gravity;

      if (x > boundxy){
        f = 50000 * (boundxy - x) / x * p->x;
        as[i].x += f;
      }

      if (y > boundxy){
        f = 50000 * (boundxy - y) / y * p->y;
        as[i].y += f;
      }

      if (p->z < 0){
        f = 50000 * (-p->z);
        as[i].z += f;
        double oldke = 0.5 * vs[i].lensqr();
        double len = vs[i].len();
        double factor = (env_vel - len) * conduct;
        V3d dv = factor * vs[i];
        vs[i].add(dv);
        double newke =  0.5 * vs[i].lensqr();
        #pragma omp critical
        {
          energy_added += newke - oldke;
        }
      }

      if (p->z > boundz) {
        f = 50000 * (boundz - p->z);
        as[i].z += f;
      }

      // if (num_inter[i] > blob) {
      //   double oldke = 0.5 * vs[i].lensqr();
      //   vs[i].mul(1.0001);
      //   double newke = 0.5 * vs[i].lensqr();
      //   energy_added += newke - oldke;
      // }
    }


    #pragma omp parallel for
    for (int i = 0; i < n; i++){
      as[i].mul(dt);
      vs[i].add(as[i]);
      as[i].reset();
    }

    #pragma omp parallel for
    for (int i = 0; i < n; i++){
      V3d d = dthalf * vs[i];
      ps[i].add(d);
    }

    if (i % 100 == 0){
      // pe = potential_energy(ps, n, gridsize);
      current_ke= temperature(vs, n);
      e = current_ke + pe;
      fprintf(kinetic, "%f %f %f %f\n", t, current_ke, e, energy_added);
      last_ke = current_ke;

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
      draw_box(boundxy, boundz, outfile);
      fprintf(outfile, "T -0.8 0.8\nt = %.2f\te = %f\n", dt * i, e);
      fprintf(outfile, "F\n");
    }
  }
  fclose(infile);

  //saving current condition to file
  double newsize = boundxy * 2 + gridsize;
  FILE *savefile = strncmp(savefilename, "stdout", 10) == 0 ? stdout :
    (strncmp(savefilename, "null", 5) == 0 ? nullptr : fopen(savefilename, "w"));
  if (savefile != nullptr) {
    for (int i = 0; i < n; i++){
      V3d *p = ps + i;
      V3d *v = vs + i;
      fprintf(savefile, "%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n", p->x, p->y, p->z, v->x, v->y, v->z, ms[i]);
    }
    fprintf(savefile, "%.15lf\n", newsize);
  }
  printf("new size is %f\n", newsize);
}
