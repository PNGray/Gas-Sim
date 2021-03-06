#include<omp.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<random>
#include "Vec_lib.hpp"
#include "Linked_lists.hpp"
#include "Grid_lib.hpp"

const double dt = 0.001;
const double dthalf = dt / 2;
extern const double r0half;

void init_vel(V3d *vel, int n, double env_vel) {
  std::default_random_engine gen;
  std::normal_distribution<double> distri(0, 1);
  double factor = env_vel / sqrt(3);
  for (int i = 0; i < n; i++) {
    vel[i].x = distri(gen) * factor;
    vel[i].y = distri(gen) * factor;
    vel[i].z = distri(gen) * factor;
  }
}

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
  if (argc < 13) {
    printf("infile outfile savefile datafile pad sizexy sizez conduct envvel numstep spf numthreads\n");
    return 0;
  }
  //inputs variables
  //file io
  char* infilename = argv[1];
  char* outfilename = argv[2];
  char* savefilename = argv[3];
  char* datafilename = argv[4];
  //env var
  double pad = atof(argv[5]);
  double sizexy = atof(argv[6]);
  double sizez = atof(argv[7]);
  double conduct = atof(argv[8]);
  double env_vel = atof(argv[9]);
  //time var
  int num_step = atoi(argv[10]);
  int steps_per_frame = atoi(argv[11]);
  omp_set_num_threads(atoi(argv[12]));

  //file io
  FILE *infile = fopen(infilename, "r");
  FILE *outfile = strncmp(outfilename, "stdout", 10) == 0 ? stdout :
    (strncmp(outfilename, "null", 5) == 0 ? nullptr : fopen(outfilename, "w"));
  FILE *data = strncmp(datafilename, "stdout", 10) == 0 ? stdout :
    (strncmp(datafilename, "null", 5) == 0 ? nullptr : fopen(datafilename, "w"));

  //prints arguments data
  fprintf(data,"#");
  for (int i = 0; i < argc; i++) fprintf(data, "%s ", argv[i]);
  fprintf(data, "\n");
  fflush(data);

  //particle variables
  int n = count_len(infile, sizexy);
  V3d ps[n];
  V3d vs[n];
  V3d as[n];
  double ms[n];
  int num_inter[n];
  int zone[n];
  double p_rad = 0.075;
  double particle_volume = 4.0 / 3.0 * M_PI * p_rad * p_rad * p_rad * n;

  //animation variables
  bool color = false;
  int blob = 5;

  //box variables
  double bound_pad = 1;
  V3d origin(-sizexy / 2, -sizexy / 2, -bound_pad / 2);
  int sidexy = (int)(sizexy / pad / r0);
  double gridsize = sizexy / sidexy;
  double gs3 = gridsize * gridsize * gridsize;
  int sidez = (int)(sizez / gridsize);
  iV3d b_side(sidez * sidexy, sidez, sidexy * sidexy * sidez);
  double boundxy = (sizexy - bound_pad) / 2.0;
  double boundz = sizez - bound_pad / 2.0;
  double rboundxy = boundxy - r0half, rboundz = boundz - r0half;
  sLink *box = new sLink[b_side.z];
  sLink pss[n];
  printf("%f %f %d %d %f\n", boundxy, boundz, sidexy, sidez, gridsize);

  //environment variable
  double env_ke = n / 2.0 * env_vel * env_vel;
  double gravity = -0.98;
  double e = 0;
  double ke = 0;
  double pe = 0;
  double energy_added = 0;
  double impulse = 0;
  double sample_t = 0;
  double sample_pe = 0;
  double pressure = 0;
  int p_sample = 5000;
  int t_sample = p_sample / steps_per_frame;
  double area = 8 * boundxy * boundxy + 8 * boundxy * boundz;
  double volume = 4 * boundxy * boundxy * boundz;
  double PV = 0;
  double NkT = 0;
  double ratio = 0;
  double n_per_v = n / volume; //density
  double beta = 0;

  //binding energy calculation
  double g2 = gridsize * gridsize;
  double g2u = g2 / r02;
  double g6u = g2u * g2u * g2u;
  double zero_point = (4 * r0 * (1 / g6u / g6u - 1 / g6u));
  printf("%f %f\n", zero_point * n, env_ke);

  //initialization
  init_ps_links(pss, n);
  generate(ps, vs, ms, n, infile);
  init_grid(box, pss, ps, zone, b_side, gridsize, origin, n);
  // init_vel(vs, n, env_vel);

  //main loop
  for (int i = 0; i <= num_step; i++){
    double t = i * dt;

    // #pragma omp parallel for
    for (int i = 0; i < n; i++){
      V3d d = dthalf * vs[i];
      energy_added += gravity * d.z;
      ps[i].add(d);
    }

    check_grid(box, ps, zone, b_side, n,  gridsize, origin);

    #pragma omp parallel for
    for (int i = 0; i < n; i++){
      apply(box, ps, as, num_inter, i, zone, b_side, gridsize);
    }

    if (i % steps_per_frame == 0)
      sample_pe += potential_energy(box, ps, zone, b_side, n, g2, zero_point);

    // compress
    // if (t > 0) {
    //   if (boundxy > 0.8) {boundxy -= 0.00005; rboundxy = boundxy - r0half;}
    //   if (boundz > 1.26953125) {boundz -= 0.0002; rboundz = boundz - r0half;}
    // }
    // if (t == 20) conduct /= 2;
    // if (ke < 11600 && t > 10) conduct = 0;

    // #pragma omp parallel for
    for (int i = 0; i < n; i++){
      V3d *p = ps + i;
      double x = fabs(p->x);
      double y = fabs(p->y);
      double f;

      as[i].z += gravity;

      if (x > rboundxy /* - 0.01 * r0 */){
        f = 10000 * (rboundxy - x);
        impulse -= f;
        f /= x / p->x;
        as[i].x += f;
      }

      if (y > rboundxy /* - 0.01 * r0 */){
        f = 10000 * (rboundxy - y);
        impulse -= f;
        f /= y / p->y;
        as[i].y += f;
      }

      if (p->z < r0half /* + 0.01 * r0 */){
        f = 10000 * (r0half - p->z);
        impulse += f;
        as[i].z += f;
        // double oldke = 0.5 * vs[i].lensqr();

        double len = vs[i].len();
        double factor = (env_vel - len) * conduct; // inject heat if bounce of bottom
        V3d dv = factor / len * vs[i];
        vs[i].add(dv);

        // double newke =  0.5 * vs[i].lensqr();
        // #pragma omp critical
        // {
        //   energy_added += newke - oldke;
        // }
      }

      if (p->z > rboundz /* - 0.01 * r0 */) {
        f = 10000 * (rboundz - p->z);
        impulse -= f;
        as[i].z += f;
        // double oldke = 0.5 * vs[i].lensqr();
        // double len = vs[i].len();
        // double factor = (0.2 - len) * 0.2;
        // V3d dv = factor * vs[i];
        // vs[i].add(dv);
        // double newke =  0.5 * vs[i].lensqr();
        // #pragma omp critical
        // {
        //   energy_added += newke - oldke;
        // }
      }

    }
    // if (ke > env_ke || t < 25)
    // for (int i = 0; i < n; i++){
    //   double oldke = 0.5 * vs[i].lensqr();
    //   double len = vs[i].len();
    //   double factor = (env_vel - len) * conduct;
    //   V3d dv = factor * vs[i];
    //   vs[i].add(dv);
    //   double newke =  0.5 * vs[i].lensqr();
    //   energy_added += newke - oldke;
    // }

    #pragma omp parallel for
    for (int i = 0; i < n; i++){
      as[i].mul(dt);
      vs[i].add(as[i]);
      as[i].reset();
    }

    // #pragma omp parallel for
    for (int i = 0; i < n; i++){
      V3d d = dthalf * vs[i];
      energy_added += gravity * d.z;
      ps[i].add(d);
    }

    //calculate average pressure
    if (i % p_sample == 0) {
      area =  8 * boundxy * boundxy + 8 * boundxy * boundz;
      pressure = impulse / p_sample / area;
      impulse = 0;

      ke = sample_t / t_sample;
      sample_t = 0;

      pe = sample_pe / t_sample;
      sample_pe = 0;
    }

    if (i % steps_per_frame == 0){
      volume = 4 * boundxy * boundxy * boundz;
      n_per_v = n / volume;
      sample_t += kinetic_energy(vs, n);
      e = ke + pe;
      PV = pressure * volume;
      NkT = ke * 2.0 / 3.0;
      beta = -pressure * (volume - particle_volume) + NkT;
      ratio = PV / NkT;
      if (data != nullptr) fprintf(data, "%f %f %f %f %f %f %f\n", t, ke, pe, energy_added, e, ratio, beta);

      printf("%d\n", i);
      if (outfile == nullptr) continue;
      for (int j = 0; j < n; j++){
        V3d *p = ps + j;
        if (num_inter[j] > blob && !color) {
          fprintf(outfile, "C 0 1 1\n");
          color = true;
        }
        else if (num_inter[j] < blob && color) {
          fprintf(outfile, "C 1 1 1\n");
          color = false;
        }
        fprintf(outfile, "c3 %f %f %f 0.07\n", p->x, p->y, p->z);
      }
      draw_box(boundxy, boundz, outfile);
      fprintf(outfile, "T -0.8 0.8\nt = %.2f\te = %f\tP = %f\tV = %f\nT -0.8 0.7\nPV = %f\t NkT = %f\t ratio = %f\n", dt * i, e, pressure, volume, PV, NkT, ratio);
      fprintf(outfile, "T -0.8 0.6\nKE = %f\tPE = %f\tbeta = %f\n", ke, pe, beta);
      fprintf(outfile, "F\n");
    }
  }
  fclose(infile);

  //saving current condition to file
  double newsize = boundxy * 2 + bound_pad;
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
