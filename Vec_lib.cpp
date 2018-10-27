#include "Vec_lib.hpp"

V3d::V3d(){
  x = 0;
  y = 0;
  z = 0;
}

V3d::V3d(double val){
  x = val;
  y = val;
  z = val;
}

V3d::V3d(double x_, double y_, double z_){
  x = x_;
  y = y_;
  z = z_;
}

void V3d::add(V3d &v){
  x += v.x;
  y += v.y;
  z += v.z;
}

void V3d::sub(V3d &v){
  x -= v.x;
  y -= v.y;
  z -= v.z;
}

void V3d::mul(double c){
  x *= c;
  y *= c;
  z *= c;
}

void V3d::div(double c){
  x /= c;
  y /= c;
  z /= c;
}

double V3d::operator*(V3d &v){
  return x * v.x + y * v.y + z * v.z;
}

V3d V3d::operator^(V3d &v){
  double x3 = y * v.z - v.y * z;
  double y3 = z * v.x - v.z * x;
  double z3 = x * v.y - v.x * y;
  V3d v3(x3, y3, z3);
  return v3;
}

V3d V3d::operator+(V3d &v){
  V3d v3(x + v.x, y + v.y, z + v.z);
  return v3;
}

V3d V3d::operator-(V3d &v){
  V3d v3(x - v.x, y - v.y, z - v.z);
  return v3;
}

double V3d::lensqr(){
  return x * x + y * y + z * z;
}

double V3d::len(){
  return sqrt(lensqr());
}

int V3d::quadrant(V3d &v){
  if (z < v.z){
    if (x < v.x){
      if (y < v.y){
        return 7;
      }
      else {
        return 6;
      }
    }
    else {
      if (y < v.y) {
        return 8;
      }
      else {
        return 5;
      }
    }
  }
  else {
    if (x < v.x){
      if (y < v.y){
        return 3;
      }
      else {
        return 2;
      }
    }
    else {
      if (y < v.y) {
        return 4;
      }
      else {
        return 1;
      }
    }
  }
}

void V3d::show(){
  printf("V3d: x: %f y:%f z: %f\n", x, y, z);
}

void V3d::reset(){
  x = 0;
  y = 0;
  z = 0;
}

V3d operator*(double c, V3d &v) {
  V3d v2(c * v.x, c * v.y, c * v.z);
  return v2;
}

iV3d::iV3d(){
  x = 0;
  y = 0;
  z = 0;
}

iV3d::iV3d(int val){
  x = val;
  y = val;
  z = val;
}

iV3d::iV3d(int valx, int valy, int valz){
  x = valx;
  y = valy;
  z = valz;
}

bool iV3d::operator==(iV3d &v){
  return x == v.x && y == v.y && z == v.z;
}

int generate(V3d* ps, V3d* vs, double* ms, double n, FILE *file){
  double x, y, z, vx, vy, vz, m;
  V3d *p, *v;
  int scan_result = 0, index = 0;
  for (int i = 0; i < n; i++){
    fscanf(file, "%lf", &x);
    fscanf(file, "%lf", &y);
    fscanf(file, "%lf", &z);
    fscanf(file, "%lf", &vx);
    fscanf(file, "%lf", &vy);
    fscanf(file, "%lf", &vz);
    scan_result = fscanf(file, "%lf", &m);
      p = ps + index;
      v = vs + index;
      p->x = x;
      p->y = y;
      p->z = z;
      v->x = vx;
      v->y = vy;
      v->z = vz;
      ms[index] = m;
      index++;
  }
  return index;
}
