#pragma once
#include<stdio.h>
#include<math.h>
#include<tuple>

class V3d {
public:
  double x;
  double y;
  double z;

  V3d();

  V3d(double);

  V3d(double, double, double);

  void add(V3d&);

  void sub(V3d&);

  void mul(double);

  void div(double);

  //dot product
  double operator*(V3d&);

  //cross product
  V3d operator^(V3d&);

  V3d operator+(V3d&);

  V3d operator-(V3d&);

  double lensqr();

  double len();

  int quadrant(V3d&);

  void show();

  void reset();
};

V3d operator*(double, V3d&);

int generate(V3d*, V3d*, double*, FILE*);

class iV3d{
public:
  int x;
  int y;
  int z;

  iV3d();

  iV3d(int);

  iV3d(int, int, int);

  bool operator==(iV3d&);
};
