#pragma once

class sLink{
public:
  int val;
  sLink *next;

  sLink();

  sLink(int);

  sLink(int, sLink*);

  void add(int);

  void add(sLink*);

  sLink* remove_next();

  sLink* remove(int);
};

class dLink{
public:
  int val;
  dLink *prev;
  dLink *next;

  dLink();

  dLink(int);

  void add(int);

  void add(dLink*);

  void remove_self();

  dLink* remove(int);

  void concat(dLink*);
};
