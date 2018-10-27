#include "Linked_lists.hpp"

sLink::sLink(){
  val = -1;
  next = nullptr;
}

sLink::sLink(int v) {
  val = v;
  next = nullptr;
}

sLink::sLink(int v, sLink* l){
  val = v;
  next = l;
}

void sLink::add(int v) {
  sLink *new_link = new sLink(v, next);
  next = new_link;
}

void sLink::add(sLink* l) {
  l->next = next;
  next = l;
}
sLink* sLink::remove_next(){
  if (next == nullptr) return nullptr;
  sLink *n = next;
  next = next->next;
  n->next = nullptr;
  return n;
}

sLink* sLink::remove(int v) {
  sLink *current = this;
  while (current->next != nullptr){
    if (current->next->val == v){
      return current->remove_next();
    }
    current = current->next;
  }
  return nullptr;
}

dLink::dLink(){
  val = -1;
  prev = this;
  next = this;
}

dLink::dLink(int v){
  val = v;
  prev = this;
  next = this;
}

void dLink::add(int val){
  dLink *new_link = new dLink(val);
  dLink *nextnext = next;
  next = new_link;
  new_link->prev = this;
  new_link->next = nextnext;
  nextnext->prev = new_link;
}

void dLink::add(dLink* new_link){
  dLink *nextnext = next;
  next = new_link;
  new_link->prev = this;
  new_link->next = nextnext;
  nextnext->prev = new_link;
}

void dLink::remove_self(){
  next->prev = prev;
  prev->next = next;
  next = this;
  prev = this;
}

dLink* dLink::remove(int v){
  dLink *current = next;
  while (current != this){
    if (current->val == v){
      current->remove_self();
      return current;
    }
  }
  return nullptr;
}

void dLink::concat(dLink* l) {
  if (l->next == l) return;
  prev->next = l->next;
  l->next->prev = prev;
  prev = l->prev;
  l->prev->next = this;
  l->next = l;
  l->prev = l;
}
