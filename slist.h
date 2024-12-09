#ifndef SLIST_H
#define SLIST_H

#include "node.h"
#include <string>
using namespace std;

class slist {
public:
    Node* head;

    slist();
    ~slist();
    void add(Airport aVal);
    void clear();
    bool equals(slist lVal);
    Node* get(int val);
    Airport getAirport(int val);
    void insert(int val, Airport aVal);
    void exchg(int val1, int val2);
    void swap(int val1, int val2);
    bool isEmpty();
    void remove(int val);
    void set(int val, Airport aVal);
    int size();
    slist subList(int start, int length);
    string toString();
};

#endif
