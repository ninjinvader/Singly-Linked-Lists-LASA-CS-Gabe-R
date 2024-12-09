#include "airport.h"
#include <iostream>

#ifndef NODE_H
#define NODE_H

class Node
{
public:
    Airport data;
    Node *next;
    Node()
    {
        next = NULL;
    }
    Node(Airport newA)
    {
        data = newA;
    }
};
#endif
