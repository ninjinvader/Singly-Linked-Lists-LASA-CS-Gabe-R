#include "slist.h"
#include <iostream>
#include <stdexcept>
using namespace std;

slist::slist() {
    head = nullptr;
}

slist::~slist() {
    clear();
}

void slist::swap(int index1, int index2) {
    if (index1 == index2) return;

    Node* node1 = get(index1);
    Node* node2 = get(index2);

    if (node1 == nullptr || node2 == nullptr) {
        throw out_of_range("Index out of bounds");
    }

  
    Airport temp = node1->data;
    node1->data = node2->data;
    node2->data = temp;
}

bool slist::equals(slist a) {
    Node *current = head;
    Node *aCurrent = a.head;
    while (current != nullptr && aCurrent != nullptr) {
        if (current->data.code != aCurrent->data.code || 
            current->data.longitude != aCurrent->data.longitude || 
            current->data.latitude != aCurrent->data.latitude) {
            return false;
        }
        current = current->next;
        aCurrent = aCurrent->next;
    }
    return (current == nullptr && aCurrent == nullptr);
}

void slist::add(Airport newA) {
    Node *pntr = new Node(newA);
    pntr->next = nullptr;

    if (head == nullptr) {
        head = pntr;
    } else {
        Node *current = head;
        while (current->next != nullptr) {
            current = current->next;
        }
        current->next = pntr;
    }
}


Node* slist::get(int index) {
    Node *current = head;
    int count = 0;

    while (current != nullptr && count < index) {
        current = current->next;
        count++;
    }
    return (count == index) ? current : nullptr;
}


Airport slist::getAirport(int index) {
    Node* node = get(index);  
    if (node != nullptr) {
        return node->data;  
    } else {
        throw std::out_of_range("Index out of bounds");
    }
}

void slist::clear() {
    Node *next;
    while (head != nullptr) {
        next = head->next;
        delete head;
        head = next;
    }
}

void slist::insert(int index, Airport aVal) {
    if (index < 0) throw out_of_range("Index cannot be negative");
    if (index == 0) {
        Node* newNode = new Node(aVal);
        newNode->next = head;
        head = newNode;
        return;
    }

    Node* current = head;
    for (int i = 0; current != nullptr && i < index - 1; ++i) {
        current = current->next;
    }

    if (current == nullptr) {
        throw out_of_range("Index is out of bounds");
    }

    Node* newNode = new Node(aVal);
    newNode->next = current->next;
    current->next = newNode;
}


bool slist::isEmpty() {
    return head == nullptr;
}


void slist::remove(int index) {
    if (index < 0 || head == nullptr) throw out_of_range("Index out of bounds");

    Node* temp;
    if (index == 0) {
        temp = head;
        head = head->next;
        delete temp;
        return;
    }

    Node* current = head;
    for (int i = 0; current != nullptr && i < index - 1; ++i) {
        current = current->next;
    }

    if (current == nullptr || current->next == nullptr) {
        throw out_of_range("Index out of bounds");
    }

    temp = current->next;
    current->next = current->next->next;
    delete temp;
}


void slist::set(int index, Airport aVal) {
    Node* node = get(index);
    if (node != nullptr) {
        node->data = aVal;
    } else {
        throw out_of_range("Index out of bounds");
    }
}


int slist::size() {
    int count = 0;
    Node* current = head;
    while (current != nullptr) {
        count++;
        current = current->next;
    }
    return count;
}

slist slist::subList(int start, int length) {
    slist sublist;
    Node* current = get(start);
    int count = 0;

    while (current != nullptr && count < length) {
        sublist.add(current->data);
        current = current->next;
        count++;
    }
    return sublist;
}

string slist::toString() {
    string result = "[";
    Node* current = head;
    while (current != nullptr) {
        result += current->data.toString(); // Assuming `Airport` has a `toString` method
        if (current->next != nullptr) {
            result += ", ";
        }
        current = current->next;
    }
    result += "]";
    return result;
}