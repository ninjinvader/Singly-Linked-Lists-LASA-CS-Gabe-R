#include <iostream>
#include <fstream>
#include <cmath>
#include "slist.h"

using namespace std;


#define pi 3.14159265358979323846
#define earthRadiusKm 6371.0

// Function declarations
void mergeSort(slist* s);
Node* merge(Node* left, Node* right);
Node* mergeSortRecursive(Node* head);
Node* getMiddle(Node* head); 
bool compare(Airport a, Airport b);
double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d); 
void within100Miles(slist* airports, double refLat, double refLon) {
    Node* current = airports->head;
    while (current != nullptr) {
        double distKm = distanceEarth(current->data.latitude, current->data.longitude, refLat, refLon);
        double dist = distKm * 0.621371;
        if (dist <= 100.0) {
            cout << current->data.code << " is <= 100 miles away. The distance is " << dist << " miles." << endl;
        }
        current = current->next;
    }
}
int main()
{
    ifstream infile;
    int i = 0;
    char cNum[10];
    slist airports;
    int airportNum;

    infile.open("./USAirportCodes.csv", ifstream::in);
    if (infile.is_open())
    {
        int c = 0;
        while (infile.good())
        {
            Airport current;
            infile.getline(current.code, 256, ',');
            infile.getline(cNum, 256, ',');
            current.latitude = atof(cNum);
            infile.getline(cNum, 256, '\n');
            current.longitude = atof(cNum);

            i++;
            c++;
            airports.add(current);
        }
        airportNum = c;
        infile.close();
    }


    
   mergeSort(&airports);

    for (int c = 0; c < airportNum; c++)
    {
       double lat =  airports.getAirport(c).latitude;
       double longi =  airports.getAirport(c).longitude;
       cout << airports.getAirport(c).code << " long: " << airports.getAirport(c).longitude
             << " lat: " << airports.getAirport(c).latitude << " dis: " << distanceEarth(30.1944, 97.6700,lat,longi) << endl;
    }
    cout << endl;
    cout << airports.getAirport(airportNum-3).code << " long: " << airports.getAirport(airportNum-3).longitude << " lat: " << airports.getAirport(airportNum-3).latitude << " dis: " << distanceEarth(airports.getAirport(airportNum-3).latitude, airports.getAirport(airportNum-3).longitude, 30.1944, 97.6700) << endl ;
    within100Miles(&airports, 30.1944, 97.6700);
}

Node* getMiddle(Node* head)
{
    if (!head)
        return nullptr;

    Node* slow = head;
    Node* fast = head;

    while (fast->next && fast->next->next)
    {
        slow = slow->next;
        fast = fast->next->next;
    }

    return slow;
}

bool compare(Airport a, Airport b)
{
    double latitude = 30.1944;
    double longitude = 97.6700;
    return distanceEarth(a.latitude, a.longitude, latitude, longitude) < distanceEarth(b.latitude, b.longitude, latitude, longitude);
}

void mergeSort(slist* airport)
{
    airport->head = mergeSortRecursive(airport->head);
}

Node* mergeSortRecursive(Node* head)
{

    if (!head || !head->next)
        return head;

    // Found this technique online finding middle node
    Node* center = getMiddle(head);
    Node* center1 = center->next;
    center->next = nullptr; 
    Node* l = mergeSortRecursive(head);
    Node* r = mergeSortRecursive(center1);

    return merge(l, r);
}

Node* merge(Node* l, Node* r)
{

    if (!l)
        return r;
    if (!r)
        return l;

    if (compare(l->data, r->data))
    {
        l->next = merge(l->next, r);
        return l;
    }
    else
    {
        r->next = merge(l, r->next);
        return r;
    }
}


double deg2rad(double deg) {
  return (deg * pi / 180);
}

double rad2deg(double rad) {
  return (rad * 180 / pi);
}

/**
 * Returns the distance between two points on the Earth.
 * Direct translation from http://en.wikipedia.org/wiki/Haversine_formula
 * @param lat1d Latitude of the first point in degrees
 * @param lon1d Longitude of the first point in degrees
 * @param lat2d Latitude of the second point in degrees
 * @param lon2d Longitude of the second point in degrees
 * @return The distance between the two points in kilometers
 */
double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d) {
  double lat1r, lon1r, lat2r, lon2r, u, v;
  lat1r = deg2rad(lat1d);
  lon1r = deg2rad(lon1d);
  lat2r = deg2rad(lat2d);
  lon2r = deg2rad(lon2d);
  u = sin((lat2r - lat1r)/2);
  v = sin((lon2r - lon1r)/2);
  return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}
