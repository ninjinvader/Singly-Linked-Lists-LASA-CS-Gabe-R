#include <string>
using namespace std;

class Airport{
    public:
        char    code[5];
        double   longitude;
        double   latitude;
    
        bool operator == (Airport other){
            return this->code == other.code && this->longitude == other.longitude && this->latitude==other.latitude;
        }

        string toString(){
            return this->code;
        }

};
