#include "library.h"
using namespace std;

void Fusion::hello() {
    cout << "Hello, World, I am "<<name<<"!" << endl;
}

Fusion::Fusion(string const& name) {
    this->name=name;
}

Fusion::~Fusion() {

}
