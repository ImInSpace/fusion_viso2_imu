#include<iostream>
#include "library.h"
using namespace std;
int main(){
    string name = "a";
    Fusion a = Fusion(name);
    name = "b";
    Fusion b = Fusion(name);
    a.hello();
    b.hello();
}
