#ifndef FUSION_VISO2_IMU_LIBRARY_H
#define FUSION_VISO2_IMU_LIBRARY_H
#include <Eigen/Dense>
#include <iostream>

using namespace std;

class Fusion{
private:
    string name;
public:
    Fusion(string const& name);

    ~Fusion();


    void hello();
};
#endif //FUSION_VISO2_IMU_LIBRARY_H
