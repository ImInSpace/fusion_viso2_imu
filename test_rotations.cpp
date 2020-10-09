//
// Created by Inigo on 08/10/2020.
//
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace std;
using namespace Eigen;

int main(){
    VectorXd x(7);
    //x = VectorXd::Random(7);

    IOFormat singleLine(StreamPrecision, DontAlignCols, ",\t", ";\t", "", "", "[", "]");
    Quaterniond q= Quaterniond::Identity();
    cout<<q.coeffs().format(singleLine)<<endl;
}

