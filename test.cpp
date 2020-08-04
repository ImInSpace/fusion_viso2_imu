#include<iostream>
#include "library.h"
#include "multivar_noise.h"
#include <fstream>

using namespace std;

int main() {
    //Setup output format and file
    IOFormat singleLine(FullPrecision,DontAlignCols,", ","; ","","","[","]");
    IOFormat csv(FullPrecision,DontAlignCols,",",",","","","","");
    ofstream file;
    file.open("data.csv", ios::trunc);

    //Initialize kalman filter
    int N = 4;
    double dt = 0.1;
    Fusion f = Fusion(N);
    f.setConstantDt(dt);

    //Randomize x, Q and R
    srand(time(NULL));
    VectorXd x(N);
    x << 10 * VectorXd::Random(N / 2), VectorXd::Random(N / 2);
    MatrixXd Q = MatrixXd::Random(N, N)/100;
    MatrixXd R = MatrixXd::Random(N / 2, N / 2)/10;

    //Make Q and R into positive semi-definite matrices
    Q=Q.transpose()*Q;
    R=R.transpose()*R;

    //Initialize noise generators for v and w
    normal_random_variable v { Q };
    normal_random_variable w { R };


    //Start simulation
    for (int i = 0; i < 100; i++) {
        //Simulate movement with simulated process noise and try to predict it
        x.head(N / 2) += x.tail(N / 2) * dt;
        x += v()*dt;
        f.predict(VectorXd::Zero(4),Q);

        //Update kalman filter with simulated measurement noise
        VectorXd z=x.head(N/2)+w();
        f.update(z,R);

        //Print information
        file<<z.head(2).format(csv)<<","<<f.getx().head(2).format(csv)<<","<<f.getP().topLeftCorner(2,2).format(csv)<<endl;
    }
    file.close();
}
