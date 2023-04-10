// Basic code that solves a 2nd order differential equation
// using 4th order Runge Kutta with constant step size

// This particular example solves a damped harmonic oscillator
// but can be easily adapted to any 2nd order ode.

// The differential equation being solved is:
//  d²x/dt² + 2 delta_0 * omega_0 * dx/dt + omega_o² * x = 0
//
//  where:
//      omega_0 := undamped angular frequency of harmonic oscillator
//      delta_0 := damping ratio
//      x := position of oscillating particle as a function of time
//      t := time
//
//  We can define a new quantity: y = dx/dt
//  Therefore: dy/dt = d²x/dt²
//
//  So we have to solve two equations:
//  First: y = dx/dt
//  Second: dy/dt = -2 delta_0 * omega_0 * y - omega_0² * x
//
//  These two equations are represented by functions:
//  First: f
//  Second: g


#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;


double f(double t, double x, double y);
double g(double t, double x, double y);
void runge_kutta4(double h, double t, double &xi, double &yi, double(*func1)(double, double, double), double(*func2)(double, double, double));
void creating_output_file(vector <double> &time, vector <double> &solution1, vector <double> &solution2, int n);


int main(){

    //Begin code by defining initial conditions and vectors
    //to store solutions:

    //Beggining time
    double a {0.0};
    //Ending time
    double b {5.0};
    //Number of steps
    int N {1000};

    //Stepsize
    double step_size {(b-a)/N};
    
    //Initial conditions
    double xi {1};
    double yi {1};
    double time {a};

    //Vectors to store solutions
    vector<double> solution1 {xi};
    vector<double> solution2 {yi};
    vector<double> time_vec {time};


    for(int i = 0; i < N; i++){
        cout << xi << endl;
        time_vec.push_back(time);
        time += i*step_size;
        runge_kutta4(step_size, time, xi, yi, &f, &g);
        solution1.push_back(xi);
        solution2.push_back(yi);

    }

    creating_output_file(time_vec, solution1, solution2, N);


    return 0;
}





double f(double t, double x, double y){
    return y;
}

double g(double t, double x, double y){
    double omega_0 {10};
    double delta_0 {0.01};

    return -2*delta_0*omega_0*y - pow(omega_0, 2)*x;
}


void runge_kutta4(double h, double t, double &xi, double &yi, double(*func1)(double, double, double), double(*func2)(double, double, double)){

    double k0 {h*func1(t, xi, yi)};
    double l0 {h*func2(t, xi, yi)};
    double k1 {h*func1(t + 0.5*h, xi + 0.5*k0, yi + 0.5*l0)};
    double l1 {h*func2(t + 0.5*h, xi + 0.5*k0, yi + 0.5*l0)};
    double k2 {h*func1(t + 0.5*h, xi + 0.5*k1, yi + 0.5*l1)};
    double l2 {h*func2(t + 0.5*h, xi + 0.5*k1, yi + 0.5*l1)};
    double k3 {h*func1(t + h, xi + k2, yi + l2)};
    double l3 {h*func2(t + h, xi + k2, yi + l2)};

    xi += (k0 + 2*k1 + 2*k2 + k3)/6;
    yi += (l0 + 2*l1 + 2*l2 + l3)/6;
}


void creating_output_file(vector <double> &time, vector <double> &solution1, vector <double> &solution2, int n){

    ofstream result_file;
    result_file.open("result.txt");
    for(int i = 0; i < n; i++){
        result_file << time.at(i) << ", ";
        result_file << solution1.at(i) << ", ";
        result_file << solution2.at(i) << ",\n";
    }
    result_file.close();

}