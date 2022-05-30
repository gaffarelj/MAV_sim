#include <iostream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <vector>

// SRM geometry class
class SRM_geometry
{
public:
    double _R_o;
    double _R_i;
    double _N_f;
    double _w_f;
    double _L_f;
    double _L;
    SRM_geometry(
        double R_o,
        double R_i,
        double N_f,
        double w_f,
        double L_f,
        double L
    )
    {
        this->_R_o = R_o;
        this->_R_i = R_i;
        this->_N_f = N_f;
        this->_w_f = w_f;
        this->_L_f = L_f;
        this->_L = L;
    }

    // Function that returns the burning surface as a function of the burnt distance
    double burningS(double b) {

        if ( _R_i + b >= _R_o ) {
            return 0.0;
        }
        double P_tube = 2*M_PI * (_R_i+b);
        double P_fin = 0.0;
        if (b < _w_f/2) {
            P_fin = 2 * _N_f * _L_f;
        }
        return (P_tube + P_fin) * _L;
    }
};
     
// SRM geometry class
class SRM_thrust_model
{
public:
    SRM_geometry _current_geo = SRM_geometry(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    double _Ve = 0.0;
    double _r = 1e-9;
    double _b = 0.0;
    double _rho_p = 1854.5;
    double _gamma = 1.125;
    double _Ra = 8.314;
    double _M = 0.02414;
    double _eta_Isp = 0.95;
    double _eta_c = 0.93;
    double _eta_F_T = 0.95;
    double _A_t = 0.065;
    double _epsilon = 45;
    double _T_c = 3645;
    double _p_a = 650*0.4;
    double _a = 0.004202/pow(10,6*0.31);
    double _n = 0.31;
    double _p_e = 0.0;
    double _A_e = 0.0;
    double _c_real = 0.0;
    double _Gamma = 0.0;
    double _p_c = 0.0;
    double _last_t = 0.0;
    double _p_ratio = -1.0;

    SRM_thrust_model(
        SRM_geometry current_geo
    )
    {
        this->_current_geo = current_geo;
        _A_e = _A_t * _epsilon;
        _Gamma = sqrt(_gamma) * pow((2/(_gamma+1)),((_gamma+1)/(2*(_gamma-1))));
        double c_ideal = sqrt(_Ra/_M * _T_c) / _Gamma;
        _c_real = c_ideal * _eta_c;
    }

    double a_ratio(double x) {
        return _Gamma / sqrt(2*_gamma/(_gamma-1) *  pow(x/_p_c, 2/_gamma) * (1- pow(x/_p_c, (_gamma-1)/_gamma))) - _A_e/_A_t;
    }

    double bisection(double a, double b) {
        double c = a;

        while ((b-a) >= 1.0e-15) {
            // Find middle point
            c = (a+b)/2;
            // Check if middle point is root
            if (a_ratio(c) == 0.0)
                break;
            // Decide the side to repeat the steps
            else if (a_ratio(c)*a_ratio(a) < 0)
                b = c;
            else
                a = c;
        }
        return c;
    }

    double * magnitude(double t, double y[3]) {
        _b = y[1];

        double S = _current_geo.burningS(_b);
        if (S == 0.0) {
            // Return vector of zeros
            double * zeros = new double[3];
            for (int i = 0; i < 3; i++) {
                zeros[i] = 0.0;
            }
            return zeros;
        }
        double mdot = S * _r * _rho_p;
        _p_c = pow(_c_real * _rho_p * _a * S / _A_t, 1/(1-_n));
        if (_p_ratio == -1.0) {
            _p_e = bisection(1.0, 1.0e5);
            _Ve = sqrt(2*_gamma/(_gamma-1) * _Ra/_M * _T_c * (1-pow((_p_e/_p_c),((_gamma-1)/_gamma))));
            _p_ratio = _p_e / _p_c;
        }
        else {
            _p_e = _p_ratio * _p_c;
        }
        double F_T = mdot * _Ve + (_p_e - _p_a) * _A_e;

        F_T *= _eta_F_T;
        _r = _a * pow(_p_c, _n);  

        double * ret_vec = new double[3];
        ret_vec[0] = -mdot;
        ret_vec[1] = _r;
        ret_vec[2] = F_T;
        return ret_vec;
    }

    std::pair<double*, double*> rk4(double dt, double y[3]) {
        double t = _last_t;
        double * k1 = magnitude(t, y);
        // Multiply all elements of k1 by dt and add to y
        double * y1 = new double[3];
        for (int i = 0; i < 3; i++) { y1[i] = y[i] + k1[i] * dt/2; }
        double * k2 = magnitude(t + dt/2, y1);
        double * y2 = new double[3];
        for (int i = 0; i < 3; i++) { y2[i] = y[i] + k2[i] * dt/2; }
        double * k3 = magnitude(t + dt/2, y2);
        double * y3 = new double[3];
        for (int i = 0; i < 3; i++) { y3[i] = y[i] + k3[i] * dt; }
        double * k4 = magnitude(t + dt, y3);
        double * y_der = new double[3];
        for (int i = 0; i < 3; i++) { y_der[i] = (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6; }
        double * y_next = new double[3];
        for (int i = 0; i < 2; i++) { y_next[i] = y[i] + dt * y_der[i]; }
        y_next[2] = y_der[2];
        _last_t += dt;
        return std::make_pair(y_next, y_der);
    }

    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
    sim_full_burn(double dt) {
        // declare vector of undefined size to save thrust
        std::vector<double> thrust_times;
        std::vector<double> thrust_magnitudes;
        std::vector<double> mass_flows;
        std::pair <double*, double*> integ_res;
        double * y = new double[3];
        double * y_der;
        y[0] = 1e3;
        int zero_thrust = 0;
        while ((zero_thrust != 2)) {
            integ_res = rk4(dt, y);
            y = integ_res.first;
            y_der = integ_res.second;
            // Print the time, thrust magnitude, mass flow rate
            // std::cout << std::setprecision(7) << _last_t
            //    << "~" << std::setprecision(7) << y_der[2]
            //    << "~" << std::setprecision(7) << y_der[0] << std::endl;
            if (y_der[2] == 0.0) {
                zero_thrust++;
            }
            else {
                thrust_times.push_back(_last_t);
                thrust_magnitudes.push_back(y_der[2]);
                mass_flows.push_back(y_der[0]);
            }
               
        }
        return std::make_tuple(thrust_times, thrust_magnitudes, mass_flows);
    }

};

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> 
run_sim(double Ro, double Ri, double Nf, double wf, double Lf, double L, double dt) {
    SRM_geometry SRM = SRM_geometry(Ro, Ri, Nf, wf, Lf, L);
    SRM_thrust_model SRM_thrust = SRM_thrust_model(SRM);
    return SRM_thrust.sim_full_burn(dt);
}

// int main(int argc, char *argv[]) {
//     // Arguments: Ro, Ri, Nf, wf, Lf, L, dt
//     double Ro = std::stod(argv[1]);
//     double Ri = std::stod(argv[2]);
//     double Nf = std::stod(argv[3]);
//     double wf = std::stod(argv[4]);
//     double Lf = std::stod(argv[5]);
//     double L = std::stod(argv[6]);
//     double dt = std::stod(argv[7]);
//     // std::cout << argc << " arguments" << std::endl;
//     // Extract value of arguments
//     // double b = std::stod(argv[1]);
//     SRM_geometry SRM = SRM_geometry(Ro, Ri, Nf, wf, Lf, L);
//     SRM_thrust_model SRM_thrust = SRM_thrust_model(SRM);
//     SRM_thrust.sim_full_burn(dt);
//   return 0;
// }

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(SRM_cpp, m) {
    m.def("run_sim", &run_sim,
        py::arg("Ro"),
        py::arg("Ri"),
        py::arg("Nf"),
        py::arg("wf"),
        py::arg("Lf"),
        py::arg("L"),
        py::arg("dt")
    );
}

// To compile, run:
// g++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) multi_fin.cpp -o SRM_cpp.so