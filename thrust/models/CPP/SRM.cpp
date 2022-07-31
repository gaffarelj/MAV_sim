#include <iostream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <chrono>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// SRM geometry class
class SRM_geometry
{
public:

    std::string _type;
    double _R_o;
    double _R_i;
    double _R_mid;
    double _L;
    int _N_f;
    double _w_f;
    double _L_f;
    double _w;
    double _r_f;
    double _delta_s;
    int _N_a;
    double _w_t;
    double _w_ds;
    double _R_e;

    SRM_geometry(
        std::string type = "",
        double R_o = 0.0,
        double R_i = 0.0,
        double R_mid = 0.0,
        double L = 0.0,
        int N_f = 0,
        double w_f = 0.0,
        double L_f = 0.0,
        double w = 0.0,
        double r_f = 0.0,
        double delta_s = 0.0,
        int N_a = 0
    )
    {
        _type = type;
        _R_o = R_o;
        _R_i = R_i;
        _R_mid = R_mid;
        _L = L;
        _N_f = N_f;
        _w_f = w_f;
        _L_f = L_f;
        _w = w;
        _r_f = r_f;
        _delta_s = delta_s;
        _N_a = N_a;
        _w_t = sqrt((_w + _r_f)*(_w + _r_f) + 2 * _R_o*_R_o - 2 * _R_o * sqrt(_R_o*_R_o - 2 * _R_o * (_w + _r_f)) + _w + _r_f) - _r_f;
        _w_ds = (_R_i*_R_i + ( 2 * _w + _r_f - sqrt(_R_i*_R_i + 2 * _R_i * (2 * _w + _r_f) + 3 * _w*_w + 2 * _r_f * _w) ) * _R_i + 2 * _w * (_w + _r_f)) / (sqrt(_R_i*_R_i + 2 * _R_i * (2 * _w + _r_f) + 3 * _w*_w + 2 * _r_f * _w) + _r_f - _R_i);
        _R_e = _R_i/2;
    };

    // Function that returns the burning surface as a function of the burnt distance
    double burningS(double b) {
        if (_type == "tubular")
        {
            if ( _R_i + b >= _R_o ) {
                return 0.0;
            }
            double P_tube = 2*M_PI * (_R_i+b);
            return P_tube * _L;
        }
        else if (_type == "rod and tube")
        {
            // std::cout << "R_o " << _R_o << std::endl;
            // std::cout << "R_i " << _R_i << std::endl;
            // std::cout << "R_mid " << _R_mid << std::endl;
            // std::cout << "L " << _L << std::endl;
            // std::cout << "b " << b << std::endl;
            double P_tube = 0.0;
            double P_rod = 0.0;
            if (_R_mid + b < _R_o) {
                P_tube = 2*M_PI * (_R_mid+b);
            }
            if (_R_i - b > 0) {
                P_rod = 2*M_PI * (_R_i-b);
            }
            return (P_tube + P_rod) * _L;
        }
        else if (_type == "multi fin")
        {
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
        else if (_type == "anchor")
        {
            double P1 = 0.0, P2 = 0.0, P3 = 0.0, P4 = 0.0, P5 = 0.0, P6 = 0.0, P7 = 0.0;
            // Perimeter 1
            if (0 <= b && b <= _w)
            {
                P1 = (_R_i + b) * (M_PI/_N_a - asin((_delta_s+2*b)/(2*(_R_i+b))));
            }
            else if (_w < b && b <= _w_ds)
            {
                P1 = (_R_i+b) * ( asin( (_r_f+_w)/(_R_i+2*_w+_r_f) ) -
                    acos( (_R_i*_R_i+(2*_w+b+_r_f)*_R_i+(2*_w-b)*_r_f+2*_w*_w)/((_R_i+2*_w+_r_f)*(_R_i+b)) ) );
            }
            // Perimeter 2
            if (b <= _w)
            {
                P2 = sqrt(pow(_R_i+2*_w-b,2)-pow(_delta_s/2+b,2)) - sqrt(pow(_R_i+b,2)-pow(_delta_s/2+b,2));
            }
            // Perimeter 3
            if (b <= _w)
            {
                P3 = (_R_i+2*_w-b)*( M_PI/_N_a - asin( (_delta_s+2*b)/(2*(_R_i+2*_w-b)) ) -
                    asin( (_r_f+_w)/(_R_i+2*_w+_r_f) ) );
            }
            // Perimeter 4
            if (b <= _w)
            {
                P4 = (_r_f+b) * acos( (_r_f+_w)/(_R_i+2*_w+_r_f) );
            }
            else if (b <= _w_ds)
            {
                P4 = (_r_f+b) * (acos( (_r_f+_w)/(_R_i+2*_w+_r_f) ) -
                    acos( (_r_f+_w)/(_r_f+b) ) - acos( (pow(_R_i+2*_w+_r_f,2)+pow(_r_f+b,2)-pow(_R_i+b,2))/(2*(_R_i+2*_w+_r_f)*(_r_f+b)) ) );
            }
            // Perimeter 5
            if (b <= _w)
            {
                P5 = sqrt(pow(_R_o-_w-_r_f,2)-pow(_r_f+_w,2)) - sqrt(pow(_R_i+2*_w+_r_f,2)-pow(_r_f+_w,2));
            }
            // Perimeter 6
            if (b <= _w)
            {
                P6 = (_r_f+b)* (M_PI/2 + asin( (_r_f+_w)/(_R_o-_w-_r_f) ));
            }
            else if (b <= _w_t)
            {
                P6 = (_r_f+b)* (acos( (pow(_R_o-_w-_r_f,2)+pow(_r_f+b,2)-pow(_R_o,2))/(2*(_R_o-_w-_r_f)*(_r_f+b)) ) -
                    acos( (_r_f+_w)/(_r_f+b) ) - acos( (_r_f+_w) / (_R_o-_w-_r_f) ) );
                if (P6 < 0.0)
                {
                    P6 = 0.0;
                }
            }
            // Perimeter 7
            if (b <= _w)
            {
                P7 = (_R_o-_w+b) * (M_PI/_N_a - asin( (_r_f+_w)/(_R_o-_w-_r_f) ));
            }
            double S = (P1 + P2 + P3 + P4 + P5 + P6 + P7) * 2 * _N_a * _L;
            return S;
        }
        else if (_type == "spherical")
        {
            if ( _R_i + b >= _R_o ) {
                return 0.0;
            }
            double S = 4 * M_PI * pow(_R_i+b, 2) - M_PI * pow(_R_e+b, 2) + 2*M_PI * (_R_e+b) * (_R_o - _R_i - b);
            return S;
        }
        else
        {
            throw std::invalid_argument("Invalid type of SRM geometry");
        }
    }
};
     
class SRM_thrust_model
{
public:
    SRM_geometry _current_geo;

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
        _current_geo = current_geo;
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
            double p_e_start = pow(10, 7);
            _p_e = p_e_start;
            // Try with different initial guess for exhaust pressure until value found (or stop when P_e below 1 Pa)
            while (_p_e == p_e_start) {
                p_e_start = _p_e/10;
                _p_e = bisection(1.0, p_e_start);
                if (_p_e < 1.0) {
                    throw std::invalid_argument("Could not find p_e");
                }
            }
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
        delete(k1);
        delete(k2);
        delete(k3);
        delete(k4);
        delete(y1);
        delete(y2);
        delete(y3);
        return std::make_pair(y_next, y_der);
    }

    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
    sim_full_burn(double dt, double timeout) {
        // declare vector of undefined size to save thrust
        std::vector<double> thrust_times;
        std::vector<double> thrust_magnitudes;
        std::vector<double> mass_flows;
        // Measure the current system time
        std::chrono::steady_clock::time_point clock_begin = std::chrono::steady_clock::now();
        double * y = new double[3];
        y[0] = 1e3;
        int zero_thrust_count = 0;
        double last_saved_dt = 0.0;
        while ((zero_thrust_count != 2)) {
            if (timeout != -1.0) {
                // Measure the time elapsed in seconds
                std::chrono::steady_clock::time_point clock_end = std::chrono::steady_clock::now();
                double dt_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(clock_end - clock_begin).count() / 1.0e3;
                // Check if timeout has been reached
                if (dt_elapsed > timeout) {
                    std::cout << "*** Timeout reached, stopping burn simulation ***" << std::endl;
                    std::vector<double> neg_vec(1, -1.0);
                    return std::make_tuple(neg_vec, neg_vec, neg_vec);
                }
            }
            std::pair <double*, double*> integ_res = rk4(dt, y);
            delete(y);
            y = integ_res.first;
            double * y_der = integ_res.second;
            last_saved_dt += dt;
            // Print the time, thrust magnitude, mass flow rate
            // std::cout << std::setprecision(7) << _last_t
            //    << "~" << std::setprecision(7) << y_der[2]
            //    << "~" << std::setprecision(7) << y_der[0] << std::endl;
            if (y_der[2] == 0.0) {
                zero_thrust_count++;
            }
            else if (last_saved_dt >= 0.0001) {
                last_saved_dt = 0.0;
                thrust_times.push_back(_last_t);
                thrust_magnitudes.push_back(y_der[2]);
                mass_flows.push_back(y_der[0]);
            }
            delete(y_der);
        }
        return std::make_tuple(thrust_times, thrust_magnitudes, mass_flows);
    }

};


std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> 
run_sim(SRM_geometry SRM, double dt, double timeout=-1.0) {
    SRM_thrust_model SRM_thrust = SRM_thrust_model(SRM);
    return SRM_thrust.sim_full_burn(dt, timeout);
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

namespace py = pybind11;

PYBIND11_MODULE(SRM_cpp, m) {

    m.def("run_sim",
        &run_sim,
        py::arg("SRM_model"),
        py::arg("dt"),
        py::arg("timeout") = -1.0
    );

    py::class_<SRM_geometry>(m, "SRM_geometry")
        .def(
            py::init<
                const std::string &,
                const float,
                const float,
                const float,
                const float,
                const int,
                const float,
                const float,
                const float,
                const float,
                const float,
                const int
            >(),
            py::arg("type"),
            py::arg("R_o") = 0,
            py::arg("R_i") = 0,
            py::arg("R_mid") = 0,
            py::arg("L") = 0,
            py::arg("N_f") = 0,
            py::arg("w_f") = 0,
            py::arg("L_f") = 0,
            py::arg("w") = 0,
            py::arg("r_f") = 0,
            py::arg("delta_s") = 0,
            py::arg("N_a") = 0
    );
}

// To compile, run:
// g++ -O3 -Wall -shared -std=c++14 -fPIC $(python3 -m pybind11 --includes) SRM.cpp -o SRM_cpp.so