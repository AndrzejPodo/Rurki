#include <bits/stdc++.h>

using namespace std;

struct fcn_and_deriv{
    function<double (double)> f;
    function<double (double)> df;

    fcn_and_deriv(function<double (double)> f, function<double (double)> df){
        this->f = f;
        this->df = df;
    }
};

function<double (double)> multiply_fcn(function<double (double)> f1, function<double (double)> f2){
    return [f1, f2](double x){ return f1(x)*f2(x);};
}

double integrate(double a, double b, function<double (double)> func, int n){
    double result = 0;
    double delta = (b-a)/n;
    double x1, x2, y1, y2;
    for(int i = 0; i <= n-1; i++){
        x1 = a + i*delta;
        x2 = a + (i+1)*(delta);
        y1 = func(x1);
        y2 = func(x2);
        result += delta*(y1+y2)/2;
    } 
    return result; 
}

double B(fcn_and_deriv* u, fcn_and_deriv* v){
    return u->df(1)*v->f(1) - u->f(0)*v->f(0) + integrate(0,1,multiply_fcn(u->df, v->df), 1000) + 2*integrate(1,2,multiply_fcn(u->df, v->df), 1000);
}

double L(fcn_and_deriv* v){
    return (-1)*20*v->f(0);
}

fcn_and_deriv* create_base_function(double a, double b, int N, int i){
    double L = b - a;
    double pointA = (-1)*L/N + a + i*L/N;
    double pointB = a+i*L/N;
    double pointC = L/N + a + i*L/N;
    function<double (double)> f = [=](double x){
        if(x < a or x < pointA or x > pointC or i == N){
            return 0.0;
        }
        else if(x >= pointA and x <= pointB){
            return N/L * (x - a - i*L/N)+1;
        }
        else{
            return (-1)*N/L * (x - a - i*L/N)+1;
        }
    };

    function<double (double)> df = [=](double x){
        if(x < a or x < pointA or x > pointC or i == N){
            return 0.0;
        }
        else if(x >= pointA and x <= pointB){
            return N/L;
        }
        else{
            return (-1)*N/L;
        }
    };

    return new fcn_and_deriv(f, df);
}

vector<double> MES(int n, double a, double b){
    vector<double> coefficients;
    vector<fcn_and_deriv*> base_functions;
    for(int i = 0; i <= n; i++){
        base_functions.push_back(create_base_function(a,b, n, i));
    }
    vector<vector<double>> b_matrix;
    for(int i = 0; i <= n; i++){
        vector<double> temp(n+1);
        b_matrix.push_back(temp);
    }
    vector<double> v_matrix(n+1);

    for(int i = 0; i <= n; i++){
        v_matrix[i] = L(base_functions[i]);
        for(int j = 0; j <= n; j++){
            b_matrix[i][j] = B(base_functions[i], base_functions[j]);
        }
    }

    for(int i = 0; i <= n; i++){
        for(int j = 0; j <= n; j++){
            cout<<b_matrix[i][j]<<" ";
        }
        cout<<endl;
    }

    cout<<endl;
    for(int i = 0; i <= n; i++){
        cout<<v_matrix[i]<<" ";
    }
    cout<<endl;

    return {-1,-1};
}

int main(){

    MES(5, 0, 2);


}