#include <bits/stdc++.h>

using namespace std;

const double epsilon = 0.0002;

vector<double> gauss(vector<vector<double>>&, vector<double>&);

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
    double x = u->df(1+epsilon)*v->f(1+epsilon);
    double y = u->f(0+epsilon)*v->f(0+epsilon);
    double z = integrate(0+epsilon,1-epsilon,multiply_fcn(u->df, v->df), 1000);
    double r = 2*integrate(1+epsilon,2-epsilon,multiply_fcn(u->df, v->df), 1000);
    return x - y + z + r;
}

double L(fcn_and_deriv* v){
    return (-1)*20*v->f(0+epsilon);
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
            return ((1.0)*N/L)*(x - a - i*L/N)+1;
        }
        else if(x >= pointB and x <= pointC){
            return ((-1.0)*N/L)*(x - a - i*L/N)+1;
        }
        return 0.0;
    };

    function<double (double)> df = [=](double x){
        if(x < a or x < pointA or x > pointC or i == N){
            return 0.0;
        }
        else if(x >= pointA and x <= pointB){
            return (1.0)*N/L;
        }
        else if(x > pointB and x <= pointC){
            return (-1.0)*N/L;
        }
        return 0.0;
    };

    return new fcn_and_deriv(f, df);
}

vector<double> MES(int n, double a, double b){
    vector<double> coefficients(n);
    vector<fcn_and_deriv*> base_functions;
    for(int i = 0; i <= n; i++){
        base_functions.push_back(create_base_function(a,b, n, i));
    }
    //base_functions.push_back(new fcn_and_deriv([](double x){return 0.0;}, [](double x){return 0.0;}));
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

    // coefficients = gauss(b_matrix, v_matrix);

    
    return gauss(b_matrix, v_matrix);;
}

vector<double> gauss(vector<vector<double>>& B, vector<double>& L){
    B.pop_back();
    int n = B.size();
    for(int i = 0; i < B.size(); i++){
        B[i].pop_back();
        B[i].push_back(L[i]);
    }

    for(int i=0;i<n;i++) 
    {                   
        for(int j=i+1;j<n;j++)
        {
            if(abs(B[i][i]) < abs(B[j][i]))
            {
                for(int k=0;k<n+1;k++)
                {
                    B[i][k]=B[i][k]+B[j][k];
                    B[j][k]=B[i][k]-B[j][k];
                    B[i][k]=B[i][k]-B[j][k];
                }
            }
      }
    }
    for(int i=0;i<n-1;i++)
    {
        for(int j=i+1;j<n;j++)
        {
            double f=B[j][i]/B[i][i];
            for(int k=0;k<n+1;k++)
            {
              B[j][k]=B[j][k]-f*B[i][k];
            }
        }
    }
    vector<double> res(n, 0);

    for(int i = n-1; i>=0; i--)          
    {                     
        res[i]=B[i][n];
                    
        for(int j=i+1;j<n;j++)
        {
          if(i!=j)
          {
              res[i]=res[i]-B[i][j]*res[j];
            }          
        }
        res[i]=res[i]/B[i][i];  
    }

    return res;
}


int main(){
    vector<double> output = MES(12, 0, 2);
    for(int i = 0; i < output.size(); i++){
        cout<< output[i] <<" "<<endl;
    }
    cout<<endl;
}