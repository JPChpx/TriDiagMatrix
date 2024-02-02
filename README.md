# TriDiagMatrix
standalone easy to use C++ Class for fast Tridigonal matrix inversion and A*X=Y resolution using Thomas algoritm

![alt text](https://external-content.duckduckgo.com/iu/?u=https%3A%2F%2Fdocs.oracle.com%2Fcd%2FE77782_01%2Fhtml%2FE77802%2Ffigures%2Fequation1211.png&f=1&nofb=1&ipt=2fc7441ade3d5bb0c9e3a5a167cab031d5bd5ff728053c6e0f38495a5ba08424&ipo=images)

 ---------------------------------------------------
 coded by J.-P. Champeaux (2024)

 a[i] : matrix diagonal  [N]
 
 b[i] : matrix upper diagonal [N-1]
 
 c[i] : matrix lower diagonal [N-1]
 
 - Solving system TDmatrix * X = Y
 
 [ref] Thomas algorithm
 
 - Inversing TDmatrix :
 
 [ref] R. Usmani, Inversion of a tridiagonal Jacobi matrix, Linear Algebra Appl. 212/213 (1994) 413-414
 
 -------------------------------------------------

 Exemple of Use : (dont forget to include "TDMatrix.h") 
      
     #include <iostream>
#include <vector>
#include <cmath>
#include "TDMatix.h"

using namespace std;


    int main()
    {

    /// The tridiagonal Matrix
    vector<double> diag = {5.0, 6.0, 7.0, 8.0};
    vector<double> upperdiag = {9.0, 10.0, 11.0};
    vector<double> lowerdiag = {1.0, 2.0, 3.0};

    /// Call of class TDMatrix
    TDMatrix<double> M(diag, upperdiag,lowerdiag);
    cout<<M;
    
    /// Evaluate determinant
    M.Eval_Teta();
    cout<<"det = "<<M.get_det()<<endl;

    /// Evaluate M^-1
    vector<vector<double>> Ainv;
    Ainv = M.InvertTridiagonal();

    /// Display the inverted matrix
    std::cout << "Inverted Matrix:\n";
    for (int i = 0; i < M.get_size(); ++i)
    {
        for (int j = 0; j < M.get_size(); ++j)
        {
            std::cout << Ainv[i][j] << "\t";
        }
        std::cout << "\n";
    }
    Ainv.clear();

    /// Solve linear system
    vector<double> Y = {1,2,3,4};
    vector<double> X;
    X = M.solve_Ax_(Y);
    cout<<"X ={ ";
    for (auto x:X) cout<<x<<" , ";
    cout<<" }"<<endl;
    return 0;
    
    }

    (°)~(°)  
    (=0.0=)  
      (w)   
    (:) (:)___°  
