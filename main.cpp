/// ---------------------------------------------------
/// Class for fast tridiagonal matrix solve and inverse
/// J.-P. Champeaux (2024)
/// a[i] : matrix diagonal  [N]
/// b[i] : matrix upper diagonal [N-1]
/// c[i] : matrix lower diagonal [N-1]
/// Solving system TDmatrix * X = Y
/// [ref] Thomas algorithm
/// Inversing TDmatrix :
/// [ref] R. Usmani, Inversion of a tridiagonal Jacobi matrix, Linear Algebra Appl. 212/213 (1994) 413-414
///
/// -------------------------------------------------

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>

#include "TDMatrix.h"

using namespace std;


int main()
{
    clock_t t1, t2;
    t1=clock();
    std::cout<<"-----------------------------------------------\n";
    std::cout<<"    TDMatrix Class RELEASE 02/02/2024 \n";
    std::cout<<"    Coded by JP CHAMPEAUX  (WORKING VERSION)   \n";
    std::cout<<"-----------------------------------------------\n";


    /// The tridiagonal Matrix
    vector<double> diag = {5.0, 6.0, 7.0, 8.0};
    vector<double> upperdiag = {9.0, 10.0, 11.0};
    vector<double> lowerdiag = {1.0, 2.0, 3.0};

    /// Call of class TDMatrix
    TDMatrix<double> M(diag, upperdiag,lowerdiag);

    cout<<"The tridiagonal matrix :\n"<< M <<endl;

    /// Evaluate determinant
    M.Eval_Teta();
    cout<<"Determinant = "<<M.get_det()<<endl;

    /// Evaluate Invert Matrix
    vector<vector<double>> Ainv;
    Ainv = M.InvertTridiagonal();

    /// Display the inverted matrix
    std::cout << "\nInverted Matrix:\n";
    for (int i = 0; i < M.get_size(); ++i)
    {
        for (int j = 0; j < M.get_size(); ++j)
        {
            std::cout<<std::scientific << Ainv[i][j] << "\t";
        }
        std::cout << "\n";
    }
    Ainv.clear();

    /// Solve linear system
    vector<double> Y = {1,2,3,4};
    vector<double> X;

    cout<<"\nSolving Mx=Y :"<<endl;
    X = M.solve_Ax_(Y);
    cout<<"X ={ ";
    for (auto x:X) cout<<x<<" , ";

    t2=clock();
    float temp=(float) (t2-t1)/CLOCKS_PER_SEC;
    cout<<" }"<<endl;
    cout<<"\nAll done in "<<temp<<" sec"<<endl<<endl;
    cout<<"Have Fun World is a Playground"<<endl;
    cout<<endl;
    cout<<" (°)~(°) \n";
    cout<<" (=0.0=) \n";
    cout<<"   (w)  \n";
    cout<<" (:) (:)___°\n";
    return 0;


}
