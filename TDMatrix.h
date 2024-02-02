/// ---------------------------------------------------
/// Class for fast tridiagonal matrix solve and inverse
/// J.-P. Champeaux (02/2024)
/// a[i] : matrix diagonal  [N]
/// b[i] : matrix upper diagonal [N-1]
/// c[i] : matrix lower diagonal [N-1]
/// Solving system TDmatrix * X = Y
/// [ref] Thomas algorithm
/// Inversing TDmatrix :
/// [ref] R. Usmani, Inversion of a tridiagonal Jacobi matrix, Linear Algebra Appl. 212/213 (1994) 413-414
///
/*
(°)~(°)
(=0.0=)
  (w)
(:) (:)___°
*/
/// -------------------------------------------------


#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

template <typename T>
class TDMatrix {
private :
    std::vector<T> a; // diag
    std::vector<T> b; // upper diag
    std::vector<T> c; // lower diag
    std::vector<T> theta;
    std::vector<T> phi;
    size_t N;
    T det;

public :
    // Constructor
    TDMatrix(const std::vector<T>& diag, const std::vector<T>& upper, const std::vector<T>& lower)
        : a(diag), b(upper), c(lower) { N = a.size(); }

// ----------------------------------------
// Tridiag Det and theta vector calulation
// ----------------------------------------
void Eval_Teta()
{
    theta.resize(N);
    theta[0] = a[0];
    theta[1] = a[1] * theta[0] - b[0] * c[0];
    for (int i = 2; i < N; i++)
    {
        theta[i] = a[i] * theta[i - 1] - b[i - 1] * c[i - 1] * theta[i - 2];
    }
    det = theta[N-1];
}

// ------------------------------
// Tridiag Phi vector calculation
// ------------------------------

void Eval_Phi()
{
    phi.resize(N);
    phi[N - 1] = a[N - 1];
    phi[N - 2] = a[N - 2] * phi[N - 1] - b[N - 2] * c[N - 2];
    for (int i = N - 3; i >= 0; i--)
    {
        phi[i] = a[i] * phi[i + 1] - b[i] * c[i] * phi[i + 2];
    }

}

T get_det() { return det; }
size_t get_size() { return N; }


// -----------------------------
// Tridiagonal Thomas Solver...
// -----------------------------
vector<T> solve_Ax_(vector<T> Y) {
    int n = a.size()-1;
    vector<double> lower_copy = this->c;
    vector<double> Y_copy = Y;
    lower_copy[0] /= a[0];
    Y_copy[0] /= a[0];
    for (int i = 1; i < n; i++) {
        lower_copy[i] /= a[i] - b[i] * lower_copy[i - 1];
        Y_copy[i] = (Y_copy[i] - b[i] * Y_copy[i - 1]) / (a[i] - b[i] * lower_copy[i - 1]);
    }
    Y_copy[n] = (Y_copy[n] - b[n] * Y_copy[n - 1]) / (a[n] - b[n] * lower_copy[n - 1]);
    for (int i = n; i-- > 0;) {
        Y_copy[i] -= lower_copy[i] * Y_copy[i + 1];
    }
    return Y_copy;
}

// -----------------------------
// Tridiagonal Matrix invesion
// -----------------------------

std::vector<std::vector<T>> InvertTridiagonal( )
{

    std::vector<std::vector<T>> Ainv(N, std::vector<T>(N));
    this->Eval_Teta();
    this->Eval_Phi();

    // -----------------------------------------------------------------------
    // T-1 (i,j) =
    // if(i<j) -1^{i+j} * b[i] *...* b[j-1] * teta[i-1] * Phi[j+1] / teta[n-1]
    // if(i=j) teta[i-1] * phi(j+1) /teta(n-1)
    // if(i>j) -1^{i+j} * c[j] *...* c[i-1] * teta[j-1] * Phi[i+1] / teta[n-1]
    // -----------------------------------------------------------------------
    for (int i = 0; i < N; i++)
    {
        double th_im1 = (i == 0) ? 1.0 : theta[i - 1];
        double ph_ip1 = (i == N - 1) ? 1.0 : phi[i + 1];
        for (int j = 0; j < N; j++)
        {
            double th_jm1 = (j == 0) ? 1.0 : theta[j - 1];
            double ph_jp1 = (j == N - 1) ? 1.0 : phi[j + 1];

            if (i < j)
            {
                double bprod = 1.0;
                for (int k = i; k <= j - 1; k++)
                {
                    bprod *= b[k];
                }
                Ainv[i][j] = ((i + j) % 2 == 0 ? 1.0 : -1.0) * bprod * th_im1 * ph_jp1 / theta[N - 1];
            }
            else if (i == j)
            {
                Ainv[i][j] = th_im1 * ph_jp1 / theta[N - 1];
            }
            else
            {
                double cprod = 1.0;
                for (int k = j; k <= i - 1; k++)
                {
                    cprod *= c[k];
                }
                Ainv[i][j] = ((i + j) % 2 == 0 ? 1.0 : -1.0) * cprod * th_jm1 * ph_ip1 / theta[N - 1];
            }
        }
    }
    return Ainv;
}
// ---------------
// print M
// ---------------
 // Opérateur de sortie <<
    friend std::ostream& operator<<(std::ostream& os, const TDMatrix& matrix) {
        int size = matrix.a.size();

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (i == j) {
                    os << matrix.a[i] << "\t";
                } else if (i == j - 1) {
                    os << matrix.b[j - 1] << "\t";
                } else if (i - 1 == j) {
                    os << matrix.c[i - 1] << "\t";
                } else {
                    os << "0\t";
                }
            }
            os << "\n";
        }

        return os;
    }

};
