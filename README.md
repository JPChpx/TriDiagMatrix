# ⚡ Fast Tridiagonal Matrix Solver and Inverter
![C++](https://img.shields.io/badge/C%2B%2B-17-blue) ![Status](https://img.shields.io/badge/status-stable-green) ![Version](https://img.shields.io/badge/version-1.0.0-blueviolet)


>Une classe C++ performante pour résoudre et inverser des matrices tridiagonales, basée sur l’algorithme de **Thomas** pour la résolution et la méthode de **Usmani** pour l’inversion.
>standalone easy to use C++ Class for fast Tridigonal matrix inversion and A*X=Y resolution using Thomas algoritm

> 📌 Auteur : J.-P. Champeaux  
> 📅 Date : Février 2024  

 ![alt text](https://external-content.duckduckgo.com/iu/?u=https%3A%2F%2Fdocs.oracle.com%2Fcd%2FE77782_01%2Fhtml%2FE77802%2Ffigures%2Fequation1211.png&f=1&nofb=1&ipt=2fc7441ade3d5bb0c9e3a5a167cab031d5bd5ff728053c6e0f38495a5ba08424&ipo=images)

---

## 🚀 Fonctionnalités

- ⚡ **Solveur rapide** pour systèmes tridiagonaux `Ax = b` (algorithme de Thomas)
- 🧠 **Calcul du déterminant** de la matrice
- 🔁 **Inversion complète** de la matrice tridiagonale (`A⁻¹`) via la méthode analytique de Usmani
- 🧾 Opérateur de sortie `<<` pour afficher joliment la matrice complète

---

## 📚 Références

- **Algorithme de Thomas** (méthode directe pour résoudre un système tridiagonal)
- **Usmani, R. A.**, *Inversion of a tridiagonal Jacobi matrix*, Linear Algebra and its Applications, Vol. 212/213, 1994, pp. 413–414

---
 
### 📦 Compilation

```bash
g++ -std=c++17 -O3 -o tdsolver main.cpp
```

---
## 🔧 Utilisation

```cpp
#include "TDMatrix.h"  // Inclure la classe optimisée

int main() {
    std::vector<double> diag = {4, 4, 4};
    std::vector<double> upper = {1, 1};
    std::vector<double> lower = {1, 1};

    TDMatrix<double> td(diag, upper, lower);

    std::vector<double> Y = {7, 8, 7};
    std::vector<double> X = td.solve_Ax_(Y);

    std::cout << "Solution X:\n";
    for (auto val : X) std::cout << val << " ";
    std::cout << "\n";

    std::cout << "Matrice inverse A⁻¹ :\n";
    auto Ainv = td.InvertTridiagonal();
    for (const auto& row : Ainv) {
        for (auto val : row) std::cout << val << "\t";
        std::cout << "\n";
    }

    return 0;
}
```
