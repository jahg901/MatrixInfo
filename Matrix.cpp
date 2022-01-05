#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
using namespace std;

const int MAX_ARRAY_SIZE = 10;
const int ASCII_ZERO = 48;

int gcd(int a, int b) {
  if (b == 0)
    return a;
  return gcd(b, a % b);
}

string toString(int x) {
  if (x == 0) return "0";
  
  string sign = "";
  if (x < 0) {
    x = -x;
    sign = "-";
  }
  
  string s = "";
  
  int digit = 1;
  while (digit <= x) {
    s = char((x % (digit * 10)) / digit + ASCII_ZERO) + s;
    digit *= 10;
  }
  
  s = sign + s;
  return s;
}

class Fraction {
  private:
    int num;
    int denom;
    
  public:
    Fraction(int numerator = 0, int denominator = 1) {
      int divisor = gcd(numerator, denominator);
      num = numerator / divisor;
      denom = denominator / divisor;
    }
    
    Fraction operator - () {
      Fraction negative(-num, denom);
      return negative;
    }
    
    Fraction operator + (const Fraction& f) {
      Fraction sum(num * f.denom + denom * f.num, denom * f.denom);
      return sum;
    }
    
    Fraction operator - (const Fraction& f) {
      Fraction diff(num * f.denom - denom * f.num, denom * f.denom);
      return diff;
    }
    
    Fraction operator * (const Fraction& f) {
      Fraction prod(num * f.num, denom * f.denom);
      return prod;
    }
    
    Fraction operator / (const Fraction& f) {
      Fraction quot(num * f.denom, denom * f.num);
      return quot;
    }
    
    operator string() {
      if (denom < 0) {
        num = -num;
        denom = -denom;
      }
      if (denom == 0)
        return "?";
        
      string s = toString(num);
      if (denom != 1)
        s += "/" + toString(denom);
        
      return s;
    }
    
    operator double() {
      return 1.0 * num / denom;
    }
    
    Fraction recip() {
      Fraction rcp(denom, num);
      return rcp;
    }
};

ostream& operator << (ostream& out, Fraction f) { return out << string(f); }

class Matrix {
  public:
    Fraction arr[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE];
    int numRows, numCols;
    
    Matrix(int r = 0, int c = 0) {
      numRows = r;
      numCols = c;
    }
    
    void print() {
      string sMatrix[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE];
      int maxRowSizes[MAX_ARRAY_SIZE] = {0};
      
      for (int row = 0; row < numRows; row++)
        for (int col = 0; col < numCols; col++) {
          sMatrix[row][col] = string(arr[row][col]);
          
          if (sMatrix[row][col].size() > maxRowSizes[col])
            maxRowSizes[col] = sMatrix[row][col].size();
        }
      
      cout << endl;
       for (int row = 0; row < numRows; row++) {
        for (int col = 0; col < numCols; col++)
          cout << setw(maxRowSizes[col] + 1 - (col == 0)) << sMatrix[row][col];
        cout << endl;
      }
    }
    
    Matrix swapRows(int rowA, int rowB) {
      rowA--;
      rowB--;
      Fraction temp;
      for (int col = 0; col < numCols; col++) {
        temp = arr[rowA][col];
        arr[rowA][col] = arr[rowB][col];
        arr[rowB][col] = temp;
      }
      return *this;
    }
    
    Matrix multiply(int row, Fraction scalar) {
      row--;
      for (int col = 0; col < numCols; col++) {
        arr[row][col] = arr[row][col] * scalar;
      }
      return *this;
    }
    
    Matrix addScalarMultiple(int rowA, int rowB, Fraction scalar) {
      rowA--;
      rowB--;
      for (int col = 0; col < numCols; col++) {
        arr[rowA][col] = arr[rowA][col] + (arr[rowB][col] * scalar);
      }
      return *this;
    }
    
    Matrix RREF() {
      Matrix M = *this;
      int freeVars = 0;
      for (int i = 1; i <= M.numRows && i + freeVars <= M.numCols; i++) {
        while (i + freeVars <= M.numCols
               && double(M.arr[i - 1][i - 1 + freeVars]) == 0) {
          for (int k = i; k < numRows; k++) {
            if (double(M.arr[k][i - 1 + freeVars]) != 0) {
              M.addScalarMultiple(i, k + 1, 1);
              break;
            }
          }
          if (double(M.arr[i - 1][i - 1 + freeVars]) == 0) freeVars++;
        }
        if (i + freeVars > M.numCols) break;
        
        if (double(M.arr[i - 1][i - 1 + freeVars]) != 1)
          M.multiply(i, M.arr[i - 1][i - 1 + freeVars].recip());
        for (int j = 1; j <= M.numRows; j++) {
          if (j != i)
            M.addScalarMultiple(j, i, -M.arr[j - 1][i - 1 + freeVars]
                                    / M.arr[i - 1][i - 1 + freeVars]);
        }
      }
      return M;
    }
    
    Matrix exclude(int row, int col)
    {
      if (row >= 0 && row < numRows && col >= 0 && col < numCols)
      {
        Matrix m(numRows - 1, numCols - 1);
        for (int rowIndex = 0; rowIndex < numRows; rowIndex++)
          for (int colIndex = 0; colIndex < numCols; colIndex++)
            if (rowIndex != row && colIndex != col)
            {
              int realRowIndex = rowIndex, realColIndex = colIndex;
              if (rowIndex > row)
                realRowIndex--;
                
              if (colIndex > col)
                realColIndex--;
              
              m.arr[realRowIndex][realColIndex] = arr[rowIndex][colIndex];
            }
        return m;
      }
      return *this;
    }
    
    Fraction cofactor(int row, int col)
    {
      //cout << "Exclude " << row << " " << col << endl;
      return pow(-1, row + col) * exclude(row, col).determinant();
    }
    
    Fraction determinant()
    {
      //cout << "Size: " << numRows << endl;
      if (numRows != numCols)
      {
        return -1000000;
      }
      
      if (numRows == 1)
        return arr[0][0];
      
      if (numRows == 2)
        return (arr[0][0] * arr[1][1]) - (arr[0][1] * arr[1][0]);
      
      Fraction d;
      for (int i = 0; i < numRows; i++)
      {
        d = d + arr[0][i] * cofactor(0, i);
        //cout << d << endl;
      }
      
      return d;
    }
};

int main()
{
  int temp = 0;
  Matrix M;
  
  cout << "Enter the number of rows in your matrix: ";
  cin >> M.numRows;
  
  cout << "Enter the number of columns in your matrix: ";
  cin >> M.numCols;
  
  cout << "Enter the matrix:" << endl;
  for (int row = 0; row < M.numRows; row++)
    for (int col = 0; col < M.numCols; col++) {
      cin >> temp;
      M.arr[row][col] = temp;
  }
  
  M.print();
  M.RREF().print();
  
  cout << endl << "Determinant: " << M.determinant();
  return 0;
}
