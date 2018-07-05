
module matrix_dynamic;

import std.stdio;
import std.math;
import std.conv;
import std.random;
import std.datetime;

T[][] dup2D(T)(const ref T[][] arr) {
    T[][] copy;
    copy.length = arr.length;
    for (size_t r = 0; r < arr.length; r++) {
        copy[r].length = arr[r].length;
        copy[r] = arr[r].dup;
    }
    return copy;
}

class Matrix(T) {

public:

    this() {
        mat.length = 0;
    }

    this(size_t rows, size_t cols, T initialValue = T(0)) {
        if (rows > 0) {
            mat.length = rows;
            foreach (ref T[] c; mat) {
                c.length = cols;
                c[0..cols] = initialValue;
            }
        }
    }

    this(const ref Matrix!(T) rhs) {
        mat = dup2D!(T)(rhs.mat);
    }

    // for some reason the .dup of a 2D array is a reference and not a copy
    // so I can't just do mat = m.dup
    this(const ref T[][] m) {
        mat = dup2D!(T)(m);
    }

    static Matrix!(T) identity(size_t size) {
        Matrix!(T) ident = new Matrix!(T)(size, size);
        for (size_t i = 0; i < size; i++) {
            ident[i, i] = 1;
        }
        return ident;
    }

    static Matrix!(T) random(size_t size) {
        Matrix!(T) random = new Matrix!(T)(size, size);
        foreach (ref T[] r; random.mat) {
            foreach (ref T v; r) {
                v = cast(T)(uniform(0, 100) / 100.0f);
            }
        }
        return random;
    }

    // override Matrix[r, c] op= value
    // must be defined before Matrix[r, c] operator
    T opIndexOpAssign(string op)(T value, size_t r, size_t c)
    in {
        assert(r < mat.length && r >= 0);
        if (mat.length > 0) {
            assert(c < mat[0].length && c >= 0);
        }
    }
    do {
        static if (op == "+") {
            mat[r][c] += value;
        } else static if (op == "-") {
            mat[r][c] -= value;
        } else static if (op == "*") {
            mat[r][c] *= value;
        } else static if (op == "/") {
            mat[r][c] /= value;
        }
        return mat[r][c];
    }

    // override Matrix[r, c] = value
    // must be defined before Matrix[r, c] operator
    T opIndexAssign(T value, size_t r, size_t c)
    in {
        assert(r < mat.length && r >= 0);
        if (mat.length > 0) {
            assert(c < mat[0].length && c >= 0);
        }
    }
    do {
        mat[r][c] = value;
        return mat[r][c];
    }

    // override Matrix[r, c] operator
    T opIndex(size_t r, size_t c) const
    in {
        assert(r < mat.length && r >= 0);
        if (mat.length > 0) {
            assert(c < mat[0].length && c >= 0);
        }
    }
    do {
        return mat[r][c];
    }

    // return number of rows in matrix
    size_t numrows() const {
        return mat.length;
    }

    // return number of columns in matrix
    size_t numcols() const {
        return (mat.length > 0 ? mat[0].length : 0);
    }

    // prints the matrix
    void print() const {
        foreach (const T[] r; mat) {
            writeln(r);
        }
    }

    // returns the matrix as a human-readable string
    override string toString() const {
        string result = "";
        foreach (const T[] r; mat) {
            result ~= to!(string)(r) ~ "\n";
        }
        return result;
    }

    private T[][] mat;
    private T det_aux();

    // returns the sum of two matrices
    Matrix!(T) add(const ref Matrix!(T) rhs) const
    in {
        assert(this.numrows() == rhs.numrows() && this.numcols() == rhs.numcols());
    } 
    do {
        T[][] sum = dup2D!(T)(mat);
        foreach (size_t i; 0..sum.length) {
            foreach (size_t j; 0..sum[0].length) {
                sum[i][j] += rhs.mat[i][j];
            }
        }
        return new Matrix!(T)(sum);
    }

    // returns the product of two matrices
    Matrix!(T) multiply(const ref Matrix!(T) rhs) const
    in {
        assert(this.numcols() == rhs.numrows());
    } 
    do {
        Matrix!(T) product = new Matrix!(T)(this.numrows(), rhs.numcols());

        // naive multiplication
        size_t lhsNumRows = this.numrows();
        size_t lhsNumCols = this.numcols();
        size_t rhsNumCols = rhs.numcols();
        /*
        foreach (size_t i; 0..lhsNumRows) {
        	foreach (size_t j; 0..rhsNumCols) {
        		foreach (size_t k; 0..lhsNumCols) {
        			product.mat[i][j] += this.mat[i][k] * rhs.mat[k][j];
                }
            }
        }
        */
        for (size_t i = 0; i < lhsNumRows; i++) {
        	for (size_t j = 0; j < rhsNumCols; j++) {
        		for (size_t k = 0; k < lhsNumCols; k++) {
        			product.mat[i][j] += this.mat[i][k] * rhs.mat[k][j];
                }
            }
        }

        return product;
    }

    Matrix!T deleteRow(size_t row) const
    in {
        assert(row < mat.length && row >= 0);
    }
    do {
        const T[][] copy = mat[0..row]~mat[row+1..mat.length];
        return new Matrix!(T)(copy);
    }

    Matrix!T deleteCol(size_t col) const 
    in {
        if (mat.length > 0) {
            assert(col < mat[0].length && col >= 0);
        }
    }
    do {
        T[][] copy;
        copy.length = mat.length;
        for (size_t r = 0; r < mat.length; r++) {
            copy[r].length = mat[r].length - 1;
            for (size_t c = 0; c < mat[0].length; c++) {
                if (c < col) {
                    copy[r][c] = mat[r][c];
                }
                else if (c > col) { // c > col
                    copy[r][c-1] = mat[r][c];
                }
            }
            
        }
        return new Matrix!(T)(copy);
    }

    T determinant() const {
        T det = T(0);
        if (numrows() == 0) {
            return 0;
        }
        else if (numrows() == 1) {
            return mat[0][0];
        }
        else if (numrows() == 2) {
            return mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
        }
        else {
            for (size_t k = 0; k < numcols(); k++) {
                det += (k % 2 == 0 ? 1 : -1) * mat[0][k] * this.deleteRow(0).deleteCol(k).determinant();
            }
        }
        return det; 
    }

    // override +=, *= operators
    Matrix!(T) opOpAssign(string op)(Matrix!(T) rhs) {
        static if (op == "+") { 
            this.mat = this.add(rhs).mat;
            return this;
        }
        else static if (op == "*") {
            this.mat = this.multiply(rhs).mat;
            return this;
        }
        else { 
            static assert(0, "Operator "~op~" not valid");
        }
    }

    // override +, * operators
    Matrix!(T) opBinary(string op)(Matrix!(T) rhs) const {
        static if (op == "+") { 
            return this.add(rhs);
        }
        else static if (op == "*") {
            return this.multiply(rhs);
        }
        else { 
            static assert(0, "Operator "~op~" not valid");
        }
    }

}

void main() {
    
    float[][] matArr = [ [-5, 0, -1],
                         [1, 2, -1],
                         [-3, 4, 1] ];
    
    Matrix!(float) mat1 = new Matrix!(float)(matArr);
    auto mat2 = new Matrix!(float)(matArr);

    writefln("%d x %d matrix: ", mat1.numrows(), mat1.numcols(), mat1);
    mat1.print();

    writeln("mat2 added to mat");
    mat2 = mat1.add(mat2);

    writeln("mat1:");
    mat1.print();

    writeln("mat2:");
    mat2.print();

    mat2[2, 2] = 3;
    writeln(mat2[2, 2]);
    mat2.print();

    auto identity = Matrix!(float).identity(3);

    writefln("\nmat1 + mat2 = \n%s", mat1 + mat2);
    writefln("\nmat1 * mat2 = \n%s", mat1 * mat2);
    writefln("\nmat1 * identity = \n%s\n", mat2 * identity);

    
    for (size_t i = 0; i < 32; i++) {
        writefln("%dx%d identity * identity = \n%s", i, i, Matrix!(float).identity(i) * Matrix!(float).identity(i));
    }
    

    size_t matSize = 3;
    Matrix!(float) a = new Matrix!(float)(matArr);
    Matrix!(float) b = new Matrix!(float)(matSize, matSize, 1);

    // determinant timing
    SysTime started = Clock.currTime();
    float det = a.determinant();
    writefln("The determinant of \n%sis %f\n", a, det);

    // multiplication timing
    matSize = 10;
    a = new Matrix!(float)(matSize, matSize, 1);
    b = new Matrix!(float)(matSize, matSize, 1);
    Matrix!(float) c;

    if (matSize <= 50) {
        writefln("Multiplying\n%sby\n%s", a, b);
    }
    started = Clock.currTime();
    c = a * b;
    writefln("Finished multiplying after %s", Clock.currTime() - started);
    if (matSize <= 50) {
        writefln("The product is:\n%s\n", c);
    }
}