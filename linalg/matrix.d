
module linalg.matrix;

import std.stdio : writefln, writeln;
import std.conv : to;
import std.random : uniform;
import std.datetime : SysTime, Clock;

/// Returns a deep copy of the given dynamic 2-dimensional array
T[][] dup2D(T)(const ref T[][] arr) {
    T[][] copy;
    copy.length = arr.length;
    for (size_t r = 0; r < arr.length; r++) {
        copy[r].length = arr[r].length;
        copy[r] = arr[r].dup;
    }
    return copy;
}

/**
 * Template parameterized Matrix class
 * Params:
 * T = the type for each matrix cell
 * rows = the number of rows in the matrix
 * cols = the number of columns in the matrix (default = rows)
 * See_Also:
 * dup2D
 */
class Matrix(T, size_t _rows, size_t _cols = _rows) {

public:

    /// default constructor
    this() {

    }

    /// constructs the matrix with all cells set to the parameter 'initialVal'
    this(T initialVal) {
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                mat[i][j] = initialVal;
            }
        }
    }

    /// copy constructor: constructs a copy of the given matrix
    this(const ref Matrix!(T, rows, cols) rhs) {
        mat = multidup!(T)(rhs.mat);
    }

    /// constructs a matrix given a 2-dimensional static array
    this(const ref T[cols][rows] a) {
        mat = multidup!(T)(a);
    }

    /// override Matrix[r, c] op= value
    // must be defined before Matrix[r, c] operator
    T opIndexOpAssign(string op)(T value, size_t r, size_t c)
    in {
        assert(r < rows && r >= 0, "Matrix::opIndexOpAssign Row: "~to!string(r)~" not within valid range [0.."~to!string(rows)~"]");
        if (rows > 0) {
            assert(c < cols && c >= 0,"Matrix::opIndexOpAssign Column: "~to!string(c)~" not within valid range [0.."~to!string(cols)~"]");
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


    /// override Matrix[r, c] = value
    // must be defined before Matrix[r, c] operator
    T opIndexAssign(T value, size_t r, size_t c)
    in {
        assert(r < rows && r >= 0, "Matrix::opIndexAssign Row: "~to!string(r)~" not within valid range [0.."~to!string(rows)~"]");
        if (rows > 0) {
            assert(c < cols && c >= 0, "Matrix::opIndexAssign Column: "~to!string(c)~" not within valid range [0.."~to!string(cols)~"]");
        }
    }
    do {
        mat[r][c] = value;
        return mat[r][c];
    }

    /// override Matrix[r, c] operator
    T opIndex(size_t r, size_t c) const
    in {
        assert(r < rows && r >= 0, "Matrix::opIndex Row: "~to!string(r)~" not within valid range [0.."~to!string(rows)~"]");
        if (rows > 0) {
            assert(c < cols && c >= 0, "Matrix::opIndex Column: "~to!string(c)~" not within valid range [0.."~to!string(cols)~"]");
        }
    }
    do {
        return mat[r][c];
    }

    /// override +=, *= operators for scalars
    Matrix!(T, rows, cols) opOpAssign(string op)(T scalar) {
        static if (op == "+") { 
            this.mat = this.add(rhs).mat;
            return this;
        }
        else static if (op == "*") {
            this.mat = this.multiply(scalar).mat;
            return this;
        }
        else { 
            static assert(0, "Operator "~op~" not valid");
        }
    }


    /** 
    override +=, *= operators for matrices
    only works for multiplications resulting in a matrix of the same dimensions as the lhs
    */
    Matrix!(T, rows, cols) opOpAssign(string op)(const ref Matrix!(T, rows, cols) rhs) {
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

    /// override * operator for scalars
    Matrix!(T, rows, cols) opBinary(string op : "*")(T scalar) const {
        static if (op == "*") {
            return this.multiply(scalar);
        }
    }

    /// override +, * operators for matrices
    Matrix!(T, rows, M.cols) opBinary(string op, M)(M rhs) const {
        static if (op == "+") { 
            return this.add(rhs);
        }
        else static if (op == "*") {
            return this.multiply!(M)(rhs);
        }
        else { 
            static assert(0, "Operator "~op~" not valid");
        }
    }
    
    /// return number of rows in matrix
    static size_t numrows() {
        return rows;
    }

    /// return number of columns in matrix
    static size_t numcols() {
        return cols;
    }

    /// prints the matrix
    void print() const {
        foreach (const T[] r; mat) {
            writeln(r);
        }
    }
    
    /// sets all matrix indices to the given value
    void fill(T val) {
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                mat[i][j] = val;
            }
        }
    }

    /// returns the matrix as a human-readable string
    override string toString() const {
        string result = "";
        foreach (const T[] r; mat) {
            result ~= to!(string)(r) ~ "\n";
        }
        return result;
    }

    /** 
    * sets the matrix to be the identity matrix of size rows x rows.
    * only available for matrices where rows == cols
    */
    static if (rows == cols) {
        void identify() {
            fill(0);
            for (size_t i = 0; i < rows; i++) {
                mat[i][i] = 1;
            }
        }
    }

    static if (rows > 0 && cols > 0 && rows <= 4 && cols <= 4) {
        /// returns a copy of the matrix with the specified row removed
        Matrix!(T, rows-1, cols) deleteRow(size_t row) const
        in {
            assert(row < rows && row >= 0, "Matrix::deleteRow Column: "~to!string(row)~" not within valid range [0.."~to!string(rows)~"]");
        }
        do {
            T[cols][rows-1] copy;
            for (size_t r = 0; r < rows; r++) {
                for (size_t c = 0; c < cols; c++) {
                    if (r < row) {
                        copy[r][c] = mat[r][c];
                    }
                    else if (r > row) {
                        copy[r-1][c] = mat[r][c];
                    }
                }
                
            }
            return new Matrix!(T, rows-1, cols)(copy);
        }

        /// returns a copy of the matrix with the specified column removed
        Matrix!(T, rows, cols-1) deleteCol(size_t col) const 
        in {
            if (rows > 0) {
                assert(col < cols && col >= 0, "Matrix::deleteCol Column: "~to!string(col)~" not within valid range [0.."~to!string(cols)~"]");
            }
        }
        do {
            T[cols-1][rows] copy;
            for (size_t r = 0; r < rows; r++) {
                for (size_t c = 0; c < cols; c++) {
                    if (c < col) {
                        copy[r][c] = mat[r][c];
                    }
                    else if (c > col) {
                        copy[r][c-1] = mat[r][c];
                    }
                }
                
            }
            return new Matrix!(T, rows, cols-1)(copy);
        }
    }

    static if (rows <= 4 && rows == cols) {
        /** 
        returns the determinant of the matrix.
        limited to 4x4 matrices and smaller
        */
        T determinant() const {
            T det = T(0);
            static if (rows == 0) {
                return 0;
            }
            else static if (rows == 1) {
                return mat[0][0];
            }
            else static if (rows == 2) {
                return mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
            }
            else {
                for (size_t k = 0; k < cols; k++) {
                    det += (k % 2 == 0 ? 1 : -1) * mat[0][k] * this.deleteRow(0).deleteCol(k).determinant();
                }
            }
            return det; 
        }
    }

    /** 
    transposes the matrix.
    limited to matrices where rows == cols
    */
    void transpose() {
        T holder;
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < i+1; j++) {
                holder = mat[i][j];
                mat[i][j] = mat[j][i];
                mat[j][i] = holder;
            }
        }
    }

    /// returns a transposed copy of the matrix
    Matrix!(T, cols, rows) transposed() const {
        T[rows][cols] t;
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                t[i][j] = mat[j][i];
            }
        }
        return new Matrix!(T, cols, rows)(t);
    }

    /// returns a rows x rows identity matrix
    static Matrix!(T, rows) identity() {
        Matrix!(T, rows) ident = new Matrix!(T, rows)(0);
        for (size_t i = 0; i < rows; i++) {
            ident[i, i] = 1;
        }
        return ident;
    }

    /// returns a rows x cols matrix of random
    static Matrix!(T, rows, cols) random() {
        Matrix!(T, rows, cols) random = new Matrix!(T, rows, cols);
        foreach (ref T[cols] r; random.mat) {
            foreach (ref T v; r) {
                v = cast(T)(uniform(0, 100) / 10.0f);
            }
        }
        return random;
    }

private:

    T[cols][rows] mat;
    static const size_t rows = _rows;
    static const size_t cols = _cols;

    // returns the sum of two matrices
    Matrix!(T, rows, cols) add(const ref Matrix!(T, rows, cols) rhs) const
    in {
        assert(rows == rhs.numrows() && cols == rhs.numcols(), "Matrix::op+(Matrix) Invalid matrix addition: lhs and rhs do not have the same dimensions");
    } 
    do {
        T[cols][rows] sum = multidup!(T)(mat);
        foreach (size_t i; 0..sum.length) {
            foreach (size_t j; 0..sum[0].length) {
                sum[i][j] += rhs.mat[i][j];
            }
        }
        return new Matrix!(T, rows, cols)(sum);
    }

    // returns the product of two matrices
    Matrix!(T, rows, M.cols) multiply(M)(const ref M rhs) const
    in {
        assert(cols == M.rows, "Matrix::op*(Matrix) Invalid matrix multiplication: lhs.numcols() != rhs.numrows()");
    } 
    do {
        auto product = new Matrix!(T, rows, M.cols)(0);

        // naive multiplication
        // generate code for multiplication of matrices up to 4x4
        static if (rows <= 4 && cols <= 4 && cols == M.rows) {
            static foreach (size_t i; 0..rows) {
                static foreach (size_t j; 0..M.cols) {
                    static foreach (size_t k;  0..cols) {
                        product[i, j] += this.mat[i][k] * rhs.mat[k][j];
                    }
                }
            }
        } else {
            /*
            // basic naive multiplication
            for (size_t i = 0; i < rows; i++) {
                for (size_t j = 0; j < M.cols; j++) {
                    for (size_t k = 0; k < cols; k++) {
                        product[i, j] += this.mat[i][k] * rhs.mat[k][j];
                    }
                }
            }
            */

            size_t i, j, k;
            for ( i = 0; i < rows; i++ )
            {  
                T Aik = this.mat[i][0];
                for ( j = 0; j < rhs.cols; j++ ) {
                    product[i, j] = Aik * rhs.mat[0][j];
                }
                for ( k = 1; k < cols; k++ )
                {  
                    Aik = this.mat[i][k];
                    for ( j = 0; j < rhs.cols; j++ ) {
                        product[i, j] += Aik * rhs.mat[k][j];
                    }
                }
            }


        }

        return product;
    }

    // returns the product of this matrix and a scalar
    Matrix!(T, rows, cols) multiply(const ref T scalar) const {
        auto product = new Matrix!(T, rows, cols)(mat);

        for (size_t i = 0; i < rows; i++) {
        	for (size_t j = 0; j < cols; j++) {
                product[i, j] *= scalar;
            }
        }

        return product;
    }

    // returns a copy of a multidimensional array
    T[cols][rows] multidup(T)(const ref T[cols][rows] arr) const {
        T[cols][rows] copy;
        for (size_t r = 0; r < arr.length; r++) {
            copy[r] = arr[r].dup;
        }
        return copy;
    }

}

void main() {

    void timeTo(string StmtString)(string label = "") { 
        SysTime start, end;
        start = Clock.currTime();
            mixin(StmtString);
        end = Clock.currTime();
        writefln("%s:\n%s\n", label, end-start);
    }

    /*
    float[3][3] matArr3x3 = [ [-5, 0, -1],
                              [1,  2, -1],
                              [-3, 4,  1] ]; // determinant = -40

    float[4][4] matArr = [ [-5, 0, -1, 2],
                           [1,  2, -1, 4],
                           [-3, 4,  1, 6],
                           [0,  4,  7, 8] ];
    
    int[4][3] matArr3by4 = [  [-5, 0, -1, 2],
                                [1,  2, -1, 4],
                                [-3, 4,  1, 6] ];

    int[3][4] matArr4by3 = [  [-5, 0, -1],
                                [1,  2, -1],
                                [-3, 4,  1],
                                [0,  4,  7] ];
    
    auto mat3x4 = new Matrix!(int, 3, 4)(matArr3by4);
    auto mat4x3 = new Matrix!(int, 4, 3)(matArr4by3);

    auto mat4x4 = mat4x3 * mat3x4;
    writefln("%s\n",mat4x4);

    alias mat3 = Matrix!(float, 3, 3);

    mat3 m = new mat3(matArr3x3);

    writefln("Determinant of:\n%sis %f\n", m, m.determinant());
    */

    int[2][2] aArr= [
        [1, 2],
        [3, 4]
    ];

    int[3][2] bArr = [
        [5, 6, 7],
        [8, 9, 10]
    ];

    timeTo!("
    auto A = new Matrix!(int, 2)(aArr);
    auto B = new Matrix!(int, 2, 3)(bArr);

    auto C = A * B;

    C.print();
    ");
    const auto n = new Matrix!(float, 2)(3);
    auto m = new Matrix!(float, 2, 2)(1);

    auto p = m * n;

    writeln(p);
    writeln(m);
    writeln(n);

    p.fill(-1);
    writeln(p);

    p.identify();
    writeln(p);
    p[0, 1] = 9;
    writeln(p);
    writeln(p.transposed());

    p = p.transposed();
    writeln(p);

    writeln(p.deleteCol(0));
    writeln(p.deleteCol(1));

    writeln(p.deleteRow(0));
    writeln(p.deleteRow(1));

    writeln(p);

    const size_t matSize = 500;
    alias matType = float;

    Matrix!(matType, matSize) a;
    Matrix!(matType, matSize) b;

    timeTo!(`a = Matrix!(float, matSize, matSize).identity();`)("matrix a init");
    timeTo!(`b = Matrix!(float, matSize, matSize).identity();`)("matrix b init");
    
    if (matSize <= 50) {
        writefln("a:\n%sb:\n%s\n", a, b);
    }

    Matrix!(matType, matSize, matSize) c;
    timeTo!(`c = a * b;`)("c = a * b");

    //writefln("Finished multiplying a %dx%d matrix by a %dx%d matrix after %s", a.numrows(), a.numcols(), b.numrows(), b.numcols(), ended - started);

    if (matSize <= 50) {
        writefln("The product is:\n%s\n", c);
        writefln("Transposed:\n%s\n", c.transposed());
        writefln("product * 0.1f:\n%s\n", c*0.1f);
    }

}