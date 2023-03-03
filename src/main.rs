use axb::Matrix;

fn main() {
    let test: Matrix<3, 3> = Matrix::<3, 3>::new([
        [ 2.0, -1.0,  2.0],
        [ 2.0,  2.0, -1.0],
        [-1.0,  3.0,  2.0]
    ]);

    test.print_matrix();
}