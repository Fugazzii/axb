use axb::matrix::Matrix;
use axb::matrix::transform::mult;

fn main() {
    let test: Matrix<3, 3> = Matrix::<3, 3>::new([
        [ 2.0, -1.0,  2.0],
        [ 2.0,  2.0, -1.0],
        [-1.0,  3.0,  2.0]
    ]);

    let vector = Matrix::<3, 1>::new([[1.0], [2.0], [3.0]]);
    
    println!("{:?}", mult(&test.matrice, &vector.matrice));
}