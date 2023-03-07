use axb::matrix::Matrix;
use axb::matrix::transform::mult;

fn main() {
    let test: Matrix<10, 10> = Matrix::<10, 10>::new([
        [9.0, 1.0, 6.0, 7.0, 1.0, 8.0, 2.0, 2.0, 5.0, 1.0],
        [5.0, 9.0, 2.0, 4.0, 6.0, 4.0, 5.0, 8.0, 4.0, 7.0],
        [4.0, 9.0, 4.0, 4.0, 3.0, 3.0, 6.0, 7.0, 5.0, 2.0],
        [3.0, 7.0, 3.0, 8.0, 6.0, 1.0, 7.0, 3.0, 8.0, 7.0],
        [5.0, 7.0, 6.0, 9.0, 9.0, 1.0, 6.0, 7.0, 2.0, 8.0],
        [6.0, 1.0, 5.0, 2.0, 9.0, 9.0, 5.0, 4.0, 6.0, 4.0],
        [8.0, 4.0, 6.0, 4.0, 4.0, 9.0, 2.0, 3.0, 3.0, 1.0],
        [8.0, 7.0, 4.0, 1.0, 9.0, 4.0, 9.0, 2.0, 3.0, 3.0],
        [3.0, 6.0, 1.0, 5.0, 7.0, 2.0, 9.0, 5.0, 4.0, 6.0],
        [1.0, 5.0, 7.0, 6.0, 2.0, 7.0, 1.0, 8.0, 6.0, 5.0]
    ]);

    // test.power().print();
    // println!("{:?}", test.determinant());
    // let vector = Matrix::<3, 1>::new([[1.0], [2.0], [3.0]]);
    
    // println!("{:?}", mult(&test.matrice, &vector.matrice));
}