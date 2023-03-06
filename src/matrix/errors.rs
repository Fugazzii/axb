pub fn non_vector_error(c: &usize) {
    if *c != 1 {
        panic!("Check dimensions of matrix, this is not a vector")
    }
}

pub fn non_square_matrix_error(r: &usize, c: &usize) {
    if *r != *c {
        panic!("Can not make operation. Matrix is non-square")
    }
}