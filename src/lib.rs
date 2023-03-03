#[derive(Debug)]
pub struct Matrix<const R:usize, const C: usize> {
    matrice: Vec<Vec<f32>>
}

impl <const R:usize, const C:usize> Matrix<R, C> {
    pub fn new(matrice: [[f32; C]; R]) -> Self {
        let mut matrix : Vec<Vec<f32>> = vec![];

        for row_idx in 0..matrice.len() {
            /* Current row */
            let mut temp_row: Vec<f32> = vec![];
            
            /* Add elements in current Row */
            for col_idx in 0..matrice[0].len() {
                temp_row.push(matrice[row_idx][col_idx]);
            }
            
            /* Add current row */
            matrix.push(temp_row);
        }
        Self { matrice: matrix }
    }

    pub fn new_from_vec(v: Vec<Vec<f32>>) -> Self {
        let mut vc: Vec<Vec<f32>> = vec![];

        for row_idx in 0..v.len() {
            let mut temp_row: Vec<f32> = vec![];
            for col_idx in 0..v.len() {
                temp_row.push(v[row_idx][col_idx]);
            }
            vc.push(temp_row);
        }

        Self { matrice: vc }
    }

    pub fn transpose(&self) -> Matrix::<C, R> {
        let mut vc: Vec<Vec<f32>> = vec![];

        for row_idx in 0..self.matrice.len() {
            let mut v: Vec<f32> = vec![];
            for j in 0..self.matrice[row_idx].len() {
                v.push(self.matrice[j][row_idx]);
            }
            vc.push(v);
        }

        Matrix::<C, R>::new_from_vec(vc)
    }

    pub fn print_matrix(&self) -> () {
        for row in &self.matrice {
            for el in row {
                match el < &0.0 {
                    true => print!("   {}", el),
                    false => print!("    {}", el)
                }
            }
            println!();
        }
    }

    pub fn is_square_matrix(&self) -> bool {
        R == C
    }

    pub fn determinant(&self) -> f32 {
        let mut upper_triangular_matrix = self.matrice.clone();
        let n = self.matrice.len();
    
        for i in 0..n {
            let pivot = upper_triangular_matrix[i][i];
            if pivot == 0.0 {
                return 0.0;
            }
            for j in (i+1)..n {
                let multiplier = upper_triangular_matrix[j][i] / pivot;
                for k in i..n {
                    upper_triangular_matrix[j][k] -= multiplier * upper_triangular_matrix[i][k];
                }
            }
        }
    
        let mut det = 1.0;
        for i in 0..n {
            det *= upper_triangular_matrix[i][i];
        }
        det
    }

    pub fn inverse(&self) -> Matrix::<R, C> {
        // Create an augmented matrix by concatenating the identity matrix
        // to the right of the original matrix
        let mut augmented_matrix = self.matrice.clone();
        let n = self.matrice.len();
        for i in 0..n {
            for j in 0..n {
                augmented_matrix[i].push(if i == j { 1.0 } else { 0.0 });
            }
        }

        // Perform row operations to transform the left side of the augmented
        // matrix into the identity matrix
        for i in 0..n {
            // If the pivot element is zero, swap the current row with a row below it
            if augmented_matrix[i][i] == 0.0 {
                let mut j = i + 1;
                while j < n && augmented_matrix[j][i] == 0.0 {
                    j += 1;
                }
                if j == n {
                    // If all remaining elements in the column are zero, the matrix is not invertible
                    panic!("Matrix is not invertible");
                }
                augmented_matrix.swap(i, j);
            }

            // Divide the pivot row by the pivot element to make the pivot element equal to 1
            let pivot = augmented_matrix[i][i];
            for j in i..(2 * n) {
                augmented_matrix[i][j] /= pivot;
            }

            // Subtract a multiple of the pivot row from all other rows to make all other
            // elements in the column equal to zero
            for j in 0..n {
                if j != i {
                    let factor = augmented_matrix[j][i];
                    for k in i..(2 * n) {
                        augmented_matrix[j][k] -= factor * augmented_matrix[i][k];
                    }
                }
            }
        }

        // Extract the right side of the augmented matrix, which is the inverse of the original matrix
        let mut inverted = Vec::new();
        for i in 0..n {
            inverted.push(augmented_matrix[i][n..].to_vec());
        }

        Matrix::<R, C>::new_from_vec(inverted)
    }

    #[allow(non_snake_case)]
    pub fn QR(&self) -> (Matrix::<R, C>, Matrix::<R, C>) {
        let m = self.matrice.len();
        let n = self.matrice[0].len();
    
        let mut q: Vec<Vec<f32>> = vec![vec![0.0; n]; m];
        let mut r: Vec<Vec<f32>> = vec![vec![0.0; n]; n];
    
        for j in 0..n {
            for i in 0..m {
                q[i][j] = self.matrice[i][j];
            }
            for k in 0..j {
                let dot_product = (0..m).fold(
                    0.0,
                    |acc,
                    l | acc + q[l][j] * r[k][l]);
                for i in 0..m {
                    q[i][j] -= dot_product * r[k][i];
                }
            }
            let norm = (0..m).fold(0.0, 
                |acc,
                l| acc + q[l][j] * q[l][j])
                .sqrt();
            for i in 0..m {
                q[i][j] /= norm;
            }
            for i in j..n {
                r[j][i] = (0..m).fold(0.0,
                    |acc, 
                     l| acc + q[l][j] * self.matrice[l][i]);
            }
        }

        (Matrix::<R, C>::new_from_vec(q), Matrix::<R, C>::new_from_vec(r))
    }

    #[allow(non_snake_case)]
    pub fn LU(&self) -> (Matrix::<R, C>, Matrix::<R, C>) {
        let n = self.matrice.len();
        let mut l = vec![vec![0.0; n]; n];
        let mut u = vec![vec![0.0; n]; n];

        for k in 0..n {
            for i in k..n {
                let sum = (0..k).map(|j| l[i][j] * u[j][k]).sum::<f32>();
                l[i][k] = if i == k { 1.0 } else { self.matrice[i][k] - sum };
            }
            for j in k..n {
                let sum = (0..k).map(|i| l[k][i] * u[i][j]).sum::<f32>();
                u[k][j] = (self.matrice[k][j] - sum) / l[k][k];
            }
        }

        (Matrix::<R, C>::new_from_vec(l), Matrix::<R, C>::new_from_vec(u))
    }

    #[allow(non_snake_case)]
    pub fn LDU(&self) -> (Matrix::<R, C>, Matrix::<R, C>, Matrix::<R, C>) {
        let n = self.matrice.len();
        let mut cloned = self.matrice.clone();
        let mut l = vec![vec![0.0; n]; n];
        let mut d = vec![vec![0.0; n]; n];
        let mut u = vec![vec![0.0; n]; n];

        for i in 0..n {
            d[i][i] = cloned[i][i];
            l[i][i] = 1.0;
        }

        for k in 0..n {
            for i in k + 1..n {
                l[i][k] = cloned[i][k] / d[k][k];
                for j in k + 1..n {
                    cloned[i][j] -= l[i][k] * d[k][j];
                }
            }
            for j in k + 1..n {
                u[k][j] = cloned[k][j] / d[k][k];
            }
        }

        (Matrix::<R, C>::new_from_vec(l), Matrix::<R, C>::new_from_vec(d), Matrix::<R, C>::new_from_vec(u))
    }

}

