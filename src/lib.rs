/// Matrix Object
#[derive(Debug)]
pub struct Matrix{
    pub matrix: Vec<Vec<f64>>,
}

impl Matrix {
    /// Allocates a new `Matrix<f64>`, and moves `initial_matrix`'s items into it
    ///
    /// `intial_matrix` is in the form of Vec<Vec<f64>>, where the inner `Vec<f64>` is
    /// each row of a matrix, and the length of `Vec<f64>` is how many columns in the `Matrix`.
    ///
    /// ### Examples
    /// ```rust
    /// use reduced_row_echelon_form::Matrix;
    /// let matrix = vec![
    ///     vec![1.0, 3.0],
    ///     vec![2.0, 1.5],
    ///     vec![-2.0, -1.5],
    /// ];
    /// let matrix = Matrix::from(matrix);
    /// ```
    pub fn from(initial_matrix: Vec<Vec<f64>>)-> Self{
        Self{ matrix: initial_matrix}
    }
    /// Returns the inverse of the `pivot point` if finite, otherwise `panics`
    #[inline]
    fn calc_inverse_pivot_point(pivot_point: f64) -> f64 {
        let row_scalar = 1.0 / pivot_point;
        if f64::is_finite(row_scalar) {
            row_scalar
        }else {
            panic!("Invalid row scalar for this value: {pivot_point}");
        }
    }
}

impl Matrix {


    /// consumes the matrix and returns its Reduced Row Echelon Form(or as close as it can)
    ///
    /// ### Algorithm
    /// Step 1. Get the first occurrence of the leftmost nonzero in the current column(giving us the target row) \
    /// Step 2. Scale the target row \
    /// Step 3. Zero out the current column based on the target rows values \
    /// Step 4. Move the target row to the "top" \
    ///
    pub fn to_reduced_row_echelon_form(mut self) -> Self{
        let mut current_col = 0;
        for current_row in 0..self.matrix.len(){
            let pivot_point = self.get_leftmost_nonzero_in_a_col(current_col);

            if pivot_point == usize::MAX{ continue; }

            self.scale_row_to_one(current_col, pivot_point);

            self.zero_a_column(current_col, pivot_point);

            if pivot_point != current_row{ self.swap_rows(pivot_point, current_row); }

            current_col += 1;
        }
        self
    }
    /// Returns the index(column) of the first `non-zero` number in the Vector, starting from the specified row
    fn get_leftmost_nonzero_in_a_row(&self, starting_row: usize) -> usize {
        for i in 0..self.matrix[0].len() {
            if self.matrix[starting_row][i] != 0.0 {
                return i
            }
        }
        usize::MAX
    }
    fn get_leftmost_nonzero_in_a_col(&self, col: usize) -> usize{
        for i in 0..self.matrix.len(){
            // if the leftmost nonzero number in a row equals the row we are in, then return
            // the row number(i)
            if self.get_leftmost_nonzero_in_a_row(i) == col{
                return i
            }
        }
        usize::MAX
    }

    fn zero_a_column(&mut self, target_column: usize, pivot_position: usize){
        for rows in 0..self.matrix.len(){
            if self.matrix[rows][target_column] != 0.0 {
                if rows != pivot_position{
                    self.replacement_addition(rows, pivot_position, target_column);
                }
            }
        }
    }
    /// Adds/Subtracts a scalar of a source row from a specified row.
    fn replacement_addition(&mut self, row_to_scale: usize, row_source: usize, starting_col: usize){
        let row_scalar: f64 = self.matrix[row_to_scale][starting_col];
        for i in 0..self.matrix[0].len() {
            self.matrix[row_to_scale][i] -= row_scalar * self.matrix[row_source][i];
        }
    }
    /// Swaps two specified rows of the internal `Matrix`
    fn swap_rows(&mut self, from_row: usize, to_row: usize) {
        //Guard clause
        if from_row == to_row{ return; }
        for i in 0..self.matrix[0].len() {
            let temp = self.matrix[to_row][i];
            self.matrix[to_row][i] = self.matrix[from_row][i];
            self.matrix[from_row][i] = temp;
        }
    }
    /// Scales a whole row of a matrix to one, starting from the specified column.
    fn scale_row_to_one(&mut self, pivot_column: usize, row_to_scale: usize) {

        let row_scalar = Self::calc_inverse_pivot_point(self.matrix[row_to_scale][pivot_column]);
        for i in pivot_column..self.matrix[row_to_scale].len(){
            self.matrix[row_to_scale][i] *= row_scalar;
        }
    }

}
#[cfg(test)]
mod test{
    use super::*;
    #[test]
    fn reduced_row_echelon_form(){
        // | 0.0 | 10.0 | 0.0 |
        // | 0.0 | 5.0 | 2.5 |
        // | 2.0 | 0.0 | 0.0 |
        let matrix = Matrix{ matrix: vec![
            vec![0.0,10.0,0.0],
            vec![0.0, 5.0,2.5],
            vec![2.0, 0.0, 0.0]]
        };

        // | 1.0 | 0.0 | 0.0 |
        // | 0.0 | 1.0 | 0.0 |
        // | 0.0 | 0.0 | 1.0 |
        let in_form_matrix = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0]
        ];
        assert_eq!(matrix.to_reduced_row_echelon_form().matrix, in_form_matrix);

    }

    #[test]
    fn rref_non_square_matrix(){
        // | 2.0 | 2.0 | 0.0 |
        // | 0.0 | 0.0 | 1.0 |
        let matrix = Matrix{ matrix: vec![
            vec![2.0, 2.0, 0.0],
            vec![0.0 , 0.0, 1.0],
        ]};

        // | 1.0 | 1.0 | 0.0 |
        // | 0.0 | 0.0 | 1.0 |
        let in_form_matrix = vec![
            vec![1.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];

        assert_eq!(matrix.to_reduced_row_echelon_form().matrix, in_form_matrix);
    }
    #[test]
    fn rref_zero_matrix(){
        // | 0.0 | 0.0 |
        // | 0.0 | 0.0 |
        let matrix = Matrix{ matrix: vec![vec![0.0, 0.0], vec![0.0, 0.0]]};

        let in_form_matrix =vec![vec![0.0, 0.0], vec![0.0, 0.0]];

        assert_eq!(matrix.to_reduced_row_echelon_form().matrix, in_form_matrix);
    }

    #[test]
    fn zero_first_column(){
        // | 2.0 | 1.0 |
        // | 1.0 | 2.0 |
        let mut matrix = Matrix{matrix: vec![
            vec![2.0, 1.0],
            vec![1.0, 2.0]]
        };
        matrix.zero_a_column(0,1);

        // | 0.0 | -3.0 |
        // | 1.0 |  2.0 |
        let target_matrix = vec![
            vec![0.0, -3.0],
            vec![1.0, 2.0]
        ];
        assert_eq!(matrix.matrix,target_matrix);
    }
    #[test]
    fn test_replacement_addition(){
        // | 2.0 | 1.0 |
        // | 1.0 | 2.0 |
        let mut matrix = Matrix{matrix: vec![
            vec![2.0, 1.0],
            vec![1.0, 2.0]
        ]};

        matrix.replacement_addition(0, 1, 0);
        // | 0.0 | -3.0 |
        // | 1.0 | 2.0 |
        let expected_matrix = vec![
            vec![0.0,-3.0],
            vec![1.0,2.0]];

        assert_eq!(matrix.matrix, expected_matrix);
    }
    #[test]
    fn calculate_row_scalar(){
        assert_eq!(Matrix::calc_inverse_pivot_point(5.0), 1.0/5.0);

    }
    #[test]
    #[should_panic]
    fn calculate_invalid_row_scalar(){
        Matrix::calc_inverse_pivot_point(0.0);
    }
    #[test]
    fn test_get_leftmost_zero(){
        // | 0.0 | 5.0 | 0.0 |
        // | 0.0 | 4.0 | 0.0 |
        // | 10.0 | 0.0 | 0.0 |
        let matrix = Matrix::from(vec![
            vec![0.0, 5.0, 0.0],
            vec![0.0, 4.0, 0.0],
            vec![10.0, 0.0, 0.0]
        ]);

        // in the first column(matrix[0]) the first occurrence of a non-zero answer is 10.0
        assert_eq!(matrix.get_leftmost_nonzero_in_a_col(0), 2);

        // because there is no leading non-zero it should return the starting row
        assert_eq!(matrix.get_leftmost_nonzero_in_a_col( 1), 0);

        // This is making sure that it is the first
        assert_eq!(matrix.get_leftmost_nonzero_in_a_col( 2), usize::MAX );
    }
    #[test]
    fn scale_row_to_one_test(){
        // | 0.0 | 5.0 | 0.0 |
        // | 10.0 | 0.0 | 2.0 |
        let mut matrix = Matrix{ matrix: vec![
            vec![0.0,5.0, 0.0],
            vec![10.0, 0.0, 2.0 ],
        ]};

        // matrix[0] so we only test the leading zeros(at most equal to number of rows)
        for i in 0..matrix.matrix.len(){
            let leftmost_nonzero = matrix.get_leftmost_nonzero_in_a_col(i);

            if leftmost_nonzero != i {
                matrix.scale_row_to_one(leftmost_nonzero, i);
            }
        }
        assert_eq!(matrix.matrix, vec![
            vec![0.0, 1.0, 0.0],
            vec![1.0, 0.0, 0.2],
        ]);
    }
    #[test]
    fn swap_rows_test(){
        // | 0.0 | 5.0 |
        // | 10.0 | 0.0 |
        let mut matrix = Matrix{ matrix: vec![
            vec![0.0, 5.0],
            vec![10.0, 0.0]] };
        matrix.swap_rows(0,1);
        // | 10.0 | 0.0 |
        // | 0.0 | 5.0 |
        let swapped_matrix =   vec![
            vec![10.0, 0.0],
            vec![0.0, 5.0]
        ];
        assert_eq!(matrix.matrix,swapped_matrix);
    }

    #[test]
    fn already_in_row_echelon(){

        let matrix = Matrix{ matrix: vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0]
        ]};

        let expected = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0]
        ];
        assert_eq!(matrix.to_reduced_row_echelon_form().matrix, expected);
   }
}