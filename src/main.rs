
fn main() {
    use reduced_row_echelon_form_jeck::Matrix;
    // Matrices are vectors of Rows, not as columns
    let matrix = vec![
        vec![1.0, 3.0],
        vec![2.0, 1.5],
        vec![-2.0, -1.5],
    ];
    let matrix = Matrix::from(matrix).to_reduced_row_echelon_form();

    let matrix_in_form = vec![
        vec![1.0, 0.0],
        vec![0.0, 1.0],
        vec![0.0, 0.0],
    ];
    assert_eq!(matrix.matrix, matrix_in_form);
    println!("{}", matrix);
}
#[test]
fn reduced_row_echelon(){
    use reduced_row_echelon_form_jeck::Matrix;
    // Matrices are vectors of Rows, not as columns
    let matrix = vec![
        vec![1.0, 3.0],
        vec![2.0, 1.5],
        vec![-2.0, -1.5],
    ];
    let matrix = Matrix::from(matrix).to_reduced_row_echelon_form();


    let matrix_in_form = vec![
        vec![1.0, 0.0],
        vec![0.0, 1.0],
        vec![0.0, 0.0],
    ];
    assert_eq!(matrix.matrix, matrix_in_form)
}
