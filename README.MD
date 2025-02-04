## reduced_row_echelon_form_jeck

This crate is a published version of my submission to an extra credit assignment for `UMD MATH3326 Vectors and Matrices`.


### **Documentation**
[Library documentation with examples](https://docs.rs/reduced_row_echelon_form_jeck/0.1.0).

### **Usage**
To bring this crate into scope, either add `reduced_row_echelon_form_jeck@0.1.0` to your dependencies in `cargo.toml`, or run `cargo add reduced_row_echelon_form_jeck@0.1.0`.

Here is an example that creates a new Rust project, adds `reduced_row_echelon_form_jeck` , and shows how to use the api to convert a `Matrix` to RREF.

First, create a new directory and create a new cargo instance:
```bash
$ mkdir reduced_row_echelon_form_jeck_example
$ cd reduced_row_echelon_form_jeck_example
$ cargo init
```
Second, add `reduced_row_echelon_form_jeck` to dependancies
```bash
$ cargo add reduced_row_echelon_form_jeck
```
Then, go to `src/main.txt`. Replace the contents with:

```rust
fn main() {
    use reduced_row_echelon_form_jeck::Matrix;
    // Matrices are vectors of Rows
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
```
Finally, run the program with `cargo run `
The output should be, Note: subsequent runs might look a little different:
```bash
$ cargo run
   Compiling reduced_row_echelon_form_jeck v0.1.0 (C:\Users\James\RustroverProjects\reduced_row_echelon_form)
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 0.85s                                                                                                          
     Running `target\debug\reduced_row_echelon_form_jeck.exe`
| 1 0 |
| 0 1 |
| 0 0 |
```
