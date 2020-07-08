# ChemCore

A cheminformatics toolkit for Rust.

ChemCore provides primitives for working with Molecules. For details, see:*[ChemCore: A Cheminformatics Toolkit for Rust](https://depth-first.com/articles/2020/06/01/chemcore-a-cheminformatics-toolkit-for-rust/)*.

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
chemcore = "0.2"
```

# Examples

Parse cyclopropane SMILES and traverse in depth-first order:

```rust
use gamma::traversal::depth_first;
use gamma::graph::Step;
use chemcore::daylight::read;
use chemcore::molecule::Error;

fn main() -> Result<(), Error> {
    let molecule = read(&"C1CC1")?;
    let traversal = depth_first(&molecule, 0).expect("traversal error");

    assert_eq!(traversal.collect::<Vec<_>>(), vec![
        Step::new(0, 2, false),
        Step::new(2, 1, false),
        Step::new(1, 0, true)
    ]);

    Ok(())
}
```

# Versions

ChemCore is not yet stable. Patch versions never introduce breaking changes, but minor/major revisions probably will.

# License

ChemCore is distributed under the terms of the MIT License. See
[LICENSE-MIT](LICENSE-MIT) and [COPYRIGHT](COPYRIGHT) for details.