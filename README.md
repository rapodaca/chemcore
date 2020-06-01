# ChemCore

A cheminformatics toolkit for Rust.

ChemCore provides primitives for working with Molecules. For details, see:*[ChemCore: A Cheminformatics Toolkit for Rust](https://depth-first.com/articles/2020/06/01/chemcore-a-cheminformatics-toolkit-for-rust/)*.

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
chemcore = "0.1"
```

# Examples

Parse cyclopropane SMILES and traverse in depth-first order:

```rust
use gamma::traversal::depth_first;
use chemcore::daylight;

fn main() {
    let molecule = daylight::Molecule::build(&"C1CC1").unwrap();
    let traversal = depth_first(&molecule, &0).unwrap();

    assert_eq!(traversal.collect::<Vec<_>>(), vec![
        (&0, &2, false),
        (&2, &1, false),
        (&1, &0, true)
    ]);
}
```

# Versions

ChemCore is not yet stable. Patch versions never introduce breaking changes, but minor/major revisions probably will.

# License

ChemCore is distributed under the terms of the MIT License. See
[LICENSE-MIT](LICENSE-MIT) and [COPYRIGHT](COPYRIGHT) for details.