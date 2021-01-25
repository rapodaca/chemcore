# ChemCore

A cheminformatics toolkit for Rust.

ChemCore provides primitives for working with Molecules. For details, see: *[ChemCore: A Cheminformatics Toolkit for Rust](https://depth-first.com/articles/2020/06/01/chemcore-a-cheminformatics-toolkit-for-rust/)*.

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
chemcore = "0.4.0"
```

# Examples

Parse cyclopropane SMILES, traverse in depth-first order, and query its `Molecule` and `Graph` traits:

```rust
use gamma::graph::{ Graph };
use gamma::traversal::{ DepthFirst, Step };
use chemcore::daylight::{ read_smiles, SmilesInputError };
use chemcore::molecule::{ Atom, Element, Molecule, Error };

fn main() -> Result<(), SmilesInputError> {
    let molecule = read_smiles(&"C1CC1", None)?;
    let traversal = DepthFirst::new(&molecule, 0).expect("traversal error");

    assert_eq!(traversal.collect::<Vec<_>>(), vec![
        Step::new(0, 2, false),
        Step::new(2, 1, false),
        Step::new(1, 0, true)
    ]);

    // Graph trait
    assert_eq!(molecule.degree(0), Ok(2));
    assert_eq!(molecule.size(), 3);
    assert_eq!(molecule.order(), 3);

    // Molecule trait
    assert_eq!(molecule.atom(0), Ok(&Atom {
        isotope: None,
        element: Some(Element::C),
        hydrogens: 2,
        electrons: 0,
        parity: None,
    }));
    assert_eq!(molecule.charge(0), Ok(0.));
    assert_eq!(molecule.bond_order(0, 1), Ok(1.));

    Ok(())
}
```

# Versions

ChemCore is not yet stable. Patch versions never introduce breaking changes, but minor/major revisions probably will.

# License

ChemCore is distributed under the terms of the MIT License. See
[LICENSE-MIT](LICENSE-MIT) and [COPYRIGHT](COPYRIGHT) for details.