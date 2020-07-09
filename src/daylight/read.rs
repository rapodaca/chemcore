use crate::molecule::{ DefaultMolecule, Error };
use super::smiles_to_adjacency::smiles_to_adjacency;

/// Returns a molecule given a SMILES-encoded string. Errors given:
/// - syntax error
/// - hypervalence
/// - impossible isotopes
/// 
/// ```rust
/// use chemcore::daylight::read;
/// use chemcore::molecule::{ Molecule, Error };
/// 
/// fn main() -> Result<(), Error> {
///     let molecule = read(&"c1ccccc1")?;
/// 
///     assert_eq!(molecule.hydrogens(0), Ok(1));
/// 
///     Ok(())
/// }
/// ```
pub fn read(smiles: &str) -> Result<DefaultMolecule, Error> {
    let adjacency = smiles_to_adjacency(smiles)?;
    
    DefaultMolecule::from_adjacency(adjacency)
}

#[cfg(test)]
mod tests {
    use super::*;
    use purr::read::Error as PurrError;
    use gamma::graph::Graph;
    use crate::molecule::{ Molecule, Parity };

    #[test]
    fn syntax_error() {
        let result = read(&"CX");

        assert_eq!(result, Err(Error::Purr(PurrError::InvalidCharacter(1))));
    }

    #[test]
    fn can_not_kekulize() {
        let result = read(&"ccc");

        assert_eq!(result, Err(Error::CanNotKekulize));
    }

    #[test]
    fn hypervalent() {
        let result = read(&"C(C)(C)(C)(C)C");

        assert_eq!(result, Err(Error::Hypervalent(0)));
    }

    #[test]
    fn impossible_isotope() {
        let result = read(&"[2C]");

        assert_eq!(result, Err(Error::ImpossibleIsotope(0)));
    }

    #[test]
    fn methane() {
        let molecule = read(&"C").unwrap();

        assert_eq!(molecule.hydrogens(0), Ok(4));
    }

    #[test]
    fn benzene_given_aromatic_bonds() {
        let molecule = read(&"C1:C:C:C:C:C:1").unwrap();
        let bonds = molecule.edges().iter().map(|
            (sid, tid)| (*sid, *tid, molecule.bond_order(*sid, *tid).unwrap())
        ).collect::<Vec<_>>();

        assert_eq!(bonds, vec![
            (0, 5, 2f32),
            (0, 1, 1f32),
            (1, 2, 2f32),
            (2, 3, 1f32),
            (3, 4, 2f32),
            (4, 5, 1f32)
        ]);
    }

    #[test]
    fn benzene_given_aromatic_atoms() {
        let molecule = read(&"c1ccccc1").unwrap();
        let bonds = molecule.edges().iter().map(|
            (sid, tid)| (*sid, *tid, molecule.bond_order(*sid, *tid).unwrap())
        ).collect::<Vec<_>>();

        assert_eq!(bonds, vec![
            (0, 5, 2f32),
            (0, 1, 1f32),
            (1, 2, 2f32),
            (2, 3, 1f32),
            (3, 4, 2f32),
            (4, 5, 1f32)
        ]);
    }

    #[test]
    fn bromochlorofluoromethane() {
        let molecule = read(&"[C@H](Br)(Cl)F").unwrap();
        let parity = molecule.atom_parity(0);

        assert_eq!(parity, Ok(Some(Parity::Negative)));
    }
}