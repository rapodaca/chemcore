use gamma::graph::Graph;
use gamma::matching::{greedy, maximum_matching};
use purr::graph::Atom;
use purr::parts::BondKind;

use super::{pi_subgraph, Error};

pub fn kekulize(atoms: &mut Vec<Atom>) -> Result<(), Error> {
    let pi = pi_subgraph(atoms);
    let mut pairing = greedy(&pi);

    maximum_matching(&pi, &mut pairing);

    if pairing.order() != pi.order() {
        return Err(Error::Kekulization);
    }

    let mut bonds = Vec::new();

    for (sid, tid) in pairing.edges() {
        bonds.push((sid, tid));
        bonds.push((tid, sid));
    }

    for (sid, tid) in bonds {
        let source = &mut atoms[sid];
        let bond = source
            .bonds
            .iter_mut()
            .find(|bond| bond.tid == tid)
            .expect("target bond");

        bond.kind = BondKind::Double;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use purr::graph::{from_tree, Atom, Bond};
    use purr::parts::{Aromatic, AtomKind};
    use purr::read::read;

    use super::*;

    #[test]
    fn unkekulizable() {
        let mut atoms = from_tree(read("ccc").unwrap().root).unwrap();

        assert_eq!(kekulize(&mut atoms), Err(Error::Kekulization))
    }

    #[test]
    fn carbon_aromatic_carbon() {
        let mut atoms = from_tree(read("C:C").unwrap().root).unwrap();

        kekulize(&mut atoms).unwrap();

        assert_eq!(atoms, from_tree(read("C=C").unwrap().root).unwrap())
    }

    #[test]
    fn benzene_aromatic_atoms() {
        let mut atoms = from_tree(read("c1ccccc1").unwrap().root).unwrap();

        kekulize(&mut atoms).unwrap();

        assert_eq!(
            atoms,
            vec![
                Atom {
                    kind: AtomKind::Aromatic(Aromatic::C),
                    bonds: vec![
                        Bond::new(BondKind::Double, 5),
                        Bond::new(BondKind::Elided, 1)
                    ]
                },
                Atom {
                    kind: AtomKind::Aromatic(Aromatic::C),
                    bonds: vec![
                        Bond::new(BondKind::Elided, 0),
                        Bond::new(BondKind::Double, 2)
                    ]
                },
                Atom {
                    kind: AtomKind::Aromatic(Aromatic::C),
                    bonds: vec![
                        Bond::new(BondKind::Double, 1),
                        Bond::new(BondKind::Elided, 3)
                    ]
                },
                Atom {
                    kind: AtomKind::Aromatic(Aromatic::C),
                    bonds: vec![
                        Bond::new(BondKind::Elided, 2),
                        Bond::new(BondKind::Double, 4)
                    ]
                },
                Atom {
                    kind: AtomKind::Aromatic(Aromatic::C),
                    bonds: vec![
                        Bond::new(BondKind::Double, 3),
                        Bond::new(BondKind::Elided, 5)
                    ]
                },
                Atom {
                    kind: AtomKind::Aromatic(Aromatic::C),
                    bonds: vec![
                        Bond::new(BondKind::Elided, 4),
                        Bond::new(BondKind::Double, 0)
                    ]
                }
            ]
        )
    }
}
