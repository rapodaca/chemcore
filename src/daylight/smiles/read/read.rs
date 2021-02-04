use std::collections::HashMap;

use purr::read::{ read as read_tree, Reading };
use purr::parts;
use purr::graph::{ from_tree, Error as GraphError };

use crate::molecule::{ DefaultMolecule };
use super::{ to_node, kekulize, Error };

pub fn read(
    smiles: &str, map: Option<&mut HashMap<usize, u16>>
) -> Result<DefaultMolecule, Error> {
    let Reading { root, trace } = read_tree(smiles)?;
    let mut atoms = match from_tree(root) {
        Ok(atoms) => atoms,
        Err(error) => match error {
            GraphError::IncompatibleJoin(first, second) =>
                return Err(Error::IncompatibleJoin(trace[first], trace[second]))
        }
    };

    kekulize(&mut atoms)?;

    let mut nodes = Vec::new();

    if let Some(map) = map {
        for (i, atom) in atoms.iter().enumerate() {
            if let parts::AtomKind::Bracket { map: klass, .. } = &atom.kind {
                if let Some(klass) = klass {
                    map.insert(i, klass.into());
                }
            }
        }
    }

    for i in 0..atoms.len() {
        nodes.push(to_node(i, &atoms, &trace)?)
    }

    Ok(DefaultMolecule::new(nodes))
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use crate::molecule::{ Node, Atom, Bond, Element, Parity };
    use super::*;

    #[test]
    fn invalid_character() {
        assert_eq!(read("CCX", None), Err(Error::Character(2)))
    }

    #[test]
    fn end_of_line() {
        assert_eq!(read("C=", None), Err(Error::EndOfLine))
    }

    #[test]
    fn stray_directional_bond() {
        assert_eq!(read("C-C/C", None), Err(Error::BondKind(4)))
    }

    #[test]
    fn unkekulizable() {
        assert_eq!(read("ccc", None), Err(Error::Kekulization))
    }

    #[test]
    fn hypervalence() {
        assert_eq!(read("C=C(C)(C)C", None), Err(Error::Valence(2)))
    }

    #[test]
    fn incompatible_join() {
        assert_eq!(read("C-1CC=1", None), Err(Error::IncompatibleJoin(0, 4)))
    }

    #[test]
    fn carbon_with_map() {
        let mut map = HashMap::new();

        read("C-[C:42]", Some(&mut map)).unwrap();

        assert_eq!(map, vec![ (1, 42) ].into_iter().collect::<HashMap<_,_>>())
    }

    #[test]
    fn organic_star() {
        assert_eq!(read("C*", None), Ok(DefaultMolecule::new(vec![
            Node {
                atom: Atom {
                    element: Some(Element::C),
                    electrons: 0,
                    hydrogens: 3,
                    isotope: None,
                    parity: None
                },
                bonds: vec![ Bond::new(2, None, 1) ]
            },
            Node {
                atom: Atom {
                    element: None,
                    electrons: 0,
                    hydrogens: 0,
                    isotope: None,
                    parity: None
                },
                bonds: vec![ Bond::new(2, None, 0) ]
            },
        ])))
    }

    #[test]
    fn bracket_star() {
        assert_eq!(read("C[*]", None), Ok(DefaultMolecule::new(vec![
            Node {
                atom: Atom {
                    element: Some(Element::C),
                    electrons: 0,
                    hydrogens: 3,
                    isotope: None,
                    parity: None
                },
                bonds: vec![ Bond::new(2, None, 1) ]
            },
            Node {
                atom: Atom {
                    element: None,
                    electrons: 0,
                    hydrogens: 0,
                    isotope: None,
                    parity: None
                },
                bonds: vec![ Bond::new(2, None, 0) ]
            },
        ])))
    }

    #[test]
    fn ethane() {
        assert_eq!(read("CC", None), Ok(DefaultMolecule::new(vec![
            Node {
                atom: Atom {
                    element: Some(Element::C),
                    electrons: 0,
                    hydrogens: 3,
                    isotope: None,
                    parity: None
                },
                bonds: vec![ Bond::new(2, None, 1) ]
            },
            Node {
                atom: Atom {
                    element: Some(Element::C),
                    electrons: 0,
                    hydrogens: 3,
                    isotope: None,
                    parity: None
                },
                bonds: vec![ Bond::new(2, None, 0) ]
            },
        ])))
    }

    #[test]
    fn ethene_aromatic_atoms() {
        assert_eq!(read("cc", None), Ok(DefaultMolecule::new(vec![
            Node {
                atom: Atom {
                    element: Some(Element::C),
                    electrons: 0,
                    hydrogens: 2,
                    isotope: None,
                    parity: None
                },
                bonds: vec![ Bond::new(4, None, 1) ]
            },
            Node {
                atom: Atom {
                    element: Some(Element::C),
                    electrons: 0,
                    hydrogens: 2,
                    isotope: None,
                    parity: None
                },
                bonds: vec![ Bond::new(4, None, 0) ]
            },
        ])))
    }

    #[test]
    fn ethene_aromatic_bond() {
        assert_eq!(read("C:C", None), Ok(DefaultMolecule::new(vec![
            Node {
                atom: Atom {
                    element: Some(Element::C),
                    electrons: 0,
                    hydrogens: 2,
                    isotope: None,
                    parity: None
                },
                bonds: vec![ Bond::new(4, None, 1) ]
            },
            Node {
                atom: Atom {
                    element: Some(Element::C),
                    electrons: 0,
                    hydrogens: 2,
                    isotope: None,
                    parity: None
                },
                bonds: vec![ Bond::new(4, None, 0) ]
            },
        ])))
    }

    #[test]
    fn trans_butene() {
        assert_eq!(read("C/C=C/C", None), Ok(DefaultMolecule::new(vec![
            Node {
                atom: Atom {
                    element: Some(Element::C),
                    electrons: 0,
                    hydrogens: 3,
                    isotope: None,
                    parity: None
                },
                bonds: vec![ Bond::new(2, None, 1) ]
            },
            Node {
                atom: Atom {
                    element: Some(Element::C),
                    electrons: 0,
                    hydrogens: 1,
                    isotope: None,
                    parity: None
                },
                bonds: vec![
                    Bond::new(2, None, 0),
                    Bond::new(4, Some(Parity::Negative), 2)
                ]
            },
            Node {
                atom: Atom {
                    element: Some(Element::C),
                    electrons: 0,
                    hydrogens: 1,
                    isotope: None,
                    parity: None
                },
                bonds: vec![
                    Bond::new(4, Some(Parity::Negative), 1),
                    Bond::new(2, None, 3)
                ]
            },
            Node {
                atom: Atom {
                    element: Some(Element::C),
                    electrons: 0,
                    hydrogens: 3,
                    isotope: None,
                    parity: None
                },
                bonds: vec![ Bond::new(2, None, 2) ]
            }
        ])))
    }
}