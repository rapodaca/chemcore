use std::convert::TryFrom;

use purr::{graph, parts};

use super::{to_bond, Error};
use crate::molecule::{Atom, Element, Node, Parity};

pub fn to_node(
    id: usize,
    atoms: &[graph::Atom],
    trace: &[usize],
) -> Result<Node, Error> {
    let atom = &atoms[id];
    let mut bonds = Vec::new();

    for bond in &atom.bonds {
        bonds.push(to_bond(id, bond, atoms, trace)?)
    }

    match to_atom(atom) {
        Ok(atom) => Ok(Node { atom, bonds }),
        Err(error) => match error {
            AtomError::Isotope => Err(Error::Isotope(trace[id])),
            AtomError::Valence => Err(Error::Valence(trace[id])),
            AtomError::Parity => Err(Error::Parity(trace[id])),
            AtomError::ChargedStar => Err(Error::ChargedStar(trace[id])),
        },
    }
}

fn to_atom(atom: &graph::Atom) -> Result<Atom, AtomError> {
    match &atom.kind {
        parts::AtomKind::Star => star_to_atom(),
        parts::AtomKind::Aliphatic(aliphatic) => {
            bare_to_atom(aliphatic.into(), atom.subvalence(), &atom.bonds)
        }
        parts::AtomKind::Aromatic(aromatic) => {
            bare_to_atom(aromatic.into(), atom.subvalence(), &atom.bonds)
        }
        parts::AtomKind::Bracket {
            isotope,
            symbol,
            hcount,
            charge,
            parity,
            ..
        } => bracket_to_atom(
            isotope,
            symbol,
            hcount,
            charge,
            parity,
            &atom.bonds,
        ),
    }
}

fn star_to_atom() -> Result<Atom, AtomError> {
    Ok(Atom {
        element: None,
        isotope: None,
        electrons: 0,
        hydrogens: 0,
        parity: None,
    })
}

fn bare_to_atom(
    element: Element,
    subvalence: u8,
    bonds: &Vec<graph::Bond>,
) -> Result<Atom, AtomError> {
    let mut valence = subvalence;

    for bond in bonds {
        valence = match valence.checked_add(bond.order()) {
            Some(order) => order,
            None => return Err(AtomError::Valence),
        }
    }

    let electrons = match element.valence_electrons().checked_sub(valence) {
        Some(electrons) => electrons,
        None => return Err(AtomError::Valence),
    };

    Ok(Atom {
        element: Some(element),
        isotope: None,
        electrons,
        hydrogens: subvalence,
        parity: None,
    })
}

fn bracket_to_atom(
    isotope: &Option<parts::Number>,
    symbol: &parts::BracketSymbol,
    hcount: &Option<parts::VirtualHydrogen>,
    charge: &Option<parts::Charge>,
    parity: &Option<parts::Parity>,
    bonds: &Vec<graph::Bond>,
) -> Result<Atom, AtomError> {
    let charge = match charge {
        Some(charge) => charge.into(),
        None => 0,
    };
    let element: Option<Element> = match symbol {
        parts::BracketSymbol::Star => {
            if charge == 0 {
                None
            } else {
                return Err(AtomError::ChargedStar);
            }
        }
        parts::BracketSymbol::Aromatic(aromatic) => Some(aromatic.into()),
        parts::BracketSymbol::Element(element) => Some(element.into()),
    };
    let isotope = to_isotope(&element, isotope)?;
    let hydrogens = match hcount {
        Some(hcount) => hcount.into(),
        None => 0,
    };
    let electrons = to_electrons(&element, hydrogens, charge, bonds)?;
    let parity = to_parity(hydrogens, parity, bonds)?;

    Ok(Atom {
        element,
        isotope,
        electrons,
        hydrogens,
        parity,
    })
}

fn to_isotope(
    element: &Option<Element>,
    isotope: &Option<parts::Number>,
) -> Result<Option<u16>, AtomError> {
    let isotope = match isotope {
        Some(isotope) => isotope.into(),
        None => return Ok(None),
    };
    let element = match element {
        Some(element) => element,
        None => return Ok(Some(isotope)),
    };

    if isotope < u16::from(element.atomic_number()) {
        return Err(AtomError::Isotope);
    } else {
        Ok(Some(isotope))
    }
}

fn to_electrons(
    element: &Option<Element>,
    hydrogens: u8,
    ion: i8,
    bonds: &Vec<graph::Bond>,
) -> Result<u8, AtomError> {
    let element = match element {
        Some(element) => element,
        None => return Ok(0),
    };
    let mut bonding = hydrogens as i16;

    for bond_spec in bonds {
        bonding += bond_spec.order() as i16;
    }

    let mut result = element.valence_electrons() as i16;

    result -= ion as i16;
    result -= bonding;

    match u8::try_from(result) {
        Ok(result) => Ok(result),
        Err(_) => Err(AtomError::Valence),
    }
}

fn to_parity(
    hydrogens: u8,
    parity: &Option<parts::Parity>,
    bonds: &Vec<graph::Bond>,
) -> Result<Option<Parity>, AtomError> {
    let parity = match parity {
        Some(parity) => parity,
        None => return Ok(None),
    };

    if hydrogens == 0 {
        if bonds.len() == 4 {
            Ok(Some(parity.into()))
        } else {
            Err(AtomError::Parity)
        }
    } else if hydrogens == 1 {
        if bonds.len() == 3 {
            Ok(Some(parity.into()))
        } else {
            Err(AtomError::Parity)
        }
    } else {
        Err(AtomError::Parity)
    }
}

enum AtomError {
    Valence,
    Isotope,
    Parity,
    ChargedStar,
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use purr::graph::from_tree;
    use purr::read::{read, Reading};

    use super::*;
    use crate::molecule::Bond;

    #[test]
    fn bracket_star_charge() {
        let Reading { root, trace } = read("C-[*+]").unwrap();
        let atoms = from_tree(root).unwrap();
        let result = to_node(1, &atoms, &trace);

        assert_eq!(result, Err(Error::ChargedStar(2)))
    }

    #[test]
    fn lithium_dication() {
        let Reading { root, trace } = read("C-[Li+2]").unwrap();
        let atoms = from_tree(root).unwrap();
        let result = to_node(1, &atoms, &trace);

        assert_eq!(result, Err(Error::Valence(2)))
    }

    #[test]
    fn carbon_h5() {
        let Reading { root, trace } = read("C-[CH5]").unwrap();
        let atoms = from_tree(root).unwrap();
        let result = to_node(1, &atoms, &trace);

        assert_eq!(result, Err(Error::Valence(2)))
    }

    #[test]
    fn carbon_5() {
        let Reading { root, trace } = read("C-[5C]").unwrap();
        let atoms = from_tree(root).unwrap();
        let result = to_node(1, &atoms, &trace);

        assert_eq!(result, Err(Error::Isotope(2)))
    }

    #[test]
    fn carbon_dioxide_parity() {
        let Reading { root, trace } = read("[C@](=O)=O").unwrap();
        let atoms = from_tree(root).unwrap();
        let result = to_node(0, &atoms, &trace);

        assert_eq!(result, Err(Error::Parity(0)))
    }

    #[test]
    fn bare_star() {
        let Reading { root, trace } = read("C*").unwrap();
        let atoms = from_tree(root).unwrap();
        let result = to_node(1, &atoms, &trace);

        assert_eq!(
            result,
            Ok(Node {
                atom: Atom {
                    element: None,
                    isotope: None,
                    electrons: 0,
                    hydrogens: 0,
                    parity: None
                },
                bonds: vec![Bond::new(2, None, 0)]
            })
        )
    }

    #[test]
    fn aromatic_carbon() {
        let Reading { root, trace } = read("Cc").unwrap();
        let atoms = from_tree(root).unwrap();
        let result = to_node(1, &atoms, &trace);

        assert_eq!(
            result,
            Ok(Node {
                atom: Atom {
                    element: Some(Element::C),
                    isotope: None,
                    electrons: 0,
                    hydrogens: 3,
                    parity: None
                },
                bonds: vec![Bond::new(2, None, 0)]
            })
        )
    }

    #[test]
    fn aliphatic_carbon() {
        let Reading { root, trace } = read("CC").unwrap();
        let atoms = from_tree(root).unwrap();
        let result = to_node(1, &atoms, &trace);

        assert_eq!(
            result,
            Ok(Node {
                atom: Atom {
                    element: Some(Element::C),
                    isotope: None,
                    electrons: 0,
                    hydrogens: 3,
                    parity: None
                },
                bonds: vec![Bond::new(2, None, 0)]
            })
        )
    }

    #[test]
    fn bracket_star() {
        let Reading { root, trace } = read("C[*]").unwrap();
        let atoms = from_tree(root).unwrap();
        let result = to_node(1, &atoms, &trace);

        assert_eq!(
            result,
            Ok(Node {
                atom: Atom {
                    element: None,
                    isotope: None,
                    electrons: 0,
                    hydrogens: 0,
                    parity: None
                },
                bonds: vec![Bond::new(2, None, 0)]
            })
        )
    }

    #[test]
    fn bracket_star_12() {
        let Reading { root, trace } = read("[12*]").unwrap();
        let atoms = from_tree(root).unwrap();
        let result = to_node(0, &atoms, &trace);

        assert_eq!(
            result,
            Ok(Node {
                atom: Atom {
                    element: None,
                    isotope: Some(12),
                    electrons: 0,
                    hydrogens: 0,
                    parity: None
                },
                bonds: vec![]
            })
        )
    }

    #[test]
    fn bracket_methane() {
        let Reading { root, trace } = read("[CH4]").unwrap();
        let atoms = from_tree(root).unwrap();
        let result = to_node(0, &atoms, &trace);

        assert_eq!(
            result,
            Ok(Node {
                atom: Atom {
                    element: Some(Element::C),
                    isotope: None,
                    electrons: 0,
                    hydrogens: 4,
                    parity: None
                },
                bonds: vec![]
            })
        )
    }

    #[test]
    fn bracket_methyl_cation() {
        let Reading { root, trace } = read("[CH3+]").unwrap();
        let atoms = from_tree(root).unwrap();
        let result = to_node(0, &atoms, &trace);

        assert_eq!(
            result,
            Ok(Node {
                atom: Atom {
                    element: Some(Element::C),
                    isotope: None,
                    electrons: 0,
                    hydrogens: 3,
                    parity: None
                },
                bonds: vec![]
            })
        )
    }

    #[test]
    fn bracket_methyl_anion() {
        let Reading { root, trace } = read("[CH3-]").unwrap();
        let atoms = from_tree(root).unwrap();
        let result = to_node(0, &atoms, &trace);

        assert_eq!(
            result,
            Ok(Node {
                atom: Atom {
                    element: Some(Element::C),
                    isotope: None,
                    electrons: 2,
                    hydrogens: 3,
                    parity: None
                },
                bonds: vec![]
            })
        )
    }

    #[test]
    fn bracket_carbon_12() {
        let Reading { root, trace } = read("[12C]").unwrap();
        let atoms = from_tree(root);
        let result = to_node(0, &atoms.unwrap(), &trace);

        assert_eq!(
            result,
            Ok(Node {
                atom: Atom {
                    element: Some(Element::C),
                    isotope: Some(12),
                    electrons: 4,
                    hydrogens: 0,
                    parity: None
                },
                bonds: vec![]
            })
        )
    }

    #[test]
    fn bracket_star_parity() {
        let Reading { root, trace } = read("C[*@](F)(Cl)I").unwrap();
        let atoms = from_tree(root);
        let result = to_node(1, &atoms.unwrap(), &trace);

        assert_eq!(
            result,
            Ok(Node {
                atom: Atom {
                    element: None,
                    isotope: None,
                    electrons: 0,
                    hydrogens: 0,
                    parity: Some(Parity::Negative)
                },
                bonds: vec![
                    Bond::new(2, None, 0),
                    Bond::new(2, None, 2),
                    Bond::new(2, None, 3),
                    Bond::new(2, None, 4)
                ]
            })
        )
    }

    #[test]
    fn bracket_carbon_single_single_single_single() {
        let Reading { root, trace } = read("C[C@](F)(Cl)I").unwrap();
        let atoms = from_tree(root).unwrap();
        let result = to_node(1, &atoms, &trace);

        assert_eq!(
            result,
            Ok(Node {
                atom: Atom {
                    element: Some(Element::C),
                    isotope: None,
                    electrons: 0,
                    hydrogens: 0,
                    parity: Some(Parity::Negative)
                },
                bonds: vec![
                    Bond::new(2, None, 0),
                    Bond::new(2, None, 2),
                    Bond::new(2, None, 3),
                    Bond::new(2, None, 4)
                ]
            })
        )
    }

    #[test]
    fn negative_boron_tetravalent() {
        let Reading { root, trace } = read("[B-](C)(C)(C)C").unwrap();
        let atoms = from_tree(root).unwrap();
        let result = to_node(0, &atoms, &trace);

        assert_eq!(
            result,
            Ok(Node {
                atom: Atom {
                    element: Some(Element::B),
                    isotope: None,
                    electrons: 0,
                    hydrogens: 0,
                    parity: None
                },
                bonds: vec![
                    Bond::new(2, None, 1),
                    Bond::new(2, None, 2),
                    Bond::new(2, None, 3),
                    Bond::new(2, None, 4)
                ]
            })
        )
    }
}
