use purr::read::read;
use purr::valence::implicit_hydrogens;
use purr::mol::{ Atom as PurrAtom, Bond as PurrBond, Style };
use purr::valence::Error as ValenceError;
use crate::molecule::{ Element, Atom, BondOrder, Parity, Error };
use super::bond_parity::bond_parity;
use super::atom_parity::atom_parity;
use super::kekulize;

pub fn smiles_to_adjacency(
    smiles: &str
) -> Result<Vec<(Atom, Vec<(usize, BondOrder, Option<Parity>)>)>, Error> {
    let atoms = kekulize(read(smiles)?)?;
    let mut result = Vec::new();

    for id in 0..atoms.len() {
        let atom = create_atom(id, &atoms)?;
        let mut bonds = Vec::new();

        for purr_bond in &atoms[id].bonds {
            bonds.push(create_bond(id, purr_bond, &atoms)?);
        }

        result.push((atom, bonds));
    }

    Ok(result)
}

fn create_atom(id: usize, atoms: &Vec<PurrAtom>) -> Result<Atom, Error> {
    let atom = &atoms[id];
    let element = Element::from(atom.nub.element);
    let hydrogens = match atom.nub.hcount {
        Some(hcount) => hcount,
        None => match implicit_hydrogens(atom) {
            Ok(count) => count.unwrap_or_default(),
            Err(ValenceError::UnmatchableValence) =>
                return Err(Error::Hypervalent(id))
        }
    };
    let ion = atom.nub.charge.unwrap_or_default();
    let isotope = atom.nub.isotope;
    let parity = atom_parity(id, atom);

    Ok(Atom { element, hydrogens, ion, isotope, parity })
}

fn create_bond(
    sid: usize, bond: &PurrBond, atoms: &Vec<PurrAtom>
) -> Result<(usize, BondOrder, Option<Parity>), Error> {
    let PurrBond { tid, style } = bond;
    let order = match style {
        Some(Style::Double) => BondOrder::Double,
        Some(Style::Triple) => BondOrder::Triple,
        Some(Style::Quadruple) => BondOrder::Quadruple,
        Some(Style::Aromatic)=> return Err(Error::UnsupportedBond(sid, *tid)),
        _ => BondOrder::Single
    };

    let parity = match order {
        BondOrder::Double => bond_parity(sid, *tid, atoms)?,
        _ => None
    };

    Ok((*tid, order, parity))
}

#[cfg(test)]
mod tests {
    use super::*;
    use purr::read::Error as PurrError;
    use crate::molecule::Element;

    #[test]
    fn invalid_atom() {
        let result = smiles_to_adjacency(&"X");

        assert_eq!(result, Err(Error::Purr(PurrError::InvalidCharacter(0))));
    }

    #[test]
    fn methane() {
        let adj = smiles_to_adjacency(&"C").unwrap();

        assert_eq!(adj, vec![
            (Atom::new(Element::C, 4, 0, None, None), vec![ ])
        ]);
    }

    #[test]
    fn methyl_cation() {
        let adj = smiles_to_adjacency(&"[CH3+]").unwrap();

        assert_eq!(adj, vec![
            (Atom::new(Element::C, 3, 1, None, None), vec![ ])
        ]);
    }

    #[test]
    fn methane_c13() {
        let adj = smiles_to_adjacency(&"[13CH4]").unwrap();

        assert_eq!(adj, vec![
            (Atom::new(Element::C, 4, 0, Some(13), None), vec![ ])
        ]);
    }

    #[test]
    fn ethane() {
        let adj = smiles_to_adjacency(&"CC").unwrap();

        assert_eq!(adj, vec![
            (Atom::new(Element::C, 3, 0, None, None), vec![
                (1, BondOrder::Single, None)
            ]),
            (Atom::new(Element::C, 3, 0, None, None), vec![
                (0, BondOrder::Single, None)
            ])
        ]);
    }

    #[test]
    fn ethane_with_single_bond() {
        let adj = smiles_to_adjacency(&"C-C").unwrap();

        assert_eq!(adj, vec![
            (Atom::new(Element::C, 3, 0, None, None), vec![
                (1, BondOrder::Single, None)
            ]),
            (Atom::new(Element::C, 3, 0, None, None), vec![
                (0, BondOrder::Single, None)
            ])
        ]);
    }

    #[test]
    fn ethene() {
        let adj = smiles_to_adjacency(&"C=C").unwrap();

        assert_eq!(adj, vec![
            (Atom::new(Element::C, 2, 0, None, None), vec![
                (1, BondOrder::Double, None)
            ]),
            (Atom::new(Element::C, 2, 0, None, None), vec![
                (0, BondOrder::Double, None)
            ])
        ]);
    }

    #[test]
    fn ethene_with_aromatic_bond() {
        let adj = smiles_to_adjacency(&"C:C").unwrap();

        assert_eq!(adj, vec![
            (Atom::new(Element::C, 2, 0, None, None), vec![
                (1, BondOrder::Double, None)
            ]),
            (Atom::new(Element::C, 2, 0, None, None), vec![
                (0, BondOrder::Double, None)
            ])
        ]);
    }

    #[test]
    fn propene_with_conformation() {
        let adj = smiles_to_adjacency(&"C=C/C").unwrap();

        assert_eq!(adj, vec![
            (Atom::new(Element::C, 2, 0, None, None), vec![
                (1, BondOrder::Double, None)
            ]),
            (Atom::new(Element::C, 1, 0, None, None), vec![
                (0, BondOrder::Double, None),
                (2, BondOrder::Single, None)
            ]),
            (Atom::new(Element::C, 3, 0, None, None), vec![
                (1, BondOrder::Single, None)
            ])
        ]);
    }

    #[test]
    fn butene_with_trans_conformation() {
        let adj = smiles_to_adjacency(&"C/C=C/C").unwrap();

        assert_eq!(adj, vec![
            (Atom::new(Element::C, 3, 0, None, None), vec![
                (1, BondOrder::Single, None)
            ]),
            (Atom::new(Element::C, 1, 0, None, None), vec![
                (0, BondOrder::Single, None),
                (2, BondOrder::Double, Some(Parity::Negative))
            ]),
            (Atom::new(Element::C, 1, 0, None, None), vec![
                (1, BondOrder::Double, Some(Parity::Negative)),
                (3, BondOrder::Single, None)
            ]),
            (Atom::new(Element::C, 3, 0, None, None), vec![
                (2, BondOrder::Single, None)
            ])
        ]);
    }

    #[test]
    fn butene_with_cis_conformation() {
        let adj = smiles_to_adjacency(&"C/C=C\\C").unwrap();

        assert_eq!(adj, vec![
            (Atom::new(Element::C, 3, 0, None, None), vec![
                (1, BondOrder::Single, None)
            ]),
            (Atom::new(Element::C, 1, 0, None, None), vec![
                (0, BondOrder::Single, None),
                (2, BondOrder::Double, Some(Parity::Positive))
            ]),
            (Atom::new(Element::C, 1, 0, None, None), vec![
                (1, BondOrder::Double, Some(Parity::Positive)),
                (3, BondOrder::Single, None)
            ]),
            (Atom::new(Element::C, 3, 0, None, None), vec![
                (2, BondOrder::Single, None)
            ])
        ]);
    }

    #[test]
    fn bromochlorofluoromethane() {
        let adj = smiles_to_adjacency(&"[C@H](Br)(Cl)F").unwrap();

        assert_eq!(adj, vec![
            (Atom::new(Element::C, 1, 0, None, Some(Parity::Negative)), vec![
                (1, BondOrder::Single, None),
                (2, BondOrder::Single, None),
                (3, BondOrder::Single, None)
            ]),
            (Atom::new(Element::Br, 0, 0, None, None), vec![
                (0, BondOrder::Single, None)
            ]),
            (Atom::new(Element::Cl, 0, 0, None, None), vec![
                (0, BondOrder::Single, None)
            ]),
            (Atom::new(Element::F, 0, 0, None, None), vec![
                (0, BondOrder::Single, None)
            ])
        ])
    }
}