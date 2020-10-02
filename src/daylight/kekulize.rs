use purr::mol::{ Atom, Style };
use gamma::matching::{ greedy, maximum_matching };
use gamma::graph::Graph;
use gamma::matching::Pairing;
use crate::molecule::Error;
use super::pi_subgraph;

pub fn kekulize(mut atoms: Vec<Atom>) -> Result<Vec<Atom>, Error> {
    let pi = pi_subgraph(&atoms)?;
    let mut pairing = greedy(&pi);

    maximum_matching(&pi, &mut pairing);

    if pi.is_empty() {
        return Ok(atoms);
    } else if pairing.order() != pi.order() {
        return Err(Error::CanNotKekulize);
    }

    for (sid, mut atom) in atoms.iter_mut().enumerate() {
        atom.nub.aromatic = false;

        if pairing.has_node(sid) {
            restyle_bonds_from_matched(sid, atom, &pairing);
        } else {
            restyle_bonds_from_unmatched(atom);
        }
    }

    Ok(atoms)
}

fn restyle_bonds_from_matched(id: usize, atom: &mut Atom, pairing: &Pairing) {
    let mate = pairing.mate(id);

    for bond in atom.bonds.iter_mut() {
        match bond.style {
            Some(Style::Aromatic) => {
                bond.style.replace(if bond.tid == mate {
                    Style::Double
                } else {
                    Style::Single
                });
            },
            None => {
                if bond.tid == mate {
                    bond.style.replace(Style::Double);
                }
            },
            _ => ()
        }
    }
}

fn restyle_bonds_from_unmatched(atom: &mut Atom) {
    for bond in atom.bonds.iter_mut() {
        match bond.style {
            Some(Style::Aromatic) => {
                bond.style.replace(Style::Single);
            },
            _ => ()
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use super::*;
    use purr::mol::{ Atom, Bond, Element, Nub };
    use purr::read::{ read };

    #[test]
    fn methane_given_aromatic() {
        let atoms = read("c").unwrap();
        let atoms = kekulize(atoms);

        assert_eq!(atoms, Err(Error::CanNotKekulize));
    }   
    
    #[test]
    fn propenyl_radical_given_aromatic_atoms() {
        let atoms = read("ccc").unwrap();
        let atoms = kekulize(atoms);

        assert_eq!(atoms, Err(Error::CanNotKekulize))
    }
    
    #[test]
    fn propenyl_radical_given_aromatic_bonds() {
        let atoms = read("C:C:C").unwrap();
        let atoms = kekulize(atoms);

        assert_eq!(atoms, Err(Error::CanNotKekulize))
    }

    #[test]
    fn ethane() {
        let atoms = read("CC").unwrap();
        let atoms = kekulize(atoms).unwrap();

        assert_eq!(atoms, vec![
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 1, style: None }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 0, style: None }
                ]
            }
        ]);
    }

    #[test]
    fn ethene_with_aromatic_atoms() {
        let atoms = read("cc").unwrap();
        let atoms = kekulize(atoms).unwrap();

        assert_eq!(atoms, vec![
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 1, style: Some(Style::Double) }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 0, style: Some(Style::Double) }
                ]
            }
        ]);
    }

    #[test]
    fn ethene_with_aromatic_bond() {
        let atoms = read("C:C").unwrap();
        let atoms = kekulize(atoms).unwrap();

        assert_eq!(atoms, vec![
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 1, style: Some(Style::Double) }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 0, style: Some(Style::Double) }
                ]
            }
        ]);
    }

    #[test]
    fn ethene_with_aromatic_atoms_and_bond() {
        let atoms = read("c:c").unwrap();
        let atoms = kekulize(atoms).unwrap();

        assert_eq!(atoms, vec![
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 1, style: Some(Style::Double) }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 0, style: Some(Style::Double) }
                ]
            }
        ]);
    }

    #[test]
    fn butadiene_with_aromatic_bonds() {
        let atoms = read("C:C:C:C").unwrap();
        let atoms = kekulize(atoms).unwrap();
        
        assert_eq!(atoms, vec![
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 1, style: Some(Style::Double) }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 0, style: Some(Style::Double) },
                    Bond { tid: 2, style: Some(Style::Single) }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 1, style: Some(Style::Single) },
                    Bond { tid: 3, style: Some(Style::Double) }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 2, style: Some(Style::Double) }
                ]
            }
        ]);
    }

    #[test]
    fn furan_all_aromatic_bonds() {
        let atoms = read("O1:C:C:C:C:1").unwrap();
        let atoms = kekulize(atoms).unwrap();
        
        assert_eq!(atoms, vec![
            Atom {
                nub: Nub { element: Element::O, ..Default::default() },
                bonds: vec![
                    Bond { tid: 4, style: Some(Style::Single) },
                    Bond { tid: 1, style: Some(Style::Single) }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 0, style: Some(Style::Single) },
                    Bond { tid: 2, style: Some(Style::Double) }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 1, style: Some(Style::Double) },
                    Bond { tid: 3, style: Some(Style::Single) }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 2, style: Some(Style::Single) },
                    Bond { tid: 4, style: Some(Style::Double) }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 3, style: Some(Style::Double) },
                    Bond { tid: 0, style: Some(Style::Single) }
                ]
            }
        ]);
    }

    #[test]
    fn thiazole_all_aromatic_atoms() {
        let atoms = read("c1cscn1").unwrap();

        assert_eq!(kekulize(atoms).unwrap(), vec![
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 4, style: None },
                    Bond { tid: 1, style: Some(Style::Double) },
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 0, style: Some(Style::Double) },
                    Bond { tid: 2, style: None }
                ]
            },
            Atom {
                nub: Nub { element: Element::S, ..Default::default() },
                bonds: vec![
                    Bond { tid: 1, style: None },
                    Bond { tid: 3, style: None }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 2, style: None },
                    Bond { tid: 4, style: Some(Style::Double) }
                ]
            },
            Atom {
                nub: Nub { element: Element::N, ..Default::default() },
                bonds: vec![
                    Bond { tid: 3, style: Some(Style::Double) },
                    Bond { tid: 0, style: None}
                ]
            }
        ])
    }

    #[test]
    fn furan_all_aromatic_atoms() {
        let atoms = read("o1cccc1").unwrap();
        let atoms = kekulize(atoms).unwrap();

        assert_eq!(atoms, vec![
            Atom {
                nub: Nub { element: Element::O, ..Default::default() },
                bonds: vec![
                    Bond { tid: 4, style: None },
                    Bond { tid: 1, style: None }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 0, style: None },
                    Bond { tid: 2, style: Some(Style::Double) }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 1, style: Some(Style::Double) },
                    Bond { tid: 3, style: None }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 2, style: None },
                    Bond { tid: 4, style: Some(Style::Double) }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 3, style: Some(Style::Double) },
                    Bond { tid: 0, style: None }
                ]
            }
        ]);
    }

    #[test]
    fn selenophene_all_aromatic_atoms() {
        let atoms = read("[se]1cccc1").unwrap();
        let atoms = kekulize(atoms).unwrap();

        assert_eq!(atoms, vec![
            Atom {
                nub: Nub {
                    element: Element::Se,
                    hcount: Some(0),
                    charge: Some(0),
                    ..Default::default()
                },
                bonds: vec![
                    Bond { tid: 4, style: None },
                    Bond { tid: 1, style: None }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 0, style: None },
                    Bond { tid: 2, style: Some(Style::Double) }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 1, style: Some(Style::Double) },
                    Bond { tid: 3, style: None }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 2, style: None },
                    Bond { tid: 4, style: Some(Style::Double) }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 3, style: Some(Style::Double) },
                    Bond { tid: 0, style: None }
                ]
            }
        ]);
    }

    #[test]
    fn hydrogen_selenidine_all_aromatic_atoms() {
        let atoms = read("[seH]1ccccc1").unwrap();
        let atoms = kekulize(atoms).unwrap();

        assert_eq!(atoms[0], Atom {
            nub: Nub {
                element: Element::Se,
                hcount: Some(1),
                charge: Some(0),
                ..Default::default()
            },
            bonds: vec![
                Bond { tid: 5, style: Some(Style::Double) },
                Bond { tid: 1, style: None },
            ]
        })
    }

    #[test]
    fn ketene_with_aromatic_carbons() {
        let atoms = read("cc=O").unwrap();
        let atoms = kekulize(atoms).unwrap();

        assert_eq!(atoms, vec![
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 1, style: Some(Style::Double) }
                ]
            },
            Atom {
                nub: Default::default(),
                bonds: vec![
                    Bond { tid: 0, style: Some(Style::Double) },
                    Bond { tid: 2, style: Some(Style::Double) }
                ]
            },
            Atom {
                nub: Nub { element: Element::O, ..Default::default() },
                bonds: vec![
                    Bond { tid: 1, style: Some(Style::Double) }
                ]
            }
        ]);
    }
}