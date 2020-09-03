use purr::mol::{ Atom, Style };
use gamma::matching::greedy;
use gamma::graph::Graph;
use crate::molecule::Error;
use super::pi_subgraph;

pub fn kekulize(mut atoms: Vec<Atom>) -> Result<Vec<Atom>, Error> {
    let pi = pi_subgraph(&atoms)?;

    println!("{:?}", pi);
    
    if pi.is_empty() {
        return Ok(atoms);
    }
    
    let matching = greedy(&pi);

    println!(" matching {:?}", matching);
    
    if matching.len() * 2 != pi.order() {
        return Err(Error::CanNotKekulize);
    }

    for (sid, tid) in matching {
        let source = &mut atoms[sid];

        if let Some(bond) = source.bonds.iter_mut().find(|bond| bond.tid == tid) {
            bond.style = Some(Style::Double);
        } else {

        }

        let target = &mut atoms[tid];

        if let Some(bond) = target.bonds.iter_mut().find(|bond| bond.tid == sid) {
            bond.style = Some(Style::Double);
        }
    }

    for atom in atoms.iter_mut() {
        atom.nub.aromatic = false;
    }

    for atom in atoms.iter_mut() {
        for out in atom.bonds.iter_mut() {
            if let Some(Style::Aromatic) = out.style {
                out.style.replace(Style::Single);
            }
        }
    }

    Ok(atoms)
}

#[cfg(test)]
mod tests {
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
        println!("{:#?}", atoms);
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
    fn foo() {
        let atoms = read("c1cscn1").unwrap();
        let atoms = kekulize(atoms).unwrap();
    }

    #[test]
    fn bar() {
        let atoms = read("c1nc2c([nH]1)c(=O)[nH]cn2").unwrap();
        let atoms = kekulize(atoms).unwrap();
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
}