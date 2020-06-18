use gamma::graph::Graph;
use gamma::matching::greedy;
use purr::mol::{ Mol, Atom, Bond, Style };
use super::pi_subgraph;

pub fn kekulize(mol: Mol) -> Mol {
    let pi = pi_subgraph(&mol);
    let matching = greedy(&pi);

    if matching.is_empty() {
        return mol;
    }

    let mut result = Mol { atoms: Vec::new(), bonds: Vec::new() };

    for atom in mol.atoms {
        if atom.aromatic {
            result.atoms.push(Atom { aromatic: false, ..atom });
        } else {
            result.atoms.push(atom);
        }
    }

    for (sid, outs) in mol.bonds.iter().enumerate() {
        result.bonds.push(outs.iter().map(|bond| {
            let tid = bond.tid;

            if matching.has_node(&sid) && matching.has_node(&tid) {
                if matching.has_edge(&sid, &tid).unwrap() {
                    Bond { tid, style: Some(Style::Double) }
                } else {
                    Bond { tid, style: bond.style }
                }
            } else {
                Bond { tid, style: bond.style }
            }
        }).collect::<Vec<_>>());
    }
    
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use purr::mol::Style;
    use purr::read::{ read, Error };

    #[test]
    fn aromatic_propane() {
        let mol = kekulize(read("C:C:C").unwrap());

        assert_eq!(mol, Mol {
            atoms: vec![
                Atom { ..Default::default() },
                Atom { ..Default::default() },
                Atom { ..Default::default() },
            ],
            bonds: vec![ 
                vec![
                    Bond { tid: 1, style: Some(Style::Double) }
                ],
                vec![
                    Bond { tid: 0, style: Some(Style::Double) },
                    Bond { tid: 2, style: Some(Style::Aromatic) }
                ],
                vec![
                    Bond { tid: 1, style: Some(Style::Aromatic) }
                ]
            ]
        })
    }

    #[test]
    fn methane() {
        let mol = kekulize(read("C").unwrap());

        assert_eq!(mol, Mol {
            atoms: vec![
                Atom { ..Default::default() }
            ],
            bonds: vec![
                vec![ ]
            ]
        });
    }

    #[test]
    fn ethane() {
        let mol = kekulize(read("CC").unwrap());

        assert_eq!(mol, Mol {
            atoms: vec![
                Atom { ..Default::default() },
                Atom { ..Default::default() }
            ],
            bonds: vec![
                vec![ Bond { tid: 1, style: None } ],
                vec![ Bond { tid: 0, style: None } ]
            ]
        });
    }

    #[test]
    fn ethane_with_aromatic_atoms() {
        let mol = kekulize(read("cc").unwrap());

        assert_eq!(mol, Mol {
            atoms: vec![
                Atom { ..Default::default() },
                Atom { ..Default::default() }
            ],
            bonds: vec![
                vec![ Bond { tid: 1, style: Some(Style::Double) } ],
                vec![ Bond { tid: 0, style: Some(Style::Double) } ]
            ]
        });
    }

    #[test]
    fn benzene_with_aromatic_atoms() {
        let mol = kekulize(read("c1ccccc1").unwrap());

        assert_eq!(mol, Mol {
            atoms: vec![
                Atom { ..Default::default() },
                Atom { ..Default::default() },
                Atom { ..Default::default() },
                Atom { ..Default::default() },
                Atom { ..Default::default() },
                Atom { ..Default::default() },
            ],
            bonds: vec![
                vec![
                    Bond { tid: 5, style: Some(Style::Double) },
                    Bond { tid: 1, style: None }
                ],
                vec![
                    Bond { tid: 0, style: None },
                    Bond { tid: 2, style: Some(Style::Double) }
                ],
                vec![
                    Bond { tid: 1, style: Some(Style::Double) },
                    Bond { tid: 3, style: None }
                ],
                vec![
                    Bond { tid: 2, style: None },
                    Bond { tid: 4, style: Some(Style::Double) }
                ],
                vec![
                    Bond { tid: 3, style: Some(Style::Double) },
                    Bond { tid: 5, style: None }
                ],
                vec![
                    Bond { tid: 4, style: None },
                    Bond { tid: 0, style: Some(Style::Double) }
                ]
            ]
        });
    }
}