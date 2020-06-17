use purr::mol::{ Mol, Atom, Style };
use purr::valence::hypovalence;
use gamma::graph::StableGraph;

pub fn pi_subgraph(mol: &Mol) -> StableGraph<usize, ()> {
    let mut flags = vec![false; mol.atoms.len()];
    let mut edges = Vec::new();

    for (id, atom) in mol.atoms.iter().enumerate() {
        if eligible(id, &atom, &mol) {
            std::mem::replace(&mut flags[id], true);
        }
    }

    for (sid, _) in mol.atoms.iter().enumerate() {
        if !flags.get(sid).unwrap() {
            continue;
        }

        for bond in mol.bonds.get(sid).unwrap().iter() {
            let tid = bond.tid;

            if tid > sid {
                match bond.style {
                    Some(Style::Aromatic) => {
                        edges.push((sid, tid, ()));
                    },
                    None => {
                        edges.push((sid, tid, ()));
                    },
                    Some(_) => ()
                }
            }
        }
    }

    let nodes = (0..mol.atoms.len()).filter(
        |id| *flags.get(*id).unwrap()
    ).collect::<Vec<_>>();
    
    StableGraph::build(nodes, edges).unwrap()
}

fn eligible(id: usize, atom: &Atom, mol: &Mol) -> bool {
    let bonds = &mol.bonds[id];

    if !atom.aromatic {
        if !bonds.iter().any(|bond| match bond.style {
            Some(Style::Aromatic) => true,
            _ => false
        }) {
            return false;
        }
    }

    match hypovalence(atom, bonds) {
        None => false,
        Some(value) => value > 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use purr::read::read;
    use gamma::graph::Graph;

    #[test]
    fn methane() {
        let pi = pi_subgraph(&read("C").unwrap());

        assert_eq!(pi.edges().collect::<Vec<_>>(), vec![ ])
    }

    #[test]
    fn ethane() {
        let pi = pi_subgraph(&read("CC").unwrap());

        assert_eq!(pi.edges().collect::<Vec<_>>(), vec![ ]);
    }
    
    #[test]
    fn ethyne() {
        let pi = pi_subgraph(&read("C#C").unwrap());

        assert_eq!(pi.edges().collect::<Vec<_>>(), vec![ ]);
    }

    #[test]
    fn ethene_with_aromatic_atoms() {
        let pi = pi_subgraph(&read("cc").unwrap());

        assert_eq!(pi.edges().collect::<Vec<_>>(), vec![
            (&0, &1)
        ]);
    }

    #[test]
    fn ethene_with_aromatic_bond() {
        let pi = pi_subgraph(&read("C:C").unwrap());

        assert_eq!(pi.edges().collect::<Vec<_>>(), vec![
            (&0, &1)
        ]);
    }

    #[test]
    fn pyrrole_with_aromatic_atoms() {
        let pi = pi_subgraph(&read("[nH]1cccc1").unwrap());

        assert_eq!(pi.edges().collect::<Vec<_>>(), vec![
            (&1, &2),
            (&2, &3),
            (&3, &4),
        ]);
    }

    #[test]
    fn furan_with_aromatic_atoms() {
        let pi = pi_subgraph(&read("o1cccc1").unwrap());

        assert_eq!(pi.edges().collect::<Vec<_>>(), vec![
            (&1, &2),
            (&2, &3),
            (&3, &4),
        ]);
    }

    #[test]
    fn pyridine_with_aromatic_atoms() {
        let pi = pi_subgraph(&read("c1ccccc1").unwrap());

        assert_eq!(pi.edges().collect::<Vec<_>>(), vec![
            (&0, &5),
            (&0, &1),
            (&1, &2),
            (&2, &3),
            (&3, &4),
            (&4, &5)
        ]);
    }

    #[test]
    fn pyridine_with_aromatic_bonds() {
        let pi = pi_subgraph(&read("C1:C:C:C:C:C:1").unwrap());

        assert_eq!(pi.edges().collect::<Vec<_>>(), vec![
            (&0, &5),
            (&0, &1),
            (&1, &2),
            (&2, &3),
            (&3, &4),
            (&4, &5)
        ]);
    }

    #[test]
    fn pyridine_with_bracket_nitrogen() {
        let pi = pi_subgraph(&read("[n]1ccccc1").unwrap());

        assert_eq!(pi.edges().collect::<Vec<_>>(), vec![
            (&0, &5),
            (&0, &1),
            (&1, &2),
            (&2, &3),
            (&3, &4),
            (&4, &5)
        ]);
    }

    #[test]
    fn pyridine_without_bracket_nitrogen() {
        let pi = pi_subgraph(&read("n1ccccc1").unwrap());

        assert_eq!(pi.edges().collect::<Vec<_>>(), vec![
            (&0, &5),
            (&0, &1),
            (&1, &2),
            (&2, &3),
            (&3, &4),
            (&4, &5)
        ]);
    }
}