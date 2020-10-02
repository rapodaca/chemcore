use purr::mol::{ Atom, Bond, Style };
use purr::valence::{ hypovalence, Error as ValenceError };
use gamma::graph::{ Graph, DefaultGraph };
use crate::molecule::Error;

pub fn pi_subgraph(atoms: &[Atom]) -> Result<DefaultGraph, Error> {
    let mut result = DefaultGraph::new();

    for (id, atom) in atoms.iter().enumerate() {
        if atom.nub.aromatic {
            if subvalent(id, atoms)? {
                result.add_node(id).unwrap();
            }
        } else {
            let aromatic_bond = atom.bonds.iter().any(|bond| match bond.style {
                Some(Style::Aromatic) => true,
                _ => false
            });

            if aromatic_bond && subvalent(id, atoms)? {
                result.add_node(id).unwrap();
            }
        }
    }

    for (sid, atom) in atoms.iter().enumerate() {
        if !result.has_node(sid) {
            continue;
        }

        for Bond { tid, style } in atom.bonds.iter() {
            if sid > *tid {
                continue;
            }

            if result.has_node(*tid) {
                match style {
                    Some(Style::Aromatic) | None => {
                        result.add_edge(sid, *tid).unwrap();
                    },
                    _ => ()
                }
            }
        }
    } 

    Ok(result)
}

fn subvalent(id: usize, atoms: &[Atom]) -> Result<bool, Error> {
    let atom = &atoms[id];

    match hypovalence(atom) {
        Ok(result) => match result {
            Some(hypovalence) => Ok(hypovalence > 0),
            None => Ok(false)
        },
        Err(ValenceError::UnmatchableValence) =>
            return Err(Error::Hypervalent(id))
    }
}

#[cfg(test)]
mod tests {
    use std::convert::TryFrom;
    use super::*;
    use purr::read::read;

    #[test]
    fn methane() {
        let pi = pi_subgraph(&read("C").unwrap());

        assert_eq!(pi, Ok(DefaultGraph::new()))
    }

    #[test]
    fn aromatic_methane() {
        let pi = pi_subgraph(&read("c").unwrap());
        let mut expected = DefaultGraph::new();

        expected.add_node(0).unwrap();

        assert_eq!(pi, Ok(expected))
    }

    #[test]
    fn ethane() {
        let pi = pi_subgraph(&read("CC").unwrap());

        assert_eq!(pi, Ok(DefaultGraph::new()))
    }
    
    #[test]
    fn ethyne() {
        let pi = pi_subgraph(&read("C#C").unwrap());

        assert_eq!(pi, Ok(DefaultGraph::new()))
    }

    #[test]
    fn ethene_with_aromatic_atoms() {
        let pi = pi_subgraph(&read("cc").unwrap());

        assert_eq!(pi, Ok(DefaultGraph::try_from(vec![
            (0, 1)
        ]).unwrap()))
    }

    #[test]
    fn ethene_with_aromatic_atoms_and_double_bond() {
        let pi = pi_subgraph(&read("c=c").unwrap());

        assert_eq!(pi, Ok(DefaultGraph::try_from(vec![
            vec![ ],
            vec![ ]
        ]).unwrap()))
    }

    #[test]
    fn ethene_with_aromatic_atoms_and_triple_bond() {
        let pi = pi_subgraph(&read("c#c").unwrap());

        assert_eq!(pi, Ok(DefaultGraph::try_from(vec![
            vec![ ],
            vec![ ]
        ]).unwrap()))
    }

    #[test]
    fn ethene_with_aromatic_bond() {
        let pi = pi_subgraph(&read("C:C").unwrap());

        assert_eq!(pi, Ok(DefaultGraph::try_from(vec![
            (0, 1)
        ]).unwrap()))
    }

    #[test]
    fn propyl_radical_with_aromatic_atom() {
        let pi = pi_subgraph(&read("CCc").unwrap());
        let mut expected = DefaultGraph::new();

        expected.add_node(2).unwrap();

        assert_eq!(pi, Ok(expected))
    }

    #[test]
    fn propenyl_radical_with_aromatic_atoms() {
        let pi = pi_subgraph(&read("ccc").unwrap());

        assert_eq!(pi, Ok(DefaultGraph::try_from(vec![
            (0, 1),
            (1, 2)
        ]).unwrap()))
    }

    #[test]
    fn pyrrole_with_aromatic_atoms() {
        let pi = pi_subgraph(&read("[nH]1cccc1").unwrap());

        assert_eq!(pi, Ok(DefaultGraph::try_from(vec![
            (1, 2),
            (2, 3),
            (3, 4)
        ]).unwrap()))
    }

    #[test]
    fn thiazole_with_aromatic_atoms() {
        let pi = pi_subgraph(&read("c1cscn1").unwrap());

        assert_eq!(pi, Ok(DefaultGraph::try_from(vec![
            (0, 1),
            (0, 4),
            (3, 4)
        ]).unwrap()))
    }

    #[test]
    fn furan_with_aromatic_atoms() {
        let pi = pi_subgraph(&read("o1cccc1").unwrap());

        assert_eq!(pi, Ok(DefaultGraph::try_from(vec![
            (1, 2),
            (2, 3),
            (3, 4)
        ]).unwrap()))
    }

    #[test]
    fn pyridine_with_aromatic_atoms() {
        let pi = pi_subgraph(&read("c1ccccc1").unwrap());

        assert_eq!(pi, Ok(DefaultGraph::try_from(vec![
            (0, 5),
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5)
        ]).unwrap()))
    }

    #[test]
    fn pyridine_with_aromatic_bonds() {
        let pi = pi_subgraph(&read("C1:C:C:C:C:C:1").unwrap());

        assert_eq!(pi, Ok(DefaultGraph::try_from(vec![
            (0, 5),
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5)
        ]).unwrap()));
    }

    #[test]
    fn pyridine_with_bracket_nitrogen() {
        let pi = pi_subgraph(&read("[n]1ccccc1").unwrap());

        assert_eq!(pi, Ok(DefaultGraph::try_from(vec![
            (0, 5),
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5)
        ]).unwrap()));
    }

    #[test]
    fn pyridine_without_bracket_nitrogen() {
        let pi = pi_subgraph(&read("n1ccccc1").unwrap());

        assert_eq!(pi, Ok(DefaultGraph::try_from(vec![
            (0, 5),
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5)
        ]).unwrap()));
    }
}