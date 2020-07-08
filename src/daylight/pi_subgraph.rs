use purr::mol::{ Atom, Bond, Style };
use purr::valence::{ hypovalence, Error as ValenceError };
use gamma::graph::HashGraph;
use crate::molecule::Error;

pub fn pi_subgraph(atoms: &[Atom]) -> Result<HashGraph, Error> {
    let mut edges = Vec::new();
    let mut singletons = Vec::new();

    for (sid, atom) in atoms.iter().enumerate() {
        if atom.nub.aromatic && atom.bonds.is_empty() {
            singletons.push(sid);
        }
    }

    for (sid, atom) in atoms.iter().enumerate() {
        for Bond { tid, style } in atom.bonds.iter() {
            if sid < *tid && is_aromatic_bond(sid, *tid, style, atoms)? {
                edges.push((sid, *tid));
            }
        }
    }

    Ok(HashGraph::from_edges(edges, singletons).expect(
        "error creating pi subgraph from adjacency"
    ))
}

fn is_aromatic_bond(
    sid: usize, tid: usize, style: &Option<Style>, atoms: &[Atom]
) -> Result<bool, Error> {
    let source = &atoms[sid];
    let target = &atoms[tid];

    match style {
        Some(Style::Aromatic) =>
            Ok(subvalent(sid, atoms)? && subvalent(tid, atoms)?),
        None => {
            if source.nub.aromatic && target.nub.aromatic {
                Ok(subvalent(sid, atoms)? && subvalent(tid, atoms)?)
            } else {
                Ok(false)
            }
        },
        _ => Ok(false)
    }
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
    use super::*;
    use purr::read::read;

    #[test]
    fn aromatic_methane() {
        let pi = pi_subgraph(&read("c").unwrap());

        assert_eq!(pi, Ok(HashGraph::from_edges(vec![ ], vec![ 0 ]).unwrap()));
    }

    #[test]
    fn methane() {
        let pi = pi_subgraph(&read("C").unwrap());

        assert_eq!(pi, Ok(HashGraph::from_edges(vec![ ], vec![ ]).unwrap()));
    }

    #[test]
    fn ethane() {
        let pi = pi_subgraph(&read("CC").unwrap());

        assert_eq!(pi, Ok(HashGraph::from_edges(vec![ ], vec![ ]).unwrap()));
    }
    
    #[test]
    fn ethyne() {
        let pi = pi_subgraph(&read("C#C").unwrap());

        assert_eq!(pi, Ok(HashGraph::from_edges(vec![ ], vec![ ]).unwrap()));
    }

    #[test]
    fn ethene_with_aromatic_atoms() {
        let pi = pi_subgraph(&read("cc").unwrap());

        assert_eq!(pi, Ok(HashGraph::from_edges(vec![
            (0, 1)
        ], vec![ ]).unwrap()));
    }

    #[test]
    fn ethene_with_aromatic_bond() {
        let pi = pi_subgraph(&read("C:C").unwrap());

        assert_eq!(pi, Ok(HashGraph::from_edges(vec![
            (0, 1)
        ], vec![ ]).unwrap()));
    }

    #[test]
    fn propenyl_radical_with_aromatic_atoms() {
        let pi = pi_subgraph(&read("ccc").unwrap());

        assert_eq!(pi, Ok(HashGraph::from_edges(vec![
            (0, 1),
            (1, 2)
        ], vec![ ]).unwrap()));
    }

    #[test]
    fn pyrrole_with_aromatic_atoms() {
        let pi = pi_subgraph(&read("[nH]1cccc1").unwrap());

        assert_eq!(pi, Ok(HashGraph::from_edges(vec![
            (1, 2),
            (2, 3),
            (3, 4)
        ], vec![ ]).unwrap()));
    }

    #[test]
    fn furan_with_aromatic_atoms() {
        let pi = pi_subgraph(&read("o1cccc1").unwrap());

        assert_eq!(pi, Ok(HashGraph::from_edges(vec![
            (1, 2),
            (2, 3),
            (3, 4)
        ], vec![ ]).unwrap()));
    }

    #[test]
    fn pyridine_with_aromatic_atoms() {
        let pi = pi_subgraph(&read("c1ccccc1").unwrap());

        assert_eq!(pi, Ok(HashGraph::from_edges(vec![
            (0, 5),
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5)
        ], vec![ ]).unwrap()));
    }

    #[test]
    fn pyridine_with_aromatic_bonds() {
        let pi = pi_subgraph(&read("C1:C:C:C:C:C:1").unwrap());

        assert_eq!(pi, Ok(HashGraph::from_edges(vec![
            (0, 5),
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5)
        ], vec![ ]).unwrap()));
    }

    #[test]
    fn pyridine_with_bracket_nitrogen() {
        let pi = pi_subgraph(&read("[n]1ccccc1").unwrap());

        assert_eq!(pi, Ok(HashGraph::from_edges(vec![
            (0, 5),
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5)
        ], vec![ ]).unwrap()));
    }

    #[test]
    fn pyridine_without_bracket_nitrogen() {
        let pi = pi_subgraph(&read("n1ccccc1").unwrap());

        assert_eq!(pi, Ok(HashGraph::from_edges(vec![
            (0, 5),
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5)
        ], vec![ ]).unwrap()));
    }
}