use purr::graph::{ Atom, Bond };
use purr::parts::BondKind;
use gamma::graph::{ DefaultGraph, Graph };

pub fn pi_subgraph(atoms: &Vec<Atom>) -> DefaultGraph {
    let mut result = DefaultGraph::new();
    let mut subvalences = vec![ ];
    
    for (index, atom) in atoms.iter().enumerate() {
        let subvalence = atom.subvalence();

        if atom.is_aromatic() && subvalence > 0 {
            result.add_node(index).expect("add node");
        }

        subvalences.push(subvalence);
    }

    for (sid, source) in atoms.iter().enumerate() {
        for Bond { tid, kind } in source.bonds.iter() {
            if *tid < sid {
                continue
            }

            match kind {
                BondKind::Elided => {
                    if result.has_id(sid) && result.has_id(*tid) {
                        result.add_edge(sid, *tid).expect("add edge")
                    }
                },
                BondKind::Aromatic => {
                    if subvalences[sid] > 0 {
                        if !result.has_id(sid) {
                            result.add_node(sid).expect("add source");
                        }

                        if subvalences[*tid] > 0 {
                            if !result.has_id(*tid) {
                                result.add_node(*tid).expect("add target");
                            }

                            result.add_edge(sid, *tid).expect("add edge")
                        }
                    } else if subvalences[*tid] > 0 {
                        if !result.has_id(*tid) {
                            result.add_node(*tid).expect("add target");
                        }
                    }
                },
                _ => ()
            }
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use std::convert::TryFrom;
    use pretty_assertions::assert_eq;
    use purr::read::read;
    use purr::graph::from_tree;
    use super::*;

    #[test]
    fn methane() {
        let atoms = from_tree(read("C").unwrap().root).unwrap();

        assert_eq!(pi_subgraph(&atoms), DefaultGraph::new())
    }

    #[test]
    fn methane_aromatic() {
        let atoms = from_tree(read("c").unwrap().root).unwrap();

        assert_eq!(pi_subgraph(&atoms), DefaultGraph::try_from(vec![
            vec![ ]
        ]).unwrap())
    }

    #[test]
    fn ethane() {
        let atoms = from_tree(read("CC").unwrap().root).unwrap();

        assert_eq!(pi_subgraph(&atoms), DefaultGraph::new())
    }

    #[test]
    fn ethene_aromatic_atoms() {
        let atoms = from_tree(read("cc").unwrap().root).unwrap();

        assert_eq!(pi_subgraph(&atoms), DefaultGraph::try_from(vec![
            vec![ 1 ],
            vec![ 0 ]
        ]).unwrap())
    }

    #[test]
    fn propene_aromatic_atoms() {
        let atoms = from_tree(read("ccC").unwrap().root).unwrap();

        assert_eq!(pi_subgraph(&atoms), DefaultGraph::try_from(vec![
            vec![ 1 ],
            vec![ 0 ]
        ]).unwrap())
    }

    #[test]
    fn carbon_iron_aromatic() {
        let atoms = from_tree(read("C:[Fe]").unwrap().root).unwrap();

        assert_eq!(pi_subgraph(&atoms), DefaultGraph::try_from(vec![
            vec![ ]
        ]).unwrap())
    }

    #[test]
    fn iron_carbon_aromatic() {
        let atoms = from_tree(read("[Fe]:C").unwrap().root).unwrap();
        let mut result = DefaultGraph::new();

        result.add_node(1).unwrap();

        assert_eq!(pi_subgraph(&atoms), result)
    }

    #[test]
    fn furan_all_aromatic() {
        let atoms = from_tree(read("c1ccco1").unwrap().root).unwrap();

        assert_eq!(pi_subgraph(&atoms), DefaultGraph::try_from(vec![
            (0, 1),
            (1, 2),
            (2, 3)
        ]).unwrap())
    }

    #[test]
    fn pyrrole_all_aromatic() {
        let atoms = from_tree(read("c1ccc[nH]1").unwrap().root).unwrap();

        assert_eq!(pi_subgraph(&atoms), DefaultGraph::try_from(vec![
            (0, 1),
            (1, 2),
            (2, 3)
        ]).unwrap())
    }

    #[test]
    fn benzene_all_aromatic() {
        let atoms = from_tree(read("c1ccccc1").unwrap().root).unwrap();

        assert_eq!(pi_subgraph(&atoms), DefaultGraph::try_from(vec![
            vec![ 5, 1 ],
            vec![ 0, 2 ],
            vec![ 1, 3 ],
            vec![ 2, 4 ],
            vec![ 3, 5 ],
            vec![ 0, 4 ]
        ]).unwrap())
    }
}