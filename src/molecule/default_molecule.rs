use gamma::graph::{ Graph, Error as GraphError };
use super::{ Node, Molecule, Atom };

#[derive(Debug,PartialEq)]
pub struct DefaultMolecule {
    nodes: Vec<Node>,
    size: usize
}

impl DefaultMolecule {
    pub fn new(nodes: Vec<Node>) -> Self {
        let size = nodes.iter().fold(0, |sum, node| sum + node.bonds.len());

        assert!(size % 2 == 0, "odd bond count");

        DefaultMolecule {
            nodes,
            size: size / 2
        }
    }

    fn node_for(&self, id: usize) -> Result<&Node, GraphError> {
        match self.nodes.get(id) {
            Some(node) => Ok(node),
            None => Err(GraphError::UnknownId(id))
        }
    }
}

impl Graph for DefaultMolecule {
    fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    fn order(&self) -> usize {
        self.nodes.len()
    }

    fn size(&self) -> usize {
        self.size
    }

    fn ids(&self) -> Box<dyn Iterator<Item=usize> + '_> {
        Box::new(0..self.nodes.len())
    }

    fn neighbors(
        &self, id: usize
    ) -> Result<Box<dyn Iterator<Item=usize> + '_>, GraphError> {
        let node = self.node_for(id)?;

        Ok(Box::new(node.bonds.iter().map(|bond| bond.tid)))
    }

    fn has_id(&self, id: usize) -> bool {
        id < self.nodes.len()
    }

    fn degree(&self, id: usize) -> Result<usize, GraphError> {
        Ok(self.node_for(id)?.bonds.len())
    }

    fn edges(&self) -> Box<dyn Iterator<Item=(usize, usize)> + '_> {
        // let mut result = Vec::new();

        // for (sid, node) in self.nodes.iter().enumerate() {
        //     for edge in node.bonds.iter() {
        //         if sid < edge.tid {
        //             result.push((sid, edge.tid))
        //         }
        //     }
        // }

        // result
        Box::new(EdgeIterator::new(&self.nodes))
    }

    fn has_edge(&self, sid: usize, tid: usize) -> Result<bool, GraphError> {
        let source = &self.node_for(sid)?;

        if tid < self.nodes.len() {
            Ok(source.bonds.iter().any(|bond| bond.tid == tid))
        } else {
            Err(GraphError::UnknownId(tid))
        }
    }
}

impl Molecule for DefaultMolecule {
    fn atom(&self, id: usize) -> Result<&Atom, GraphError> {
        Ok(&self.node_for(id)?.atom)
    }

    fn charge(&self, id: usize) -> Result<f32, GraphError> {
        let node = self.node_for(id)?;
        let element = match &node.atom.element {
            Some(element) => element,
            None => return Ok(0f32)
        };

        let mut result = element.valence_electrons() as f32;

        result -= node.atom.hydrogens as f32;
        result -= node.bonds.iter().fold(0f32, |sum, bond| {
            sum + bond.order()
        });
        result -= node.atom.electrons as f32;

        Ok(result)
    }

    fn bond_order(&self, sid: usize, tid: usize) -> Result<f32, GraphError> {
        let source = self.node_for(sid)?;

        if tid >= self.nodes.len() {
            return Err(GraphError::UnknownId(tid))
        }

        match source.bonds.iter().find(|bond| bond.tid == tid) {
            Some(bond) => Ok(bond.order()),
            None => Ok(0f32)
        }
    }
}

struct EdgeIterator<'a> {
    nodes: &'a Vec<Node>,
    row: usize,
    col: usize
}

impl<'a> EdgeIterator<'a> {
    fn new(nodes: &'a Vec<Node>) -> Self {
        Self {
            nodes,
            row: 0,
            col: 0
        }
    }
}

impl<'a> Iterator for EdgeIterator<'a> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.nodes.get(self.row) {
                Some(node) => match node.bonds.get(self.col) {
                    Some(bond) => {
                        self.col += 1;

                        if bond.tid > self.row {
                            break Some((self.row, bond.tid))
                        }
                    },
                    None => {
                        self.col = 0;
                        self.row += 1
                    }
                },
                None => break None
            }
        }
    }
}

#[cfg(test)]
mod is_empty {
    use pretty_assertions::assert_eq;
    use super::*;

    #[test]
    fn empty() {
        let molecule = DefaultMolecule::new(vec![ ]);

        assert_eq!(molecule.is_empty(), true)
    }

    #[test]
    fn one_atom() {
        let molecule = DefaultMolecule::new(vec![
            Node { atom: Atom::default(), bonds: vec![ ] }
        ]);

        assert_eq!(molecule.is_empty(), false)
    }
}

#[cfg(test)]
mod order {
    use pretty_assertions::assert_eq;
    use super::*;

    #[test]
    fn empty() {
        let molecule = DefaultMolecule::new(vec![ ]);

        assert_eq!(molecule.order(), 0)
    }

    #[test]
    fn one_atom() {
        let molecule = DefaultMolecule::new(vec![
            Node { atom: Atom::default(), bonds: vec![ ] }
        ]);

        assert_eq!(molecule.order(), 1)
    }
}

#[cfg(test)]
mod size {
    use pretty_assertions::assert_eq;
    use super::super::{ Bond };
    use super::*;

    #[test]
    fn empty() {
        let molecule = DefaultMolecule::new(vec![ ]);

        assert_eq!(molecule.size(), 0)
    }

    #[test]
    fn one_edge() {
        let molecule = DefaultMolecule::new(vec![
            Node { atom: Atom::default(), bonds: vec![
                Bond::new(2, None, 1)
            ]},
            Node { atom: Atom::default(), bonds: vec![
                Bond::new(2, None, 0)
            ]}
        ]);

        assert_eq!(molecule.size(), 1)
    }
}

#[cfg(test)]
mod nodes {
    use pretty_assertions::assert_eq;
    use super::*;

    #[test]
    fn empty() {
        let molecule = DefaultMolecule::new(vec![ ]);

        assert_eq!(molecule.ids().collect::<Vec<_>>(), [ ])
    }

    #[test]
    fn three_atoms() {
        let molecule = DefaultMolecule::new(vec![
            Node { atom: Atom::default(), bonds: vec![ ] },
            Node { atom: Atom::default(), bonds: vec![ ] },
            Node { atom: Atom::default(), bonds: vec![ ] }
        ]);

        assert_eq!(molecule.ids().collect::<Vec<_>>(), [ 0, 1, 2 ])
    }
}

#[cfg(test)]
mod neighbors {
    use pretty_assertions::assert_eq;
    use super::super::{ Bond };
    use super::*;

    #[test]
    fn unknown_id() {
        let molecule = DefaultMolecule::new(vec![ ]);

        assert_eq!(molecule.neighbors(0).err(), Some(GraphError::UnknownId(0)))
    }

    #[test]
    fn known_id() {
        let molecule = DefaultMolecule::new(vec![
            Node { atom: Atom::default(), bonds: vec![
                Bond::new(2, None, 1),
                Bond::new(2, None, 2)
            ]},
            Node { atom: Atom::default(), bonds: vec![
                Bond::new(2, None, 0)
            ]},
            Node { atom: Atom::default(), bonds: vec![
                Bond::new(2, None, 0)
            ]}
        ]);

        assert_eq!(molecule.neighbors(0).unwrap().collect::<Vec<_>>(), [ 1, 2 ])
    }
}

#[cfg(test)]
mod has_node {
    use pretty_assertions::assert_eq;
    use super::*;

    #[test]
    fn unknown_id() {
        let molecule = DefaultMolecule::new(vec![ ]);

        assert_eq!(molecule.has_id(0), false)
    }

    #[test]
    fn known_id() {
        let molecule = DefaultMolecule::new(vec![
            Node { atom: Atom::default(), bonds: vec![ ] }
        ]);

        assert_eq!(molecule.has_id(0), true)
    }
}

#[cfg(test)]
mod degree {
    use pretty_assertions::assert_eq;
    use super::super::{ Node, Bond };
    use super::*;

    #[test]
    fn unknown_id() {
        let molecule = DefaultMolecule::new(vec![ ]);

        assert_eq!(molecule.degree(0), Err(GraphError::UnknownId(0)))
    }

    #[test]
    fn known_id() {
        let molecule = DefaultMolecule::new(vec![
            Node { atom: Atom::default(), bonds: vec![
                Bond::new(2, None, 1),
                Bond::new(2, None, 2)
            ]},
            Node { atom: Atom::default(), bonds: vec![
                Bond::new(2, None, 0)
            ]},
            Node { atom: Atom::default(), bonds: vec![
                Bond::new(2, None, 0)
            ]}
        ]);

        assert_eq!(molecule.degree(0), Ok(2))
    }
}

#[cfg(test)]
mod edges {
    use pretty_assertions::assert_eq;
    use super::super::{ Bond };
    use super::*;

    #[test]
    fn methane() {
        let molecule = DefaultMolecule::new(vec![
            Node { atom: Atom::default(), bonds: vec![ ] }
        ]);

        assert_eq!(molecule.edges().collect::<Vec<_>>(), [ ])
    }

    #[test]
    fn trimethylboron() {
        let molecule = DefaultMolecule::new(vec![
            Node { atom: Atom::default(), bonds: vec![
                Bond::new(2, None, 1),
                Bond::new(2, None, 2),
                Bond::new(2, None, 3)
            ]},
            Node { atom: Atom::default(), bonds: vec![
                Bond::new(2, None, 0)
            ]},
            Node { atom: Atom::default(), bonds: vec![
                Bond::new(2, None, 0)
            ]},
            Node { atom: Atom::default(), bonds: vec![
                Bond::new(2, None, 0)
            ]}
        ]);

        assert_eq!(molecule.edges().collect::<Vec<_>>(), [
            (0, 1), (0, 2), (0, 3)
        ])
    }
}

#[cfg(test)]
mod has_edge {
    use pretty_assertions::assert_eq;
    use super::super::{ Bond };
    use super::*;

    #[test]
    fn unknown_sid() {
        let molecule = DefaultMolecule::new(vec![
            Node { atom: Atom::default(), bonds: vec![ ] }
        ]);

        assert_eq!(molecule.has_edge(1, 0), Err(GraphError::UnknownId(1)))
    }

    #[test]
    fn unknown_tid() {
        let molecule = DefaultMolecule::new(vec![
            Node { atom: Atom::default(), bonds: vec![ ]}
        ]);

        assert_eq!(molecule.has_edge(0, 1), Err(GraphError::UnknownId(1)))
    }

    #[test]
    fn no_bond() {
        let molecule = DefaultMolecule::new(vec![
            Node { atom: Atom::default(), bonds: vec![ ] },
            Node { atom: Atom::default(), bonds: vec![ ] }
        ]);

        assert_eq!(molecule.has_edge(0, 1), Ok(false))
    }

    #[test]
    fn bond() {
        let molecule = DefaultMolecule::new(vec![
            Node { atom: Atom::default(), bonds: vec![ Bond::new(2, None, 1)] },
            Node { atom: Atom::default(), bonds: vec![ Bond::new(2, None, 0)] }
        ]);

        assert_eq!(molecule.has_edge(0, 1), Ok(true))
    }
}

#[cfg(test)]
mod atom {
    use pretty_assertions::assert_eq;
    use crate::molecule::{ Element };
    use super::*;

    #[test]
    fn unknown_id() {
        let molecule = DefaultMolecule::new(vec![ ]);

        assert_eq!(molecule.atom(0), Err(GraphError::UnknownId(0)))
    }

    #[test]
    fn methane() {
        let molecule = DefaultMolecule::new(vec![
            Node {
                atom: Atom {
                    isotope: None,
                    element: Some(Element::C),
                    hydrogens: 4,
                    electrons: 0,
                    parity: None,
                },
                bonds: vec![ ]
            }
        ]);

        assert_eq!(molecule.atom(0), Ok(&Atom {
            isotope: None,
            element: Some(Element::C),
            hydrogens: 4,
            electrons: 0,
            parity: None,
        }))
    }
}

#[cfg(test)]
mod charge {
    use pretty_assertions::assert_eq;
    use crate::molecule::{ Element, Bond };
    use super::*;

    #[test]
    fn unkown_id() {
        let molecule = DefaultMolecule::new(vec![ ]);
    
        assert_eq!(molecule.charge(0), Err(GraphError::UnknownId(0)))
    }

    #[test]
    fn ethyl_cation() {
        let molecule = DefaultMolecule::new(vec![
            Node {
                atom: Atom {
                    isotope: None,
                    element: Some(Element::C),
                    hydrogens: 3,
                    electrons: 0,
                    parity: None,
                },
                bonds: vec![ Bond::new(2, None, 1)]
            },
            Node {
                atom: Atom {
                    isotope: None,
                    element: Some(Element::C),
                    hydrogens: 2,
                    electrons: 0,
                    parity: None,
                },
                bonds: vec![ Bond::new(2, None, 0)]
            }
        ]);

        assert_eq!(molecule.charge(1), Ok(1f32))
    }

    #[test]
    fn ethyl_anion() {
        let molecule = DefaultMolecule::new(vec![
            Node {
                atom: Atom {
                    isotope: None,
                    element: Some(Element::C),
                    hydrogens: 3,
                    electrons: 0,
                    parity: None,
                },
                bonds: vec![ Bond::new(2, None, 1)]
            },
            Node {
                atom: Atom {
                    isotope: None,
                    element: Some(Element::C),
                    hydrogens: 2,
                    electrons: 2,
                    parity: None,
                },
                bonds: vec![ Bond::new(2, None, 0)]
            }
        ]);

        assert_eq!(molecule.charge(1), Ok(-1f32))
    }
}

#[cfg(test)]
mod bond_order {
    use pretty_assertions::assert_eq;
    use crate::molecule::{ Element, Bond };
    use super::*;

    #[test]
    fn unknown_sid() {
        let molecule = DefaultMolecule::new(vec![
            Node { atom: Atom::default(), bonds: vec![ ] }
        ]);

        assert_eq!(molecule.bond_order(1, 0), Err(GraphError::UnknownId(1)))
    }

    #[test]
    fn unknown_tid() {
        let molecule = DefaultMolecule::new(vec![
            Node { atom: Atom::default(), bonds: vec![ ] }
        ]);

        assert_eq!(molecule.bond_order(0, 1), Err(GraphError::UnknownId(1)))
    }

    #[test]
    fn no_bond() {
        let molecule = DefaultMolecule::new(vec![
            Node { atom: Atom::default(), bonds: vec![ ] },
            Node { atom: Atom::default(), bonds: vec![ ] }
        ]);

        assert_eq!(molecule.bond_order(0, 1), Ok(0f32))
    }

    #[test]
    fn bond() {
        let molecule = DefaultMolecule::new(vec![
            Node {
                atom: Atom {
                    isotope: None,
                    element: Some(Element::C),
                    hydrogens: 3,
                    electrons: 0,
                    parity: None,
                },
                bonds: vec![ Bond::new(2, None, 1)]
            },
            Node {
                atom: Atom {
                    isotope: None,
                    element: Some(Element::C),
                    hydrogens: 3,
                    electrons: 0,
                    parity: None,
                },
                bonds: vec![ Bond::new(2, None, 0)]
            }
        ]);

        assert_eq!(molecule.bond_order(0, 1), Ok(1f32))
    }

    #[test]
    fn half_bond() {
        let molecule = DefaultMolecule::new(vec![
            Node {
                atom: Atom {
                    isotope: None,
                    element: Some(Element::C),
                    hydrogens: 3,
                    electrons: 0,
                    parity: None,
                },
                bonds: vec![ Bond::new(1, None, 1)]
            },
            Node {
                atom: Atom {
                    isotope: None,
                    element: Some(Element::C),
                    hydrogens: 3,
                    electrons: 0,
                    parity: None,
                },
                bonds: vec![ Bond::new(1, None, 0)]
            }
        ]);

        assert_eq!(molecule.bond_order(0, 1), Ok(0.5f32))
    }
}