use gamma::graph::Graph;
use gamma::graph::Error as GraphError;

use super::spec;
use super::molecule::Molecule;
use super::element::Element;
use super::error::Error;
use super::atom::Atom;
use super::parity::Parity;

/// An implementation of the minimal molecule API concept.
pub struct DefaultMolecule {
    atoms: Vec<Atom>,
    indices: Vec<usize>,
    edges: Vec<(usize, usize)>
}

impl DefaultMolecule {
    pub fn build(spec: spec::Molecule) -> Result<Self, Error> {
        let mut atoms = Vec::new();

        for atom_spec in spec.atoms {
            match Atom::build(atom_spec) {
                Ok(node) => atoms.push(node),
                Err(error) => return Err(error)
            }
        }

        let mut edges = Vec::new();

        for spec::Bond { sid, tid, order, parity } in spec.bonds {
            match atoms.get_mut(sid) {
                Some(source) => {
                    if let Err(error) = source.add_bond(tid, order, parity) {
                        return Err(error);
                    }
                },
                None => return Err(Error::MisplacedBond)
            }

            match atoms.get_mut(tid) {
                Some(target) => {
                    if let Err(error) = target.add_bond(sid, order, parity) {
                        return Err(error);
                    }
                },
                None => return Err(Error::MisplacedBond)
            }

            edges.push((sid, tid));
        }

        let indices = (0..atoms.len()).collect::<Vec<usize>>();

        Ok(DefaultMolecule { atoms, indices, edges })
    }
}

pub struct EdgeIterator<'a> {
    cursor: usize,
    edges: &'a Vec<(usize, usize)>
}

impl<'a> EdgeIterator<'a> {
    pub fn new(edges: &'a Vec<(usize, usize)>) -> Self {
        EdgeIterator { cursor: 0, edges }
    }
}

impl<'a> Iterator for EdgeIterator<'a> {
    type Item = (&'a usize, &'a usize);

    fn next(&mut self) -> Option<(&'a usize, &'a usize)> {
        if let Some((sid, tid)) = self.edges.get(self.cursor) {
            self.cursor += 1;

            Some((&sid, &tid))
        } else {
            None
        }
    }
}

impl<'a> Graph<'a, usize> for DefaultMolecule {
    type NodeIterator = std::slice::Iter<'a, usize>;
    type NeighborIterator = std::slice::Iter<'a, usize>;
    type EdgeIterator = EdgeIterator<'a>;

    fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    fn order(&self) -> usize {
        self.atoms.len()
    }

    fn size(&self) -> usize {
        self.edges.len()
    }

    fn nodes(&'a self) -> Self::NodeIterator {
        self.indices.iter()
    }

    fn edges(&'a self) -> Self::EdgeIterator {
        EdgeIterator::new(&self.edges)
    }

    fn has_node(&self, id: &usize) -> bool {
        match self.atoms.get(*id) {
            Some(_) => true,
            None => false
        }
    }

    fn neighbors(
        &'a self, id: &usize
    ) -> Result<Self::NeighborIterator, GraphError> {
        match self.atoms.get(*id) {
            Some(atom) => Ok(atom.neighbors()),
            None => Err(GraphError::UnknownNode)
        }
    }

    fn degree(&self, id: &usize) -> Result<usize, GraphError> {
        match self.atoms.get(*id) {
            Some(atom) => Ok(atom.degree()),
            None => Err(GraphError::UnknownNode)
        }
    }

    fn has_edge(&self, sid: &usize, tid: &usize) -> Result<bool, GraphError> {
        let target = self.atoms.get(*tid);
        
        if target.is_none() {
            return Err(GraphError::UnknownNode);
        }

        match self.atoms.get(*sid) {
            Some(atom) => Ok(atom.has_edge(*tid)),
            None => Err(GraphError::UnknownNode)
        }
    }
}

impl<'a> Molecule<'a, usize> for DefaultMolecule {
    fn element(&self, id: &usize) -> Result<Element, GraphError> {
        match self.atoms.get(*id) {
            Some(node) => Ok(node.element()),
            None => Err(GraphError::UnknownNode)
        }
    }

    fn isotope(&self, id: &usize) -> Result<Option<u16>, GraphError> {
        match self.atoms.get(*id) {
            Some(atom) => Ok(atom.isotope()),
            None => Err(GraphError::UnknownNode)
        }
    }

    fn electrons(&self, id: &usize) -> Result<u8, GraphError> {
        match self.atoms.get(*id) {
            Some(atom) => Ok(atom.electrons()),
            None => Err(GraphError::UnknownNode)
        }
    }

    fn hydrogens(&self, id: &usize) -> Result<u8, GraphError> {
        match self.atoms.get(*id) {
            Some(atom) => Ok(atom.hydrogens()),
            None => Err(GraphError::UnknownNode)
        }
    }

    fn charge(&self, id: &usize) -> Result<i8, GraphError> {
        match self.atoms.get(*id) {
            Some(atom) => Ok(atom.charge()),
            None => Err(GraphError::UnknownNode)
        }
    }

    fn atom_parity(&self, id: &usize) -> Result<Option<Parity>, GraphError> {
        match self.atoms.get(*id) {
            Some(atom) => Ok(atom.atom_parity()),
            None => Err(GraphError::UnknownNode)
        }
    }

    fn bond_order(&self, sid: &usize, tid: &usize) -> Result<u8, GraphError> {
        let target = self.atoms.get(*tid);
        
        if target.is_none() {
            return Err(GraphError::UnknownNode);
        }

        match self.atoms.get(*sid) {
            Some(atom) => Ok(atom.bond_order(tid)),
            None => Err(GraphError::UnknownNode)
        }
    }

    fn bond_parity(
        &self, sid: &usize, tid: &usize
    ) -> Result<Option<Parity>, GraphError> {
        let target = self.atoms.get(*tid);

        if target.is_none() {
            return Err(GraphError::UnknownNode);
        }

        match self.atoms.get(*sid) {
            Some(atom) => Ok(atom.bond_parity(tid)),
            None => Err(GraphError::UnknownNode)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::spec;
    use crate::molecule::BondOrder;

    #[test]
    fn build_returns_error_given_too_many_hydrogens() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 5, ..Default::default()
                }
            ],
            bonds: vec![ ]
        });

        assert_eq!(molecule.err(), Some(Error::HypervalentAtom));
    }

    #[test]
    fn build_returns_error_given_invalid_bond_sid() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 2, tid: 1, ..Default::default() }
            ]
        });

        assert_eq!(molecule.err(), Some(Error::MisplacedBond));
    }

    #[test]
    fn build_returns_error_given_invalid_bond_tid() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 0, tid: 2, ..Default::default() }
            ]
        });

        assert_eq!(molecule.err(), Some(Error::MisplacedBond));
    }

    #[test]
    fn element_given_invalid_id() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 4, ..Default::default()
                }
            ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.element(&1), Err(GraphError::UnknownNode));
    }

    #[test]
    fn element_given_valid_id() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, ..Default::default()
                }
            ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.element(&0), Ok(Element::C));
    }

    #[test]
    fn isotope_given_invalid_id() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 4, ..Default::default()
                }
            ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.isotope(&1), Err(GraphError::UnknownNode));
    }

    #[test]
    fn isotope_given_valid_id() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, isotope: Some(13), ..Default::default()
                }
            ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.isotope(&0), Ok(Some(13)));
    }

    #[test]
    fn electrons_given_invalid_id() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 0, ..Default::default()
                }
            ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.electrons(&1), Err(GraphError::UnknownNode));
    }

    #[test]
    fn electrons_given_helium() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::He, hydrogens: 0, ..Default::default()
                }
            ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.electrons(&0), Ok(2));
    }

    #[test]
    fn hydrogens_given_invalid_id() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 4, ..Default::default()
                }
            ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.hydrogens(&1), Err(GraphError::UnknownNode));
    }

    #[test]
    fn hydrogens_given_methane() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 4, ..Default::default()
                }
            ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.hydrogens(&0), Ok(4));
    }

    #[test]
    fn charge_given_invalid_id() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 4, ..Default::default()
                }
            ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.charge(&1), Err(GraphError::UnknownNode));
    }

    #[test]
    fn charge_given_methyl_cation() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ion: 1, isotope: None,
                    parity: None
                }
            ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.charge(&0), Ok(1));
    }

    #[test]
    fn charge_given_foo() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 2, ..Default::default()
                },
                spec::Atom {
                    element: Element::O, hydrogens: 2, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond {
                    sid: 0, tid: 1, order: BondOrder::Double, parity: None
                }
            ]
        }).unwrap();

        assert_eq!(molecule.charge(&0), Ok(0));
    }

    #[test]
    fn atom_parity_given_invalid_id() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 4, ..Default::default()
                }
            ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.atom_parity(&1), Err(GraphError::UnknownNode));
    }

    #[test]
    fn atom_parity_given_methane() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ion: 1, isotope: None,
                    parity: None
                }
            ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.atom_parity(&0), Ok(None));
    }

    #[test]
    fn atom_parity_given_bromochlorofluoromethane() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 1, ion: 0, isotope: None,
                    parity: Some(Parity::Positive)
                },
                spec::Atom {
                    element: Element::Br, ..Default::default()
                },
                spec::Atom {
                    element: Element::Cl, ..Default::default()
                },
                spec::Atom {
                    element: Element::F, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 0, tid: 1, ..Default::default() },
                spec::Bond { sid: 0, tid: 2, ..Default::default() },
                spec::Bond { sid: 0, tid: 3, ..Default::default() }
            ]
        }).unwrap();

        assert_eq!(molecule.atom_parity(&0), Ok(Some(Parity::Positive)));
    }

    #[test]
    fn bond_order_given_invalid_sid() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 0, tid: 1, ..Default::default() }
            ]
        }).unwrap();

        assert_eq!(molecule.bond_order(&2, &0), Err(GraphError::UnknownNode))
    }

    #[test]
    fn bond_order_given_invalid_tid() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 0, tid: 1, ..Default::default() }
            ]
        }).unwrap();

        assert_eq!(molecule.bond_order(&0, &2), Err(GraphError::UnknownNode));
    }

    #[test]
    fn bond_order_given_disconnected_methanes() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 4, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 4, ..Default::default()
                }
            ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.bond_order(&0, &1), Ok(0));
    }

    #[test]
    fn bond_order_given_zero() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 2, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 2, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond {
                    sid: 0, tid: 1, order: BondOrder::Zero, parity: None
                }
            ]
        }).unwrap();

        assert_eq!(molecule.bond_order(&0, &1), Ok(0));
    }

    #[test]
    fn bond_order_given_ethane() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 0, tid: 1, ..Default::default() }
            ]
        }).unwrap();

        assert_eq!(molecule.bond_order(&0, &1), Ok(1));
    }

    #[test]
    fn bond_order_given_ethane_reversed() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 0, tid: 1, ..Default::default() }
            ]
        }).unwrap();

        assert_eq!(molecule.bond_order(&1, &0), Ok(1));
    }

    #[test]
    fn bond_parity_given_invalid_sid() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 0, tid: 1, ..Default::default() }
            ]
        }).unwrap();

        assert_eq!(molecule.bond_parity(&2, &0), Err(GraphError::UnknownNode))
    }

    #[test]
    fn bond_parity_given_invalid_tid() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 0, tid: 1, ..Default::default() }
            ]
        }).unwrap();

        assert_eq!(molecule.bond_parity(&0, &2), Err(GraphError::UnknownNode));
    }

    #[test]
    fn bond_parity_given_ethane() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond {
                    sid: 0, tid: 1, order: BondOrder::Single,
                    parity: Some(Parity::Negative)
                }
            ]
        }).unwrap();

        assert_eq!(molecule.bond_parity(&0, &1), Ok(Some(Parity::Negative)));
    }

    #[test]
    fn is_empty_given_empty() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![ ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.is_empty(), true);
    }

    #[test]
    fn is_empty_given_methane() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 4, ..Default::default()
                }
            ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.is_empty(), false);
    }

    #[test]
    fn size_given_empty() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![ ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.size(), 0);
    }

    #[test]
    fn size_given_ethane() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 0, tid: 1, ..Default::default() }
            ]
        }).unwrap();

        assert_eq!(molecule.size(), 1);
    }

    #[test]
    fn order_given_empty() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![ ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.order(), 0);
    }

    #[test]
    fn order_given_ethane() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 0, tid: 1, ..Default::default() }
            ]
        }).unwrap();

        assert_eq!(molecule.order(), 2);
    }

    #[test]
    fn nodes_given_empty() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![ ],
            bonds: vec![ ]
        }).unwrap();
        let nodes = molecule.nodes().collect::<Vec<_>>();

        assert_eq!(nodes.is_empty(), true);
    }

    #[test]
    fn nodes_given_methane() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 4, ..Default::default()
                }
            ],
            bonds: vec![ ]
        }).unwrap();
        let nodes = molecule.nodes().collect::<Vec<_>>();

        assert_eq!(nodes, vec![ &0 ]);
    }

    #[test]
    fn edges_given_empty() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![ ],
            bonds: vec![ ]
        }).unwrap();
        let edges = molecule.edges().collect::<Vec<_>>();

        assert_eq!(edges.is_empty(), true);
    }

    #[test]
    fn edges_given_ethane() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 0, tid: 1, ..Default::default() }
            ]
        }).unwrap();
        let edges = molecule.edges().collect::<Vec<_>>();

        assert_eq!(edges, vec![ ( &0, &1 ) ]);
    }

    #[test]
    fn has_node_given_empty() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![ ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.has_node(&0), false);
    }

    #[test]
    fn has_node_given_methane() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 4, ..Default::default()
                }
            ],
            bonds: vec![ ]
        }).unwrap();

        assert_eq!(molecule.has_node(&0), true);
    }

    #[test]
    fn degree_given_empty() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![ ],
            bonds: vec![ ]
        }).unwrap();
        let degree = molecule.degree(&0);

        assert_eq!(degree.err(), Some(GraphError::UnknownNode));
    }

    #[test]
    fn degree_given_ethane() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 0, tid: 1, ..Default::default() }
            ]
        }).unwrap();

        assert_eq!(molecule.degree(&0), Ok(1));
    }

    #[test]
    fn neighbors_given_empty() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![ ],
            bonds: vec![ ]
        }).unwrap();
        let neighbors = molecule.neighbors(&0);

        assert_eq!(neighbors.err(), Some(GraphError::UnknownNode))
    }

    #[test]
    fn neighbors_given_ethane() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 0, tid: 1, ..Default::default() }
            ]
        }).unwrap();
        let neighbors = molecule.neighbors(&0).unwrap().collect::<Vec<_>>();;

        assert_eq!(neighbors, vec![ &1 ]);
    }

    #[test]
    fn has_edge_given_invalid_sid() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 0, tid: 1, ..Default::default() }
            ]
        }).unwrap();

        assert_eq!(molecule.has_edge(&2, &1), Err(GraphError::UnknownNode));
    }

    #[test]
    fn has_edge_given_invalid_tid() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 0, tid: 1, ..Default::default() }
            ]
        }).unwrap();

        assert_eq!(molecule.has_edge(&0, &2), Err(GraphError::UnknownNode));
    }

    #[test]
    fn has_edge_given_no_edge() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom { 
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 2, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 0, tid: 1, ..Default::default() },
                spec::Bond { sid: 1, tid: 2, ..Default::default() }
            ]
        }).unwrap();

        assert_eq!(molecule.has_edge(&0, &2), Ok(false));
    }

    #[test]
    fn has_edge_given_edge() {
        let molecule = DefaultMolecule::build(spec::Molecule {
            atoms: vec![
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                },
                spec::Atom {
                    element: Element::C, hydrogens: 3, ..Default::default()
                }
            ],
            bonds: vec![
                spec::Bond { sid: 0, tid: 1, ..Default::default() }
            ]
        }).unwrap();

        assert_eq!(molecule.has_edge(&0, &1), Ok(true));
    }
}