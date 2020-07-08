use std::collections::HashMap;

use gamma::graph::{ ArrayGraph, Graph, Error as GraphError };
use super::{ Molecule, Atom, BondOrder, Parity, Element, Error };

/// An implementation of the minimal molecule API concept.
/// 
/// ```rust
/// use chemcore::molecule::{
///     Molecule, DefaultMolecule, Element, Atom, BondOrder, Error
/// };
/// 
/// fn main() -> Result<(), Error> {
///     let molecule = DefaultMolecule::from_adjacency(vec![
///         (Atom::new(Element::C, 3, 0, None, None), vec![
///             (1, BondOrder::Single, None)
///         ]),
///         (Atom::new(Element::C, 3, 0, None, None), vec![
///             (0, BondOrder::Single, None)
///         ])
///     ])?;
/// 
///     assert_eq!(molecule.element(0), Ok(Element::C));
/// 
///     Ok(())
/// }
/// ```
#[derive(Debug, PartialEq)]
pub struct DefaultMolecule {
    nodes: Vec<Node>,
    bonds: HashMap<(usize, usize), (BondOrder, Option<Parity>)>,
    graph: ArrayGraph
}

impl DefaultMolecule {
    pub fn new() -> Self {
        DefaultMolecule {
            nodes: Vec::new(), bonds: HashMap::new(), graph: ArrayGraph::new()
        }
    }

    pub fn from_adjacency(
        entries: Vec<(Atom, Vec<(usize, BondOrder, Option<Parity>)>)>
    ) -> Result<Self, Error> {
        let mut nodes = Vec::new();
        let mut bonds = HashMap::new();
        let mut adjacency = Vec::new();

        for (atom, outs) in entries {
            let mut neighbors = Vec::new();
            let mut orders = Vec::new();

            for (tid, order, parity) in outs {
                neighbors.push(tid);
                orders.push(order);
                bonds.insert((nodes.len(), tid), (order, parity));
            }

            match Node::build(atom, orders) {
                Ok(node) => {
                    nodes.push(node);
                },
                Err(NodeError::Hypervalent) => {
                    return Err(Error::Hypervalent(nodes.len()));
                },
                Err(NodeError::ImpossibleIsotope) => {
                    return Err(Error::ImpossibleIsotope(nodes.len()));
                },
                Err(NodeError::ParityNotAllowed) => {
                    return Err(Error::AtomParityNotAllowed(nodes.len()));
                }
            }

            adjacency.push(neighbors);
        }

        let graph = ArrayGraph::from_adjacency(adjacency).unwrap();

        Ok(DefaultMolecule { nodes, graph, bonds })
    }

    fn node_at(&self, id: usize) -> Result<&Node, GraphError> {
        match self.nodes.get(id) {
            Some(node) => Ok(node),
            None => Err(GraphError::MissingNode(id))
        }
    }

    fn edge_at(
        &self, sid: usize, tid: usize
    ) -> Result<&(BondOrder, Option<Parity>), GraphError> {
        if !self.graph.has_node(sid) {
            Err(GraphError::MissingNode(sid))
        } else if !self.graph.has_node(tid) {
            Err(GraphError::MissingNode(tid))
        } else {
            match self.bonds.get(&(sid, tid)) {
                Some(edge) => Ok(edge),
                None => Err(GraphError::MissingEdge(sid, tid))
            }
        }
    }
}

#[derive(Debug, PartialEq)]
enum NodeError {
    Hypervalent,
    ImpossibleIsotope,
    ParityNotAllowed
}

#[derive(Debug, PartialEq)]
struct Node {
    element: Element,
    electrons: u8,
    charge: i8,
    hydrogens: u8,
    isotope: Option<u16>,
    parity: Option<Parity>
}

// TODO: parity should require a total of four substituents,
// including up to one hydrogen

impl Node {
    fn build(
        atom: Atom, orders: Vec<BondOrder>
    ) -> Result<Self, NodeError> {
        let element = atom.element;
        let hydrogens = atom.hydrogens;
        let isotope = atom.isotope;
        let parity = atom.parity;
        let mut electrons = element.valence_electrons() as i32;
        let valence = orders.iter().fold(0, |valence, order| {
            valence + order.multiplicity()
        });
        let ion = atom.ion;
        
        electrons -= hydrogens as i32;
        electrons -= ion as i32;
        electrons -= valence as i32;
        
        if electrons < 0 {
            return Err(NodeError::Hypervalent);
        }

        if let Some(isotope) = isotope {
            if isotope < element.atomic_number() {
                return Err(NodeError::ImpossibleIsotope);
            }
        }

        let substitution = hydrogens as usize + orders.len();

        if parity.is_some() && substitution != 4 {
            return Err(NodeError::ParityNotAllowed);
        }
        
        Ok(Node {
            element,
            charge: ion,
            electrons: electrons as u8,
            hydrogens,
            isotope,
            parity
        })
    }
}

impl Graph for DefaultMolecule {
    fn is_empty(&self) -> bool {
        self.graph.is_empty()
    }

    fn order(&self) -> usize {
        self.graph.order()
    }

    fn size(&self) -> usize { 
        self.graph.size()
    }

    fn nodes(&self) -> &[usize] {
        self.graph.nodes()
    }

    fn neighbors(&self, id: usize) -> Result<&[usize], GraphError> {
        self.graph.neighbors(id)
    }
    
    fn has_node(&self, id: usize) -> bool {
        self.graph.has_node(id)
    }

    fn degree(&self, id: usize) -> Result<usize, GraphError> {
        self.graph.degree(id)
    }

    fn edges(&self) -> &[(usize, usize)] {
        self.graph.edges()
    }

    fn has_edge(&self, sid: usize, tid: usize) -> Result<bool, GraphError> {
        self.graph.has_edge(sid, tid)
    }
}

impl Molecule for DefaultMolecule {
    fn element(&self, id: usize) -> Result<Element, GraphError> {
        Ok(self.node_at(id)?.element)
    }

    fn isotope(&self, id: usize) -> Result<Option<u16>, GraphError> {
        Ok(self.node_at(id)?.isotope)
    }
    
    fn electrons(&self, id: usize) -> Result<u8, GraphError> {
        Ok(self.node_at(id)?.electrons)
    }

    fn hydrogens(&self, id: usize) -> Result<u8, GraphError> {
        Ok(self.node_at(id)?.hydrogens)
    }

    fn charge(&self, id: usize) -> Result<i8, GraphError> {
        Ok(self.node_at(id)?.charge)
    }

    fn atom_parity(&self, id: usize) -> Result<Option<Parity>, GraphError> {
        Ok(self.node_at(id)?.parity)
    }

    fn bond_order(&self, sid: usize, tid: usize) -> Result<f32, GraphError> {
        Ok(self.edge_at(sid, tid)?.0.multiplicity() as f32)
    }
    
    fn bond_parity(
        &self, sid: usize, tid: usize
    ) -> Result<Option<Parity>, GraphError> {
        Ok(self.edge_at(sid, tid)?.1)
    }
}

#[cfg(test)]
mod molecule_tests {
    use super::*;
    use super::super::Element;

    #[test]
    fn from_adjacency_given_atom_initially_oversaturated() {
        let molecule = DefaultMolecule::from_adjacency(vec![
            (Atom::new(Element::C, 5, 0, None, None), vec![ ])
        ]);

        assert_eq!(molecule, Err(Error::Hypervalent(0)));
    }

    #[test]
    fn from_adjacency_given_atom_hypervalent() {
        let a0 = Atom::new(Element::C, 3, 0, None, None);
        let a1 = Atom::new(Element::C, 3, 0, None, None);

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![ (1, BondOrder::Double, None) ]),
            (a1, vec![ (0, BondOrder::Double, None) ])
        ]);

        assert_eq!(molecule, Err(Error::Hypervalent(0)));
    }

    #[test]
    fn from_adjacency_given_invalid_atom_parity() {
        let a0 = Atom::new(Element::C, 3, 0, None, Some(Parity::Negative));

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![ ])
        ]);

        assert_eq!(molecule, Err(Error::AtomParityNotAllowed(0)));
    }

    #[test]
    fn is_empty_given_methane() {
        let molecule = DefaultMolecule::from_adjacency(vec![
            (Atom::new(Element::C, 4, 0, None, None), vec![ ])
        ]).unwrap();

        assert_eq!(molecule.is_empty(), false);
    }

    #[test]
    fn order_given_methane() {
        let molecule = DefaultMolecule::from_adjacency(vec![
            (Atom::new(Element::C, 4, 0, None, None), vec![ ])
        ]).unwrap();

        assert_eq!(molecule.order(), 1);
    }

    #[test]
    fn size_given_ethane() {
        let a0 = Atom::new(Element::C, 3, 0, None, None);
        let a1 = Atom::new(Element::C, 3, 0, None, None);

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![ (1, BondOrder::Single, None) ]),
            (a1, vec![ (0, BondOrder::Single, None) ])
        ]).unwrap();

        assert_eq!(molecule.size(), 1);
    }

    #[test]
    fn nodes_given_ethane() {
        let a0 = Atom::new(Element::C, 3, 0, None, None);
        let a1 = Atom::new(Element::C, 3, 0, None, None);

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![ (1, BondOrder::Single, None) ]),
            (a1, vec![ (0, BondOrder::Single, None) ])
        ]).unwrap();

        assert_eq!(molecule.nodes(), &[ 0, 1 ]);
    }

    #[test]
    fn neighbors_given_propane() {
        let a0 = Atom::new(Element::C, 3, 0, None, None);
        let a1 = Atom::new(Element::C, 2, 0, None, None);
        let a2 = Atom::new(Element::C, 3, 0, None, None);

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![ (1, BondOrder::Single, None) ]),
            (a1, vec![
                (0, BondOrder::Single, None),
                (2, BondOrder::Single, None)
            ]),
            (a2, vec![ (1, BondOrder::Single, None) ])
        ]).unwrap();

        assert_eq!(molecule.neighbors(1).unwrap(), &[ 0, 2 ]);
    }

    #[test]
    fn has_node_given_methane() {
        let molecule = DefaultMolecule::from_adjacency(vec![
            (Atom::new(Element::C, 4, 0, None, None), vec![ ])
        ]).unwrap();

        assert_eq!(molecule.has_node(0), true);
    }

    #[test]
    fn degree_given_propane_secondary() {
        let a0 = Atom::new(Element::C, 3, 0, None, None);
        let a1 = Atom::new(Element::C, 2, 0, None, None);
        let a2 = Atom::new(Element::C, 3, 0, None, None);

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![ (1, BondOrder::Single, None) ]),
            (a1, vec![
                (0, BondOrder::Single, None),
                (2, BondOrder::Single, None)
            ]),
            (a2, vec![ (1, BondOrder::Single, None) ])
        ]).unwrap();

        assert_eq!(molecule.degree(1), Ok(2));
    }

    #[test]
    fn edges_given_propane() {
        let a0 = Atom::new(Element::C, 3, 0, None, None);
        let a1 = Atom::new(Element::C, 2, 0, None, None);
        let a2 = Atom::new(Element::C, 3, 0, None, None);

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![ (1, BondOrder::Single, None) ]),
            (a1, vec![
                (0, BondOrder::Single, None),
                (2, BondOrder::Single, None)
            ]),
            (a2, vec![ (1, BondOrder::Single, None) ])
        ]).unwrap();

        assert_eq!(molecule.edges(), &[
            (0, 1),
            (1, 2)
        ]);
    }

    #[test]
    fn has_edge_given_propane() {
        let a0 = Atom::new(Element::C, 3, 0, None, None);
        let a1 = Atom::new(Element::C, 2, 0, None, None);
        let a2 = Atom::new(Element::C, 3, 0, None, None);

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![ (1, BondOrder::Single, None) ]),
            (a1, vec![
                (0, BondOrder::Single, None),
                (2, BondOrder::Single, None)
            ]),
            (a2, vec![ (1, BondOrder::Single, None) ])
        ]).unwrap();

        assert_eq!(molecule.has_edge(0, 1), Ok(true));
    }

    #[test]
    fn element() {
        let molecule = DefaultMolecule::from_adjacency(vec![
            (Atom::new(Element::C, 0, 0, None, None), vec![ ])
        ]).unwrap();

        assert_eq!(molecule.element(0), Ok(Element::C));
    }

    #[test]
    fn hydrogens() {
        let molecule = DefaultMolecule::from_adjacency(vec![
            (Atom::new(Element::C, 4, 0, None, None), vec![  ])
        ]).unwrap();

        assert_eq!(molecule.hydrogens(0), Ok(4));
    }

    #[test]
    fn charge() {
        let molecule = DefaultMolecule::from_adjacency(vec![
            (Atom::new(Element::C, 3, 1, None, None), vec![  ])
        ]).unwrap();

        assert_eq!(molecule.charge(0), Ok(1));
    }

    #[test]
    fn electrons() {
        let molecule = DefaultMolecule::from_adjacency(vec![
            (Atom::new(Element::C, 0, 0, None, None), vec![ ])
        ]).unwrap();

        assert_eq!(molecule.electrons(0), Ok(4));
    }

    #[test]
    fn isotope() {
        let molecule = DefaultMolecule::from_adjacency(vec![
            (Atom::new(Element::C, 0, 0, Some(13), None), vec![ ])
        ]).unwrap();

        assert_eq!(molecule.isotope(0), Ok(Some(13)));
    }

    #[test]
    fn atom_parity() {
        let a0 = Atom::new(Element::C, 1, 0, None, Some(Parity::Negative));
        let a1 = Atom::new(Element::C, 0, 0, None, None);
        let a2 = Atom::new(Element::C, 0, 0, None, None);
        let a3 = Atom::new(Element::C, 0, 0, None, None);

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![
                (1, BondOrder::Single, None),
                (2, BondOrder::Single, None),
                (3, BondOrder::Single, None)
            ]),
            (a1, vec! [
                (0, BondOrder::Single, None)
            ]),
            (a2, vec![
                (0, BondOrder::Single, None)
            ]),
            (a3, vec![
                (0, BondOrder::Single, None)
            ])
        ]).unwrap();

        assert_eq!(molecule.atom_parity(0), Ok(Some(Parity::Negative)));
    }

    #[test]
    fn bond_order_given_outside_sid() {
        let a0 = Atom::new(Element::C, 3, 0, None, None);
        let a1 = Atom::new(Element::C, 3, 0, None, None);

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![ (1, BondOrder::Single, None) ]),
            (a1, vec![ (0, BondOrder::Single, None) ])
        ]).unwrap();
        let result = molecule.bond_order(2, 1);

        assert_eq!(result, Err(GraphError::MissingNode(2)));
    }

    #[test]
    fn bond_order_given_outside_tid() {
        let a0 = Atom::new(Element::C, 3, 0, None, None);
        let a1 = Atom::new(Element::C, 3, 0, None, None);

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![ (1, BondOrder::Single, None) ]),
            (a1, vec![ (0, BondOrder::Single, None) ])
        ]).unwrap();
        let result = molecule.bond_order(0, 2);

        assert_eq!(result, Err(GraphError::MissingNode(2)));
    }

    #[test]
    fn bond_order_given_no_edge() {
        let a0 = Atom::new(Element::C, 3, 0, None, None);
        let a1 = Atom::new(Element::C, 3, 0, None, None);

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![ ]),
            (a1, vec![ ])
        ]).unwrap();
        let result = molecule.bond_order(0, 1);

        assert_eq!(result, Err(GraphError::MissingEdge(0, 1)));
    }

    #[test]
    fn bond_order_given_inside() {
        let a0 = Atom::new(Element::C, 3, 0, None, None);
        let a1 = Atom::new(Element::C, 3, 0, None, None);

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![ (1, BondOrder::Single, None) ]),
            (a1, vec![ (0, BondOrder::Single, None) ])
        ]).unwrap();

        assert_eq!(molecule.bond_order(0, 1), Ok(1.0f32));
    }




    #[test]
    fn bond_parity_given_outside_sid() {
        let a0 = Atom::new(Element::C, 2, 0, None, None);
        let a1 = Atom::new(Element::C, 2, 0, None, None);

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![ (1, BondOrder::Double, Some(Parity::Positive)) ]),
            (a1, vec![ (0, BondOrder::Double, Some(Parity::Positive)) ])
        ]).unwrap();
        let result = molecule.bond_parity(2, 1);

        assert_eq!(result, Err(GraphError::MissingNode(2)));
    }

    #[test]
    fn bond_parity_given_outside_tid() {
        let a0 = Atom::new(Element::C, 2, 0, None, None);
        let a1 = Atom::new(Element::C, 2, 0, None, None);

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![ (1, BondOrder::Double, Some(Parity::Positive)) ]),
            (a1, vec![ (0, BondOrder::Double, Some(Parity::Positive)) ])
        ]).unwrap();
        let result = molecule.bond_parity(0, 2);

        assert_eq!(result, Err(GraphError::MissingNode(2)));
    }

    #[test]
    fn bond_parity_given_no_edge() {
        let a0 = Atom::new(Element::C, 2, 0, None, None);
        let a1 = Atom::new(Element::C, 2, 0, None, None);

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![ ]),
            (a1, vec![ ])
        ]).unwrap();
        let result = molecule.bond_parity(0, 1);

        assert_eq!(result, Err(GraphError::MissingEdge(0, 1)));
    }

    #[test]
    fn bond_parity_given_none() {
        let a0 = Atom::new(Element::C, 2, 0, None, None);
        let a1 = Atom::new(Element::C, 2, 0, None, None);

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![ (1, BondOrder::Double, None) ]),
            (a1, vec![ (0, BondOrder::Double, None) ])
        ]).unwrap();

        assert_eq!(molecule.bond_parity(0, 1), Ok(None));
    }

    #[test]
    fn bond_parity_given_some() {
        let a0 = Atom::new(Element::C, 2, 0, None, None);
        let a1 = Atom::new(Element::C, 2, 0, None, None);

        let molecule = DefaultMolecule::from_adjacency(vec![
            (a0, vec![ (1, BondOrder::Double, Some(Parity::Positive)) ]),
            (a1, vec![ (0, BondOrder::Double, Some(Parity::Positive)) ])
        ]).unwrap();

        assert_eq!(molecule.bond_parity(0, 1), Ok(Some(Parity::Positive)));
    }
}

#[cfg(test)]
mod node_tests {
    use super::*;

    #[test]
    fn build_given_oversaturated() {
        let atom = Atom::new(Element::C, 5, 0, None, None);

        assert_eq!(Node::build(atom, vec![ ]), Err(NodeError::Hypervalent));
    }

    #[test]
    fn build_given_overcharged() {
        let atom = Atom::new(Element::C, 0, 5, None, None);

        assert_eq!(Node::build(atom, vec![ ]), Err(NodeError::Hypervalent));
    }

    #[test]
    fn build_given_impossible_isotope() {
        let atom = Atom::new(Element::C, 0, 0, Some(5), None);

        assert_eq!(Node::build(atom, vec![ ]), Err(NodeError::ImpossibleIsotope));
    }

    #[test]
    fn build_given_parity_and_trivalent() {
        let atom = Atom::new(Element::C, 0, 0, None, Some(Parity::Positive));
        let node = Node::build(atom, vec![
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single
        ]);

        assert_eq!(node, Err(NodeError::ParityNotAllowed));
    }

    #[test]
    fn neutral() {
        let atom = Atom::new(Element::C, 0, 0, None, None);
        let node = Node::build(atom, vec![ ]).unwrap();

        assert_eq!(node, Node {
            element: Element::C,
            electrons: 4,
            hydrogens: 0,
            charge: 0,
            isotope: None,
            parity: None 
        });
    }

    #[test]
    fn neutral_with_single_bond() {
        let atom = Atom::new(Element::C, 0, 0, None, None);
        let node = Node::build(atom, vec![
            BondOrder::Single
        ]).unwrap();

        assert_eq!(node, Node {
            element: Element::C,
            electrons: 3,
            hydrogens: 0,
            charge: 0,
            isotope: None,
            parity: None 
        });
    }

    #[test]
    fn neutral_with_double_bond() {
        let atom = Atom::new(Element::C, 0, 0, None, None);
        let node = Node::build(atom, vec![
            BondOrder::Double
        ]).unwrap();

        assert_eq!(node, Node {
            element: Element::C,
            electrons: 2,
            hydrogens: 0,
            charge: 0,
            isotope: None,
            parity: None 
        });
    }

    #[test]
    fn neutral_with_triple_bond() {
        let atom = Atom::new(Element::C, 0, 0, None, None);
        let node = Node::build(atom, vec![
            BondOrder::Triple
        ]).unwrap();

        assert_eq!(node, Node {
            element: Element::C,
            electrons: 1,
            hydrogens: 0,
            charge: 0,
            isotope: None,
            parity: None 
        });
    }

    #[test]
    fn neutral_with_hydrogen() {
        let atom = Atom::new(Element::C, 1, 0, None, None);
        let node = Node::build(atom, vec![ ]).unwrap();

        assert_eq!(node, Node {
            element: Element::C,
            electrons: 3,
            hydrogens: 1,
            charge: 0,
            isotope: None,
            parity: None 
        });
    }

    #[test]
    fn positively_charged() {
        let atom = Atom::new(Element::C, 0, 1, None, None);
        let node = Node::build(atom, vec![ ]).unwrap();

        assert_eq!(node, Node {
            element: Element::C,
            electrons: 3,
            hydrogens: 0,
            charge: 1,
            isotope: None,
            parity: None 
        });
    }

    #[test]
    fn positively_charged_with_hydrogen() {
        let atom = Atom::new(Element::C, 3, 1, None, None);
        let node = Node::build(atom, vec![ ]).unwrap();

        assert_eq!(node, Node {
            element: Element::C,
            electrons: 0,
            hydrogens: 3,
            charge: 1,
            isotope: None,
            parity: None 
        });
    }

    #[test]
    fn negatively_charged() {
        let atom = Atom::new(Element::C, 0, -1, None, None);
        let node = Node::build(atom, vec![ ]).unwrap();

        assert_eq!(node, Node {
            element: Element::C,
            electrons: 5,
            hydrogens: 0,
            charge: -1,
            isotope: None,
            parity: None 
        });
    }

    #[test]
    fn negatively_charged_with_hydrogen() {
        let atom = Atom::new(Element::C, 3, -1, None, None);
        let node = Node::build(atom, vec![ ]).unwrap();

        assert_eq!(node, Node {
            element: Element::C,
            electrons: 2,
            hydrogens: 3,
            charge: -1,
            isotope: None,
            parity: None 
        });
    }

    #[test]
    fn isotope() {
        let atom = Atom::new(Element::C, 0, 0, Some(13), None);
        let node = Node::build(atom, vec![ ]).unwrap();

        assert_eq!(node, Node {
            element: Element::C,
            electrons: 4,
            hydrogens: 0,
            charge: 0,
            isotope: Some(13),
            parity: None 
        });
    }

    #[test]
    fn parity() {
        let atom = Atom::new(Element::C, 0, 0, None, Some(Parity::Negative));
        let node = Node::build(atom, vec![
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single
        ]).unwrap();

        assert_eq!(node, Node {
            element: Element::C,
            electrons: 0,
            hydrogens: 0,
            charge: 0,
            isotope: None,
            parity: Some(Parity::Negative) 
        });
    }
}