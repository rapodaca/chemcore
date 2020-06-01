use purr::read::read;
use purr::mol::{ Mol, Atom, Bond, Style, Parity as PurrParity };
use purr::valence::implicit_hydrogens;
use gamma::graph::Graph;
use gamma::graph::Error as GraphError;
use crate::molecule::{ Element, Parity };
use crate::molecule::Molecule as IMolecule;

#[derive(Debug, PartialEq, Eq)]
pub enum Error {
    
}

pub struct Molecule {
    mol: Mol,
    neighbors: Vec<Vec<usize>>,
    indices: Vec<usize>,
    edges: Vec<(usize, usize)>
}

impl Molecule {
    pub fn build(smiles: &str) -> Result<Self, Error> {
        let mol = read(smiles).unwrap();
        let mut edges = vec![];
        let mut neighbors = vec![ ];

        for atom in mol.atoms.iter() {
            if atom.aromatic {
                panic!("Aromatic SMILES not yet supported");
            }
        }

        for (sid, bonds) in mol.bonds.iter().enumerate() {
            for bond in bonds.iter() {
                if sid < bond.tid {
                    edges.push((sid, bond.tid));
                }

                if let Some(style) = bond.style {
                    if style == Style::Aromatic {
                        panic!("Aromatic SMILES not yet supported");
                    }
                }
            }

            neighbors.push(
                bonds.iter().map(|bond| bond.tid).collect()
            );
        }

        let indices = (0..mol.atoms.len()).collect::<Vec<usize>>();

        Ok(Molecule {
            mol, edges, indices, neighbors
        })
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

impl<'a> Graph<'a, usize> for Molecule {
    type NodeIterator = std::slice::Iter<'a, usize>;
    type NeighborIterator = std::slice::Iter<'a, usize>;
    type EdgeIterator = EdgeIterator<'a>;

    fn is_empty(&self) -> bool {
        self.mol.atoms.is_empty()
    }

    fn order(&self) -> usize {
        self.mol.atoms.len()
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
        match self.mol.atoms.get(*id) {
            Some(_) => true,
            None => false
        }
    }

    fn neighbors(
        &'a self, id: &usize
    ) -> Result<Self::NeighborIterator, GraphError> {
        Ok(get_neighbors(id, &self.neighbors)?.iter())
    }

    fn degree(&self, id: &usize) -> Result<usize, GraphError> {
        Ok(get_neighbors(id, &self.neighbors)?.len())
    }

    fn has_edge(&self, sid: &usize, tid: &usize) -> Result<bool, GraphError> {
        match self.neighbors.get(*tid) {
            Some(_) => Ok(get_neighbors(sid, &self.neighbors)?.contains(tid)),
            None => Err(GraphError::UnknownNode)
        }
    }
}

impl<'a> IMolecule<'a, usize> for Molecule {
    fn element(&self, id: &usize) -> Result<Element, GraphError> {
        let atom = atom_at(id, &self.mol)?;

        Ok(Element::from(atom.element))
    }

    fn isotope(&self, id: &usize) -> Result<Option<u16>, GraphError> {
        let atom = atom_at(id, &self.mol)?;

        Ok(atom.isotope)
    }

    fn electrons(&self, id: &usize) -> Result<u8, GraphError> {
        let atom = atom_at(id, &self.mol)?;
        let element = Element::from(atom.element);
        let mut result = element.valence_electrons() as i8;
        
        if let Some(charge) = atom.charge {
            result = result - charge;
        }

        result = result - valence(self.mol.bonds.get(*id).unwrap()) as i8;

        if let Some(hcount) = atom.hcount {
            result = result - hcount as i8;
        } else {
            let implicit = implicit_hydrogens(id, &self.mol);
            result = result - implicit.unwrap().unwrap() as i8;
        }

        Ok(result as u8)
    }

    fn hydrogens(&self, id: &usize) -> Result<u8, GraphError> {
        if let Ok(result) = implicit_hydrogens(id, &self.mol) {
            if let Some(value) = result {
                Ok(value)
            } else {
                let atom = atom_at(id, &self.mol)?;

                match atom.hcount {
                    Some(hcount) => Ok(hcount),
                    None => Ok(0)
                }
            }
        } else {
            Err(GraphError::UnknownNode)
        }
    }

    fn charge(&self, id: &usize) -> Result<i8, GraphError> {
        let atom = atom_at(id, &self.mol)?;

        match atom.charge {
            Some(charge) => Ok(charge),
            None => Ok(0)
        }
    }

    fn atom_parity(&self, id: &usize) -> Result<Option<Parity>, GraphError> {
        let atom = atom_at(id, &self.mol)?;

        if let Some(parity) = atom.parity {
            let mut result = match parity {
                PurrParity::Counterclockwise => Parity::Negative,
                PurrParity::Clockwise => Parity::Positive
            };

            if *id != 0 && atom.hcount.is_some() {
                // result = Parity::negate(result);
                result = result.negate();
            }

            Ok(Some(result))
        } else {
            Ok(None)
        }
    }

    fn bond_order(&self, sid: &usize, tid: &usize) -> Result<u8, GraphError> {
        match bond_at(sid, tid, &self.mol)? {
            Some(bond) => Ok(multiplicity(bond)),
            None => Ok(0)
        }
    }

    fn bond_parity(&self, sid: &usize, tid: &usize) -> Result<Option<Parity>, GraphError> {
        match bond_at(sid, tid, &self.mol)?.and_then(|bond| bond.style) {
            Some(Style::Double) => (),
            _ => return Ok(None)
        }

        let left = score_half_bond(self.mol.bonds.get(*sid).unwrap());
        let right = score_half_bond(self.mol.bonds.get(*tid).unwrap());

        match left {
            Some(p1) => {
                match right {
                    Some(p2) => {
                        Ok(Some(p1.multiply(&p2)))
                    },
                    None => Ok(None)
                }
            },
            None => Ok(None)
        }
    }
}

fn get_neighbors<'a>(id: &usize, neighbors: &'a Vec<Vec<usize>>) -> Result<&'a Vec<usize>, GraphError> {
    match neighbors.get(*id) {
        Some(neighbors) => Ok(neighbors),
        None => Err(GraphError::UnknownNode)
    }
}

fn atom_at<'a>(id: &usize, mol: &'a Mol) -> Result<&'a Atom, GraphError> {
    match mol.atoms.get(*id) {
        Some(atom) => Ok(atom),
        None => Err(GraphError::UnknownNode)
    }
}

fn bond_at<'a>(sid: &usize, tid: &usize, mol: &'a Mol) -> Result<Option<&'a Bond>, GraphError> {
    if let Some(bonds) = mol.bonds.get(*sid) {
        if mol.bonds.get(*tid).is_none() {
            Err(GraphError::UnknownNode)
        } else {
            for bond in bonds.iter() {
                if &bond.tid == tid {
                    return Ok(Some(bond));
                }
            }
    
            Ok(None)
        }
    } else {
        Err(GraphError::UnknownNode)
    }
}

fn valence(bonds: &Vec<Bond>) -> u8 {
    bonds.iter().fold(0, |total, bond| total + multiplicity(bond))
}

fn score_half_bond(bonds: &Vec<Bond>) -> Option<Parity> {
    match low_bond(bonds) {
        Some(bond) => match bond.style {
            Some(Style::Up) => Some(Parity::Positive),
            Some(Style::Down) => Some(Parity::Negative),
            _ => None
        },
        None => None
    }
}

fn low_bond(bonds: &Vec<Bond>) -> Option<&Bond> {
    let mut up_and_down = bonds.iter().filter(|bond| {
        match bond.style {
            Some(Style::Up) => true,
            Some(Style::Down) => true,
            _ => false
        }
    }).collect::<Vec<&Bond>>();

    up_and_down.sort_by(|a, b| a.tid.cmp(&b.tid));

    up_and_down.get(0).cloned()
}

fn multiplicity(bond: &Bond) -> u8 {
    match bond.style {
        Some(Style::Single) => 1,
        Some(Style::Double) => 2,
        Some(Style::Triple) => 3,
        Some(Style::Quadruple) => 4,
        Some(Style::Aromatic) => 1,
        Some(Style::Up) => 1,
        Some(Style::Down) => 1,
        None => 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic(expected = "Aromatic SMILES not yet supported")]
    fn aromatic_atom() {
        Molecule::build(&"cc").unwrap();
    }

    #[test]
    #[should_panic(expected = "Aromatic SMILES not yet supported")]
    fn aromatic_bond() {
        Molecule::build(&"C:C").unwrap();
    }

    #[test]
    fn element_given_unknown_id() {
        let molecule = Molecule::build(&"C").unwrap();
        
        assert_eq!(molecule.element(&1), Err(GraphError::UnknownNode));
    }

    #[test]
    fn element_given_known_id() {
        let molecule = Molecule::build(&"N").unwrap();

        assert_eq!(molecule.element(&0), Ok(Element::N));
    }

    #[test]
    fn isotope_given_unknown_id() {
        let molecule = Molecule::build(&"C").unwrap();

        assert_eq!(molecule.isotope(&1), Err(GraphError::UnknownNode));
    }

    #[test]
    fn isotope_given_none() {
        let molecule = Molecule::build(&"C").unwrap();

        assert_eq!(molecule.isotope(&0), Ok(None));
    }

    #[test]
    fn isotope_given_some() {
        let molecule = Molecule::build(&"[13C]").unwrap();
        
        assert_eq!(molecule.isotope(&0), Ok(Some(13)));
    }

    #[test]
    fn electrons_given_unknown_id() {
        let molecule = Molecule::build(&"C").unwrap();
        
        assert_eq!(molecule.electrons(&1), Err(GraphError::UnknownNode));
    }

    #[test]
    fn electrons_given_methane() {
        let molecule = Molecule::build(&"C").unwrap();

        assert_eq!(molecule.electrons(&0), Ok(0));
    }

    #[test]
    fn electrons_given_methyl_cation() {
        let molecule = Molecule::build(&"[CH3+]").unwrap();

        assert_eq!(molecule.electrons(&0), Ok(0));
    }

    #[test]
    fn electrons_given_methyl_anion() {
        let molecule = Molecule::build(&"[CH3-]").unwrap();

        assert_eq!(molecule.electrons(&0), Ok(2));
    }

    #[test]
    fn electrons_given_ammonia() {
        let molecule = Molecule::build(&"N").unwrap();

        assert_eq!(molecule.electrons(&0), Ok(2));
    }

    #[test]
    fn electrons_given_dimethylamine() {
        let molecule = Molecule::build(&"N(C)C").unwrap();

        assert_eq!(molecule.electrons(&0), Ok(2));
    }

    #[test]
    fn electrons_given_water() {
        let molecule = Molecule::build(&"O").unwrap();

        assert_eq!(molecule.electrons(&0), Ok(4));
    }

    #[test]
    fn hydrogens_given_unknown_id() {
        let molecule = Molecule::build(&"C").unwrap();
        
        assert_eq!(molecule.hydrogens(&1), Err(GraphError::UnknownNode));
    }

    #[test]
    fn hydrogens_given_methane() {
        let molecule = Molecule::build(&"C").unwrap();

        assert_eq!(molecule.hydrogens(&0), Ok(4));
    }

    #[test]
    fn hydrogens_given_methyl_cation() {
        let molecule = Molecule::build(&"[CH3+]").unwrap();

        assert_eq!(molecule.hydrogens(&0), Ok(3));
    }

    #[test]
    fn hydrogens_given_dimethylamine() {
        let molecule = Molecule::build(&"N(C)C").unwrap();

        assert_eq!(molecule.hydrogens(&0), Ok(1));
    }

    #[test]
    fn charge_given_unknown_id() {
        let molecule = Molecule::build(&"C").unwrap();

        assert_eq!(molecule.charge(&1), Err(GraphError::UnknownNode));
    }
    
    #[test]
    fn charge_given_methane() {
        let molecule = Molecule::build(&"C").unwrap();

        assert_eq!(molecule.charge(&0), Ok(0));
    }

    #[test]
    fn charge_given_methyl_cation() {
        let molecule = Molecule::build(&"[CH3+]").unwrap();

        assert_eq!(molecule.charge(&0), Ok(1));
    }

    #[test]
    fn atom_parity_given_unknown_id() {
        let molecule = Molecule::build(&"[CH4]").unwrap();

        assert_eq!(molecule.atom_parity(&1), Err(GraphError::UnknownNode));
    }

    #[test]
    fn atom_parity_given_none() {
        let molecule = Molecule::build(&"[CH4]").unwrap();

        assert_eq!(molecule.atom_parity(&0), Ok(None));
    }

    #[test]
    fn atom_parity_given_counterclockwise_leading_with_hydrogen() {
        let molecule = Molecule::build(&"[C@H](N)(O)C").unwrap();

        assert_eq!(molecule.atom_parity(&0), Ok(Some(Parity::Negative)));
    }

    #[test]
    fn atom_parity_given_clockwise_leading_with_hydrogen() {
        let molecule = Molecule::build(&"[C@@H](N)(O)C").unwrap();

        assert_eq!(molecule.atom_parity(&0), Ok(Some(Parity::Positive)));
    }

    #[test]
    fn atom_parity_given_clockwise_trailing_with_hydrogen() {
        let molecule = Molecule::build(&"N[C@@H](O)C").unwrap();

        assert_eq!(molecule.atom_parity(&1), Ok(Some(Parity::Negative)));
    }

    #[test]
    fn atom_parity_given_counterclockwise_trailing_with_hydrogen() {
        let molecule = Molecule::build(&"N[C@H](O)C").unwrap();

        assert_eq!(molecule.atom_parity(&1), Ok(Some(Parity::Positive)));
    }

    #[test]
    fn atom_parity_given_counterclockwise_trailing_without_hydrogen() {
        let molecule = Molecule::build(&"N[C@]([H])(O)C").unwrap();

        assert_eq!(molecule.atom_parity(&1), Ok(Some(Parity::Negative)));
    }

    #[test]
    fn atom_parity_given_clockwise_trailing_without_hydrogen() {
        let molecule = Molecule::build(&"N[C@@]([H])(O)C").unwrap();

        assert_eq!(molecule.atom_parity(&1), Ok(Some(Parity::Positive)));
    }

    #[test]
    fn atom_parity_given_counterclockwise_leading_with_hydrogen_and_cut() {
        let molecule = Molecule::build(&"[C@H](C)1OC1").unwrap();

        assert_eq!(molecule.atom_parity(&0), Ok(Some(Parity::Negative)));
    }

    #[test]
    fn atom_parity_given_clockwise_leading_with_hydrogen_and_cut() {
        let molecule = Molecule::build(&"[C@@H](C)1OC1").unwrap();

        assert_eq!(molecule.atom_parity(&0), Ok(Some(Parity::Positive)));
    }

    #[test]
    fn atom_parity_given_counterclockwise_trailing_without_hydrogen_with_cut() {
        let molecule = Molecule::build(&"[H][C@]1(N)OC1").unwrap();

        assert_eq!(molecule.atom_parity(&1), Ok(Some(Parity::Negative)));
    }

    #[test]
    fn atom_parity_given_counterclockwise_trailing_without_hydrogen_with_distal_cut() {
        let molecule = Molecule::build(&"[H][C@](N)1OC1").unwrap();

        assert_eq!(molecule.atom_parity(&1), Ok(Some(Parity::Negative)));
    }

    #[test]
    fn bond_order_given_unknown_sid() {
        let molecule = Molecule::build(&"C").unwrap();

        assert_eq!(molecule.bond_order(&1, &0), Err(GraphError::UnknownNode));
    }

    #[test]
    fn bond_order_given_unknown_tid() {
        let molecule = Molecule::build(&"C").unwrap();

        assert_eq!(molecule.bond_order(&0, &1), Err(GraphError::UnknownNode));
    }

    #[test]
    fn bond_order_given_dot() {
        let molecule = Molecule::build(&"C.C").unwrap();

        assert_eq!(molecule.bond_order(&0, &1), Ok(0));
    }

    #[test]
    fn bond_order_given_none() {
        let molecule = Molecule::build(&"CC").unwrap();

        assert_eq!(molecule.bond_order(&0, &1), Ok(1));
    }

    #[test]
    fn bond_order_given_single() {
        let molecule = Molecule::build(&"C-C").unwrap();

        assert_eq!(molecule.bond_order(&0, &1), Ok(1));
    }

    #[test]
    fn bond_order_given_double() {
        let molecule = Molecule::build(&"C=C").unwrap();

        assert_eq!(molecule.bond_order(&0, &1), Ok(2));
    }

    #[test]
    fn bond_order_given_up() {
        let molecule = Molecule::build(&"C/C").unwrap();

        assert_eq!(molecule.bond_order(&0, &1), Ok(1));
    }

    #[test]
    fn bond_order_given_down() {
        let molecule = Molecule::build(&"C\\C").unwrap();

        assert_eq!(molecule.bond_order(&0, &1), Ok(1));
    }

    #[test]
    fn bond_parity_given_unknown_sid() {
        let molecule = Molecule::build(&"C").unwrap();

        assert_eq!(molecule.bond_parity(&1, &0), Err(GraphError::UnknownNode));
    }

    #[test]
    fn bond_parity_given_unknown_tid() {
        let molecule = Molecule::build(&"C").unwrap();

        assert_eq!(molecule.bond_parity(&0, &1), Err(GraphError::UnknownNode));
    }

    #[test]
    fn bond_parity_given_dot() {
        let molecule = Molecule::build(&"C.C").unwrap();

        assert_eq!(molecule.bond_parity(&0, &1), Ok(None));
    }

    #[test]
    fn bond_parity_given_single() {
        let molecule = Molecule::build(&"C-C").unwrap();

        assert_eq!(molecule.bond_parity(&0, &1), Ok(None));
    }

    #[test]
    fn bond_parity_given_double() {
        let molecule = Molecule::build(&"C=C").unwrap();

        assert_eq!(molecule.bond_parity(&0, &1), Ok(None));
    }

    #[test]
    fn bond_parity_given_double_with_bl_tl() {
        let molecule = Molecule::build(&"C/N(\\C)=N").unwrap();

        assert_eq!(molecule.bond_parity(&1, &3), Ok(None));
    }

    #[test]
    fn bond_parity_given_double_with_br_tr() {
        let molecule = Molecule::build(&"C=N(/C)\\N").unwrap();

        assert_eq!(molecule.bond_parity(&0, &1), Ok(None));
    }

    #[test]
    fn bond_parity_given_double_with_bl_tr() {
        let molecule = Molecule::build(&"C/N=N/C").unwrap();

        assert_eq!(molecule.bond_parity(&1, &2), Ok(Some(Parity::Negative)));
    }

    #[test]
    fn bond_parity_given_double_with_tl_br() {
        let molecule = Molecule::build(&"C\\N=N\\C").unwrap();

        assert_eq!(molecule.bond_parity(&1, &2), Ok(Some(Parity::Negative)));
    }

    #[test]
    fn bond_parity_given_double_with_bl_br() {
        let molecule = Molecule::build(&"C/N=N\\C").unwrap();

        assert_eq!(molecule.bond_parity(&1, &2), Ok(Some(Parity::Positive)));
    }

    #[test]
    fn bond_parity_given_double_with_tl_tr() {
        let molecule = Molecule::build(&"C\\N=N/C").unwrap();

        assert_eq!(molecule.bond_parity(&1, &2), Ok(Some(Parity::Positive)));
    }

    #[test]
    fn bond_parity_given_double_with_bl_tl_tr() {
        let molecule = Molecule::build(&"C/C(/C)=C/C").unwrap();

        assert_eq!(molecule.bond_parity(&1, &3), Ok(Some(Parity::Negative)));
    }

    #[test]
    fn bond_parity_given_double_with_bl_tl_br() {
        let molecule = Molecule::build(&"C/C(/C)=C\\C").unwrap();

        assert_eq!(molecule.bond_parity(&1, &3), Ok(Some(Parity::Positive)));
    }

    #[test]
    fn bond_parity_given_double_with_bl_tl_tl_br() {
        let molecule = Molecule::build(&"C/C(/C)=C(/C)\\C").unwrap();

        assert_eq!(molecule.bond_parity(&1, &3), Ok(Some(Parity::Negative)));
    }

    #[test]
    fn is_empty_given_methane() {
        let molecule = Molecule::build(&"C").unwrap();

        assert_eq!(molecule.is_empty(), false);
    }

    #[test]
    fn order_given_ethane() {
        let molecule = Molecule::build(&"CC").unwrap();

        assert_eq!(molecule.order(), 2);
    }

    #[test]
    fn size_given_propane() {
        let molecule = Molecule::build(&"CCC").unwrap();

        assert_eq!(molecule.size(), 2);
    }

    #[test]
    fn nodes_given_ethane() {
        let molecule = Molecule::build(&"CC").unwrap();
        let nodes = molecule.nodes().collect::<Vec<_>>();

        assert_eq!(nodes, vec![ &0, &1 ]);
    }

    #[test]
    fn edges_given_propane() {
        let molecule = Molecule::build(&"CCC").unwrap();
        let edges = molecule.edges().collect::<Vec<_>>();

        assert_eq!(edges, vec![ ( &0, &1 ), ( &1, &2 ) ]);
    }

    #[test]
    fn has_node_given_unknown() {
        let molecule = Molecule::build(&"C").unwrap();

        assert_eq!(molecule.has_node(&1), false);
    }

    #[test]
    fn has_node_given_known() {
        let molecule = Molecule::build(&"C").unwrap();

        assert_eq!(molecule.has_node(&0), true);
    }

    #[test]
    fn neighbors_given_unknown() {
        let molecule = Molecule::build(&"CCC").unwrap();
        let neighbors = molecule.neighbors(&3);

        assert_eq!(neighbors.err(), Some(GraphError::UnknownNode));
    }

    #[test]
    fn neighbors_given_known() {
        let molecule = Molecule::build(&"CCC").unwrap();
        let neighbors = molecule.neighbors(&1).unwrap().collect::<Vec<_>>();

        assert_eq!(neighbors, vec![ &0, &2 ]);
    }

    #[test]
    fn degree_given_unknown() {
        let molecule = Molecule::build(&"C").unwrap();
        let neighbors = molecule.degree(&1);

        assert_eq!(neighbors.err(), Some(GraphError::UnknownNode));
    }

    #[test]
    fn degree_given_methane() {
        let molecule = Molecule::build(&"C").unwrap();

        assert_eq!(molecule.degree(&0), Ok(0));
    }

    #[test]
    fn degree_given_propane() {
        let molecule = Molecule::build(&"C(C)C").unwrap();

        assert_eq!(molecule.degree(&0), Ok(2));
    }
    
    #[test]
    fn has_edge_given_unknown_source() {
        let molecule = Molecule::build(&"C").unwrap();

        assert_eq!(molecule.has_edge(&0, &1), Err(GraphError::UnknownNode));
    }
    
    #[test]
    fn has_edge_given_unknown_target() {
        let molecule = Molecule::build(&"C").unwrap();

        assert_eq!(molecule.has_edge(&1, &0), Err(GraphError::UnknownNode));
    }

    #[test]
    fn has_edge_given_none() {
        let molecule = Molecule::build(&"C.C").unwrap();

        assert_eq!(molecule.has_edge(&0, &1).unwrap(), false);
    }

    #[test]
    fn has_edge_given_edge() {
        let molecule = Molecule::build(&"CC").unwrap();

        assert_eq!(molecule.has_edge(&0, &1).unwrap(), true);
    }
}