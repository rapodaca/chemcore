use super::spec;
use super::element::Element;
use super::error::Error;
use super::parity::Parity;
use super::bond::Bond;
use super::bond_order::BondOrder;

/// Exposes methods to be accessed privately by DefaultMolecule.
pub struct Atom {
    element: Element,
    nonbonding_electrons: u8,
    bonding_electrons: u8,
    hydrogens: u8,
    isotope: Option<u16>,
    parity: Option<Parity>,
    bonds: Vec<Bond>,
    neighbors: Vec<usize>
}

impl<'a> Atom {
    pub fn build(spec: spec::Atom) -> Result<Self, Error> {
        let element = spec.element;
        let isotope = spec.isotope;

        if let Some(mass) = isotope {
            if mass < Element::atomic_number(&element) {
                return Err(Error::ImpossibleIsotope);
            }
        }

        let mut nonbonding_electrons =
            Element::valence_electrons(&element) as i16;
        let mut bonding_electrons = 0i16;
        let hydrogens = spec.hydrogens as i16;
        let ion = spec.ion as i16;

        if ion > nonbonding_electrons {
            return Err(Error::HypervalentAtom);
        } else {
            nonbonding_electrons -= ion;
        }

        if hydrogens > nonbonding_electrons {
            return Err(Error::HypervalentAtom);
        } else {
            nonbonding_electrons -= hydrogens;
            bonding_electrons += hydrogens;
        }

        Ok(Atom {
            element,
            nonbonding_electrons: nonbonding_electrons as u8,
            bonding_electrons: bonding_electrons as u8,
            hydrogens: hydrogens as u8,
            isotope,
            parity: spec.parity,
            bonds: Vec::new(),
            neighbors: Vec::new()
        })
    }

    pub fn element(&self) -> Element {
        self.element
    }

    pub fn isotope(&self) -> Option<u16> {
        self.isotope
    }

    pub fn electrons(&self) -> u8 {
        self.nonbonding_electrons
    }

    pub fn hydrogens(&self) -> u8 {
        self.hydrogens
    }

    pub fn charge(&self) -> i8 {
        let mut result = Element::valence_electrons(&self.element) as i16;

        result -= self.nonbonding_electrons as i16;
        result -= self.bonding_electrons as i16;

        result as i8
    }

    pub fn bond_order(&self, tid: &usize) -> u8 {
        match self.bonds.iter().find(|bond| bond.tid == *tid) {
            Some(bond) => bond.multiplicity(),
            None => 0
        }
    }

    pub fn bond_parity(&self, tid: &usize) -> Option<Parity> {
        match self.bonds.iter().find(|bond| bond.tid == *tid) {
            Some(bond) => bond.parity,
            None => None
        }
    }

    pub fn atom_parity(&self) -> Option<Parity> {
        self.parity
    }

    pub fn degree(&self) -> usize {
        self.neighbors.len()
    }

    pub fn neighbors(&'a self) -> std::slice::Iter<'a, usize> {
        self.neighbors.iter()
    }

    pub fn has_edge(&self, tid: usize) -> bool {
        self.bonds.iter().any(|bond| bond.tid == tid)
    }

    pub fn add_bond(
        &mut self, tid: usize, order: BondOrder, parity: Option<Parity>
    ) -> Result<(), Error> {
        if self.neighbors.contains(&tid) {
            return Err(Error::DuplicateBond);
        }

        let bond = Bond { tid, order, parity };
        let nonbonding_electrons =
            self.nonbonding_electrons as i16 -
            bond.multiplicity() as i16;
        let bonding_electrons =
            self.bonding_electrons as i16 +
            bond.multiplicity() as i16;

        if nonbonding_electrons < 0 {
            Err(Error::HypervalentAtom)   
        } else {
            self.bonds.push(bond);
            self.neighbors.push(tid);
            self.nonbonding_electrons = nonbonding_electrons as u8;
            self.bonding_electrons = bonding_electrons as u8;

            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn returns_error_given_overionized() {
        let atom = Atom::build(spec::Atom {
            element: Element::H, hydrogens: 0, ion: 2, isotope: None,
            parity: None
        });

        assert_eq!(atom.err(), Some(Error::HypervalentAtom));
    }

    #[test]
    fn returns_error_given_oversaturated() {
        let atom = Atom::build(spec::Atom {
            element: Element:: H, hydrogens: 2, ion: 0, isotope: None,
            parity: None
        });

        assert_eq!(atom.err(), Some(Error::HypervalentAtom));
    }

    #[test]
    fn returns_error_given_impossible_isotope() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 4, ion: 0, isotope: Some(5),
            parity: None
        });

        assert_eq!(atom.err(), Some(Error::ImpossibleIsotope));
    }

    #[test]
    fn element_given_hydrogen() {
        let atom = Atom::build(spec::Atom {
            element: Element:: H, ..Default::default()
        }).unwrap();

        assert_eq!(atom.element(), Element::H);
    }

    #[test]
    fn isotope_given_hydrogen() {
        let atom = Atom::build(spec::Atom {
            element: Element:: H, ..Default::default()
        }).unwrap();

        assert_eq!(atom.isotope(), None);
    }

    #[test]
    fn isotope_given_deuterium() {
        let atom = Atom::build(spec::Atom {
            element: Element:: H, hydrogens: 0, ion: 0, isotope: Some(2),
            parity: None
        }).unwrap();

        assert_eq!(atom.isotope(), Some(2));
    }

    #[test]
    fn electrons_given_carbon() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, ..Default::default()
        }).unwrap();

        assert_eq!(atom.electrons(), 4);
    }

    #[test]
    fn electrons_given_methyl_radical() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 3, ion: 0, isotope: None,
            parity: None
        }).unwrap();

        assert_eq!(atom.electrons(), 1);
    }

    #[test]
    fn electrons_given_methane() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 4, ion: 0, isotope: None,
            parity: None
        }).unwrap();

        assert_eq!(atom.electrons(), 0);
    }

    #[test]
    fn electrons_given_methyl_anion() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 3, ion: -1, isotope: None,
            parity: None
        }).unwrap();

        assert_eq!(atom.electrons(), 2);
    }

    #[test]
    fn electrons_given_ethane() {
        let mut atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 3, ion: 0, isotope: None,
            parity: None
        }).unwrap();

        assert_eq!(atom.add_bond(1, BondOrder::Single, None), Ok(()));
        assert_eq!(atom.electrons(), 0);
    }

    #[test]
    fn electrons_given_ethene() {
        let mut atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 2, ion: 0, isotope: None,
            parity: None
        }).unwrap();

        assert_eq!(atom.add_bond(1, BondOrder::Double, None), Ok(()));
        assert_eq!(atom.electrons(), 0);
    }

    #[test]
    fn hydrogens_given_carbon() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, ..Default::default()
        }).unwrap();

        assert_eq!(atom.hydrogens(), 0);
    }

    #[test]
    fn hydrogens_given_methyl_radical_and_single() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 3, ion: 0, isotope: None,
            parity: None
        }).unwrap();

        assert_eq!(atom.hydrogens(), 3);
    }

    #[test]
    fn hydrogens_given_methane() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 4, ion: 0, isotope: None,
            parity: None
        }).unwrap();

        assert_eq!(atom.hydrogens(), 4);
    }

    #[test]
    fn add_bond_given_methyl_radical_and_double() {
        let mut atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 3, ion: 0, isotope: None,
            parity: None
        }).unwrap();
        let err = atom.add_bond(1, BondOrder::Double, None).err();
    
        assert_eq!(err, Some(Error::HypervalentAtom));
    }

    #[test]
    fn add_bond_given_methyl_radical_and_duplicate_tid() {
        let mut atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 2, ion: 0, isotope: None,
            parity: None
        }).unwrap();

        assert_eq!(atom.add_bond(1, BondOrder::Single, None), Ok(()));
        assert_eq!(
            atom.add_bond(1, BondOrder::Single, None),
            Err(Error::DuplicateBond)
        );
    }

    #[test]
    fn charge_given_carbon() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, ..Default::default()
        }).unwrap();

        assert_eq!(atom.charge(), 0);
    }

    #[test]
    fn charge_given_methyl_radical() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 3, ..Default::default()
        }).unwrap();

        assert_eq!(atom.charge(), 0);
    }

    #[test]
    fn charge_given_methane() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 4, ..Default::default()
        }).unwrap();

        assert_eq!(atom.charge(), 0);
    }

    #[test]
    fn charge_given_methyl_cation() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 3, ion: 1, isotope: None,
            parity: None
        }).unwrap();

        assert_eq!(atom.charge(), 1);
    }

    #[test]
    fn charge_given_methyl_anion() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 3, ion: -1, isotope: None,
            parity: None
        }).unwrap();

        assert_eq!(atom.charge(), -1);
    }

    #[test]
    fn charge_given_hydroborate() {
        let atom = Atom::build(spec::Atom {
            element: Element::B, hydrogens: 4, ion: -1, isotope: None,
            parity: None
        }).unwrap();

        assert_eq!(atom.charge(), -1);
    }

    #[test]
    fn charge_given_triple_bond() {
        let mut atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 1, ..Default::default()
        }).unwrap();

        assert_eq!(atom.add_bond(1, BondOrder::Triple, None), Ok(()));
        assert_eq!(atom.charge(), 0);
    }

    #[test]
    fn bond_order_given_methane() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 3, ..Default::default()
        }).unwrap();

        assert_eq!(atom.bond_order(&1), 0);
    }

    #[test]
    fn bond_order_given_ethane() {
        let mut atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 3, ..Default::default()
        }).unwrap();

        assert_eq!(atom.add_bond(1, BondOrder::Single, None), Ok(()));
        assert_eq!(atom.bond_order(&1), 1);
    }

    #[test]
    fn bond_order_given_ethene() {
        let mut atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 2, ..Default::default()
        }).unwrap();

        assert_eq!(atom.add_bond(1, BondOrder::Double, None), Ok(()));
        assert_eq!(atom.bond_order(&1), 2);
    }

    #[test]
    fn bond_parity_given_methane() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 3, ion: 0, isotope: None,
            parity: None
        }).unwrap();

        assert_eq!(atom.bond_parity(&1), None);
    }

    #[test]
    fn bond_parity_given_ethane_with_parity() {
        let mut atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 3, ..Default::default()
        }).unwrap();

        assert_eq!(
            atom.add_bond(1, BondOrder::Single, Some(Parity::Positive)),
            Ok(())
        );
        assert_eq!(atom.bond_parity(&1), Some(Parity::Positive));
    }

    #[test]
    fn atom_parity_given_methane() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 2, ..Default::default()
        }).unwrap();

        assert_eq!(atom.atom_parity(), None);
    }

    #[test]
    fn atom_parity_given_methane_with_parity() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 2, ion: 0, isotope: None,
            parity: Some(Parity::Positive)
        }).unwrap();

        assert_eq!(atom.atom_parity(), Some(Parity::Positive));
    }

    #[test]
    fn degree_given_methane() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 4, ..Default::default()
        }).unwrap();

        assert_eq!(atom.degree(), 0);
    }

    #[test]
    fn degree_given_ethane() {
        let mut atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 3, ..Default::default()
        }).unwrap();

        assert_eq!(atom.add_bond(1, BondOrder::Single, None), Ok(()));
        assert_eq!(atom.degree(), 1);
    }

    #[test]
    fn neighbors_given_methane() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 4, ..Default::default()
        }).unwrap();
        let neighbors = atom.neighbors().collect::<Vec<_>>();

        assert_eq!(neighbors.is_empty(),true);
    }

    #[test]
    fn neighbors_given_propane() {
        let mut atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 2, ..Default::default()
        }).unwrap();

        assert_eq!(atom.add_bond(1, BondOrder::Single, None), Ok(()));
        assert_eq!(atom.add_bond(2, BondOrder::Single, None), Ok(()));

        let neighbors = atom.neighbors().collect::<Vec<_>>();

        assert_eq!(neighbors, vec![ &1, &2 ]);
    }

    #[test]
    fn has_edge_given_methane() {
        let atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 2, ..Default::default()
        }).unwrap();

        assert_eq!(atom.has_edge(1), false);
    }

    #[test]
    fn has_edge_given_ethane() {
        let mut atom = Atom::build(spec::Atom {
            element: Element::C, hydrogens: 2, ..Default::default()
        }).unwrap();

        assert_eq!(atom.add_bond(1, BondOrder::Single, None), Ok(()));
        assert_eq!(atom.has_edge(1), true);
    }
}