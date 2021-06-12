use purr::graph;
use purr::parts::BondKind;

use super::{trigonal_parity, Error};
use crate::molecule::Bond;

pub fn to_bond(
    sid: usize,
    bond: &graph::Bond,
    atoms: &[graph::Atom],
    trace: &[usize],
) -> Result<Bond, Error> {
    let (electrons, parity) = match &bond.kind {
        BondKind::Elided | BondKind::Aromatic | BondKind::Single => (2, None),
        BondKind::Up | BondKind::Down => {
            if has_double(sid, bond, atoms) {
                (2, None)
            } else {
                return Err(Error::BondKind(trace[bond.tid]));
            }
        }
        BondKind::Double => {
            let left_parity = match trigonal_parity(&atoms[sid].bonds) {
                Ok(parity) => parity,
                Err(()) => return Err(Error::BondKind(trace[sid])),
            };
            let right_parity = match trigonal_parity(&atoms[bond.tid].bonds) {
                Ok(parity) => parity,
                Err(()) => return Err(Error::BondKind(trace[bond.tid])),
            };

            if let Some(left) = left_parity {
                if let Some(right) = &right_parity {
                    let neg_right = right.negate();

                    (4, Some(left.multiply(&neg_right)))
                } else {
                    (4, None)
                }
            } else {
                (4, None)
            }
        }
        BondKind::Triple => (6, None),
        BondKind::Quadruple => (8, None),
    };

    Ok(Bond {
        electrons,
        parity,
        tid: bond.tid,
    })
}

fn has_double(sid: usize, bond: &graph::Bond, atoms: &[graph::Atom]) -> bool {
    let source = &atoms[sid];

    if source
        .bonds
        .iter()
        .any(|bond| bond.kind == BondKind::Double)
    {
        return true;
    }

    let target = &atoms[bond.tid];

    target
        .bonds
        .iter()
        .any(|bond| bond.kind == BondKind::Double)
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use purr::graph::{from_tree, Bond as PurrBond};
    use purr::read::{read, Reading};

    use super::*;
    use crate::molecule::Parity;

    #[test]
    fn elided() {
        let Reading { root, trace } = read("CC").unwrap();
        let atoms = from_tree(root).unwrap();
        let input = PurrBond::new(BondKind::Elided, 1);
        let bond = to_bond(0, &input, &atoms, &trace);

        assert_eq!(bond, Ok(Bond::new(2, None, 1)))
    }

    #[test]
    fn single() {
        let Reading { root, trace } = read("C-C").unwrap();
        let atoms = from_tree(root).unwrap();
        let input = PurrBond::new(BondKind::Single, 1);
        let bond = to_bond(0, &input, &atoms, &trace);

        assert_eq!(bond, Ok(Bond::new(2, None, 1)))
    }

    #[test]
    fn double() {
        let Reading { root, trace } = read("C=C").unwrap();
        let atoms = from_tree(root).unwrap();
        let input = PurrBond::new(BondKind::Double, 1);
        let bond = to_bond(0, &input, &atoms, &trace);

        assert_eq!(bond, Ok(Bond::new(4, None, 1)))
    }

    #[test]
    fn double_up() {
        let Reading { root, trace } = read("C=C/C").unwrap();
        let atoms = from_tree(root).unwrap();
        let input = PurrBond::new(BondKind::Double, 1);
        let bond = to_bond(0, &input, &atoms, &trace);

        assert_eq!(bond, Ok(Bond::new(4, None, 1)))
    }

    #[test]
    fn double_up_down() {
        let Reading { root, trace } = read("C=C(/C)\\C").unwrap();
        let atoms = from_tree(root).unwrap();
        let input = PurrBond::new(BondKind::Double, 1);
        let bond = to_bond(0, &input, &atoms, &trace);

        assert_eq!(bond, Ok(Bond::new(4, None, 1)))
    }

    #[test]
    fn double_up_up_elided() {
        let Reading { root, trace } = read("C=P(/C)(/C)C").unwrap();
        let atoms = from_tree(root).unwrap();
        let input = PurrBond::new(BondKind::Double, 1);
        let bond = to_bond(0, &input, &atoms, &trace);

        assert_eq!(bond, Ok(Bond::new(4, None, 1)))
    }

    #[test]
    fn up_double_up() {
        let Reading { root, trace } = read("C/C=C/C").unwrap();
        let atoms = from_tree(root).unwrap();
        let input = PurrBond::new(BondKind::Double, 2);
        let bond = to_bond(1, &input, &atoms, &trace);

        assert_eq!(bond, Ok(Bond::new(4, Some(Parity::Negative), 2)))
    }

    #[test]
    fn up_double_down() {
        let Reading { root, trace } = read("C/C=C\\C").unwrap();
        let atoms = from_tree(root).unwrap();
        let input = PurrBond::new(BondKind::Double, 2);
        let bond = to_bond(1, &input, &atoms, &trace);

        assert_eq!(bond, Ok(Bond::new(4, Some(Parity::Positive), 2)))
    }

    #[test]
    fn up_up_double_down() {
        let Reading { root, trace } = read("C/C(/C)=C\\C").unwrap();
        let atoms = from_tree(root).unwrap();
        let input = PurrBond::new(BondKind::Double, 3);
        let bond = to_bond(1, &input, &atoms, &trace);

        assert_eq!(bond, Ok(Bond::new(4, Some(Parity::Negative), 3)))
    }

    #[test]
    fn up_up_double_up() {
        let Reading { root, trace } = read("C/C(/C)=C/C").unwrap();
        let atoms = from_tree(root).unwrap();
        let input = PurrBond::new(BondKind::Double, 3);
        let bond = to_bond(1, &input, &atoms, &trace);

        assert_eq!(bond, Ok(Bond::new(4, Some(Parity::Positive), 3)))
    }

    #[test]
    fn up() {
        let Reading { root, trace } = read("C/C=C/C").unwrap();
        let atoms = from_tree(root).unwrap();
        let input = PurrBond::new(BondKind::Double, 1);
        let bond = to_bond(0, &input, &atoms, &trace);

        assert_eq!(bond, Ok(Bond::new(4, None, 1)))
    }

    #[test]
    fn down() {
        let Reading { root, trace } = read("C\\C=C/C").unwrap();
        let atoms = from_tree(root).unwrap();
        let input = PurrBond::new(BondKind::Double, 1);
        let bond = to_bond(0, &input, &atoms, &trace);

        assert_eq!(bond, Ok(Bond::new(4, None, 1)))
    }

    #[test]
    fn triple() {
        let Reading { root, trace } = read("C#C").unwrap();
        let atoms = from_tree(root).unwrap();
        let input = PurrBond::new(BondKind::Triple, 1);
        let bond = to_bond(0, &input, &atoms, &trace);

        assert_eq!(bond, Ok(Bond::new(6, None, 1)))
    }

    #[test]
    fn quadruple() {
        let Reading { root, trace } = read("C#C").unwrap();
        let atoms = from_tree(root).unwrap();
        let input = PurrBond::new(BondKind::Quadruple, 1);
        let bond = to_bond(0, &input, &atoms, &trace);

        assert_eq!(bond, Ok(Bond::new(8, None, 1)))
    }
}
