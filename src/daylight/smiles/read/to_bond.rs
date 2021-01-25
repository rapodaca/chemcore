use purr::{ graph };
use purr::parts::BondKind;

use crate::molecule::{ Bond, Parity };
use super::Error;

pub fn to_bond(
    sid: usize,
    bond: &graph::Bond,
    atoms: &[graph::Atom],
    trace: &[usize]
) -> Result<Bond, Error> {
    let (electrons, parity) = match &bond.kind {
        BondKind::Elided |
        BondKind::Aromatic |
        BondKind::Single => (2, None),
        BondKind::Up |
        BondKind::Down => if has_double(sid, bond, atoms) {
            (2, None)
        } else {
            return Err(Error::BondKind(trace[bond.tid]))
        },
        BondKind::Double => {
            let left_parity = trigonal_parity(sid, atoms, trace)?;
            let right_parity = trigonal_parity(bond.tid, atoms, trace)?;

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
        },
        BondKind::Triple => (6, None),
        BondKind::Quadruple => (8, None)
    };

    Ok(Bond {
        electrons,
        parity,
        tid: bond.tid
    })
}

fn has_double(sid: usize, bond: &graph::Bond, atoms: &[graph::Atom]) -> bool {
    let source = &atoms[sid];

    if source.bonds.iter().any(|bond| bond.kind == BondKind::Double) {
        return true
    }

    let target = &atoms[bond.tid];

    target.bonds.iter().any(|bond| bond.kind == BondKind::Double)
}

fn trigonal_parity(
    id: usize, atoms: &[graph::Atom], trace: &[usize]
) -> Result<Option<Parity>, Error> {
    let bonds = &atoms[id].bonds;
    let first = match bonds.get(0) {
        Some(first) => first,
        None => return Ok(None)
    };
    let second = match bonds.get(1) {
        Some(second) => second,
        None => return Ok(None)
    };
    let third = bonds.get(2);

    match first.kind {
        BondKind::Up => match second.kind {
            BondKind::Up => Err(Error::BondKind(trace[second.tid])),
            BondKind::Down => match third {
                Some(third) => match third.kind {
                    BondKind::Double =>  if bonds.len() == 3 {
                        Ok(Some(Parity::Positive))
                    } else {
                        Err(Error::BondKind(trace[first.tid]))
                    },
                    _ => Err(Error::BondKind(trace[first.tid]))
                },
                None => unreachable!()
            },
            BondKind::Double => match third {
                Some(third) => match third.kind {
                    BondKind::Down |
                    BondKind::Elided => if bonds.len() == 3 {
                        Ok(Some(Parity::Positive))
                    } else {
                        Err(Error::BondKind(trace[first.tid]))
                    },
                    _ => Err(Error::BondKind(trace[third.tid]))
                },
                None => Ok(Some(Parity::Positive))
            },
            BondKind::Elided => match third {
                Some(third) => match third.kind {
                    BondKind::Double => if bonds.len() == 3 {
                        Ok(Some(Parity::Positive))
                    } else {
                        Err(Error::BondKind(trace[first.tid]))
                    },
                    _ => Err(Error::BondKind(trace[third.tid]))
                },
                None => unreachable!()
            },
            _ => Err(Error::BondKind(trace[second.tid]))
        },
        BondKind::Down => match second.kind {
            BondKind::Down => Err(Error::BondKind(trace[second.tid])),
            BondKind::Up => match third {
                Some(third) => match third.kind {
                    BondKind::Double => if bonds.len() == 3 {
                        Ok(Some(Parity::Negative))
                    } else {
                        Err(Error::BondKind(trace[first.tid]))
                    },
                    _ => Err(Error::BondKind(trace[first.tid]))
                },
                None => unreachable!()
            },
            BondKind::Double => match third {
                Some(third) => match third.kind {
                    BondKind::Up |
                    BondKind::Elided => if bonds.len() == 3 {
                        Ok(Some(Parity::Negative))
                    } else {
                        Err(Error::BondKind(trace[first.tid]))
                    },
                    _ => Err(Error::BondKind(trace[third.tid]))
                },
                None => Ok(Some(Parity::Negative))
            },
            BondKind::Elided => match third {
                Some(third) => match third.kind {
                    BondKind::Double => if bonds.len() == 3 {
                        Ok(Some(Parity::Negative))
                    } else {
                        Err(Error::BondKind(trace[first.tid]))
                    },
                    _ => Err(Error::BondKind(trace[third.tid]))
                },
                None => unreachable!()
            },
            _ => Err(Error::BondKind(trace[second.tid]))

        },
        BondKind::Double => match second.kind {
            BondKind::Up => match third {
                Some(third) => match third.kind {
                    BondKind::Up => Err(Error::BondKind(trace[third.tid])),
                    BondKind::Down => if bonds.len() == 3 {
                        Ok(Some(Parity::Negative))
                    } else {
                        Err(Error::BondKind(trace[second.tid]))
                    },
                    BondKind::Elided => if bonds.len() == 3 {
                        Ok(Some(Parity::Negative))
                    } else {
                        Err(Error::BondKind(trace[bonds[3].tid]))
                    }
                    _ => Err(Error::BondKind(trace[third.tid]))
                },
                None => Ok(Some(Parity::Negative))
            },
            BondKind::Down => match third {
                Some(third) => match third.kind {
                    BondKind::Down => Err(Error::BondKind(trace[third.tid])),
                    BondKind::Up => if bonds.len() == 3 {
                        Ok(Some(Parity::Positive))
                    } else {
                        Err(Error::BondKind(trace[second.tid]))
                    },
                    BondKind::Elided => if bonds.len() == 3 {
                        Ok(Some(Parity::Positive))
                    } else {
                        Err(Error::BondKind(trace[bonds[3].tid]))
                    },
                    _ => Err(Error::BondKind(trace[third.tid]))
                },
                None => Ok(Some(Parity::Positive))
            }
            _ => match third {
                Some(third) => match third.kind {
                    BondKind::Up |
                    BondKind::Down => Err(Error::BondKind(trace[second.tid])),
                    _ => Ok(None)
                },
                None => Ok(None)
            }
        }
        BondKind::Elided => match second.kind {
            BondKind::Double => match third {
                Some(third) => match third.kind {
                    BondKind::Up => Ok(Some(Parity::Negative)),
                    BondKind::Down => Ok(Some(Parity::Positive)),
                    _ => Ok(None)
                },
                None => Ok(None)
            },
            BondKind::Up => Ok(Some(Parity::Negative)),
            BondKind::Down => Ok(Some(Parity::Positive)),
            _ => Ok(None)
        }
        _ => match second.kind {
            BondKind::Up |
            BondKind::Down => Err(Error::BondKind(trace[first.tid])),
            _ => Ok(None)
        }
    }
}