use purr::mol::{ Atom, Style };

use crate::molecule::Parity;
use crate::molecule::Error;

pub fn bond_parity(
    sid: usize, tid: usize, atoms: &[Atom]
) -> Result<Option<Parity>, Error> {
    let left = score_half_bond(sid, tid, atoms)?;
    let right = score_half_bond(tid, sid, atoms)?;

    if let Some(p1) = left {
        if let Some(p2) = right {
            Ok(Some(p1.multiply(&p2)))
        } else {
            Ok(None)
        }
    } else {
        Ok(None)
    }
}

fn score_half_bond(
    sid: usize, tid: usize, atoms: &[Atom]
) -> Result<Option<Parity>, Error> {
    let mut up = None;
    let mut down = None;

    for bond in &atoms[sid].bonds {
        match bond.style {
            Some(Style::Up) => {
                if up.replace(bond).is_some() {
                    return Err(Error::InvalidConformation(sid, tid));
                }
            },
            Some(Style::Down) => {
                if down.replace(bond).is_some() {
                    return Err(Error::InvalidConformation(sid, tid));
                }
            },
            _ => ()
        }
    }

    if let Some(up) = up {
        if let Some(down) = down {
            if up.tid < down.tid {
                Ok(Some(Parity::Positive))
            } else {
                Ok(Some(Parity::Negative))
            }
        } else {
            Ok(Some(Parity::Positive))
        }
    } else if down.is_some() {
        Ok(Some(Parity::Negative))
    } else {
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use purr::read::read;

    #[test]
    fn tl_tl() {
        let atoms = read(&"C\\C(/C)=C").unwrap();
        let parity = bond_parity(1, 3, &atoms);

        assert_eq!(parity, Err(Error::InvalidConformation(1, 3)));
    }

    #[test]
    fn bl_bl() {
        let atoms = read(&"C/C(\\C)=C").unwrap();
        let parity = bond_parity(1, 3, &atoms);

        assert_eq!(parity, Err(Error::InvalidConformation(1, 3)));
    }

    #[test]
    fn tr_tr() {
        let atoms = read(&"C=C(/C)/C").unwrap();
        let parity = bond_parity(0, 1, &atoms);

        assert_eq!(parity, Err(Error::InvalidConformation(1, 0)));
    }

    #[test]
    fn br_br() {
        let atoms = read(&"C=C(\\C)\\C").unwrap();
        let parity = bond_parity(0, 1, &atoms);

        assert_eq!(parity, Err(Error::InvalidConformation(1, 0)));
    }
    
    #[test]
    fn tl_tr() {
        let atoms = read(&"C\\C=C/C").unwrap();
        let parity = bond_parity(1, 2, &atoms);

        assert_eq!(parity, Ok(Some(Parity::Positive)));
    }

    #[test]
    fn bl_tr() {
        let atoms = read(&"C/C=C/C").unwrap();
        let parity = bond_parity(1, 2, &atoms);

        assert_eq!(parity, Ok(Some(Parity::Negative)));
    }

    #[test]
    fn tl_br() {
        let atoms = read(&"C\\C=C\\C").unwrap();
        let parity = bond_parity(1, 2, &atoms);

        assert_eq!(parity, Ok(Some(Parity::Negative)));
    }

    #[test]
    fn bl_br() {
        let atoms = read(&"C/C=C\\C").unwrap();
        let parity = bond_parity(1, 2, &atoms);

        assert_eq!(parity, Ok(Some(Parity::Positive)));
    }
}