use purr::graph::Bond;
use purr::parts::BondKind;

use crate::molecule::Parity;

fn is_directional(kind: &BondKind) -> bool {
    kind == &BondKind::Up || kind == &BondKind::Down
}

pub fn trigonal_parity(bonds: &[Bond]) -> Result<Option<Parity>, ()> {
    if bonds.len() > 3 {
        return Ok(None);
    }

    let first = match bonds.get(0) {
        Some(first) => &first.kind,
        None => return Ok(None),
    };
    let second = match bonds.get(1) {
        Some(second) => &second.kind,
        None => return Ok(None),
    };
    let third = match bonds.get(2) {
        Some(third) => Some(&third.kind),
        None => None,
    };

    if !is_directional(first) && !is_directional(second) {
        if let Some(third) = third {
            if !is_directional(third) {
                return Ok(None);
            }
        } else {
            return Ok(None);
        }
    }

    match first {
        BondKind::Double => match second {
            BondKind::Elided => match third {
                Some(BondKind::Up) => Ok(Some(Parity::Positive)),
                Some(BondKind::Down) => Ok(Some(Parity::Negative)),
                _ => unreachable!(),
            },
            BondKind::Down => match third {
                None | Some(BondKind::Up) | Some(BondKind::Elided) => Ok(Some(Parity::Positive)),
                _ => Err(()),
            },
            BondKind::Up => match third {
                None | Some(BondKind::Down) | Some(BondKind::Elided) => Ok(Some(Parity::Negative)),
                _ => Err(()),
            },
            _ => Err(()),
        },
        BondKind::Up => match second {
            BondKind::Double => match third {
                None | Some(BondKind::Down) | Some(BondKind::Elided) => Ok(Some(Parity::Positive)),
                _ => Err(()),
            },
            BondKind::Down | BondKind::Elided => match third {
                Some(BondKind::Double) => Ok(Some(Parity::Negative)),
                _ => Err(()),
            },
            _ => Err(()),
        },
        BondKind::Down => match second {
            BondKind::Double => match third {
                None | Some(BondKind::Elided) | Some(BondKind::Up) => Ok(Some(Parity::Negative)),
                _ => Err(()),
            },
            BondKind::Up | BondKind::Elided => match third {
                Some(BondKind::Double) => Ok(Some(Parity::Positive)),
                _ => Err(()),
            },
            _ => Err(()),
        },
        BondKind::Elided => match second {
            BondKind::Up => match third {
                Some(BondKind::Double) => Ok(Some(Parity::Positive)),
                _ => Err(()),
            },
            BondKind::Down => match third {
                Some(BondKind::Double) => Ok(Some(Parity::Negative)),
                _ => Err(()),
            },
            BondKind::Double => match third {
                Some(BondKind::Up) => Ok(Some(Parity::Negative)),
                Some(BondKind::Down) => Ok(Some(Parity::Positive)),
                _ => Ok(None),
            },
            _ => Err(()),
        },
        _ => unreachable!(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn double_up() {
        let bonds = vec![Bond::new(BondKind::Double, 0), Bond::new(BondKind::Up, 1)];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Negative)))
    }

    #[test]
    fn double_elided_up() {
        let bonds = vec![
            Bond::new(BondKind::Double, 0),
            Bond::new(BondKind::Elided, 1),
            Bond::new(BondKind::Up, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Positive)))
    }

    #[test]
    fn double_elided_up_elided() {
        let bonds = vec![
            Bond::new(BondKind::Double, 0),
            Bond::new(BondKind::Elided, 1),
            Bond::new(BondKind::Up, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn double_down() {
        let bonds = vec![Bond::new(BondKind::Double, 0), Bond::new(BondKind::Down, 1)];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Positive)))
    }

    #[test]
    fn double_down_elided() {
        let bonds = vec![
            Bond::new(BondKind::Double, 0),
            Bond::new(BondKind::Down, 1),
            Bond::new(BondKind::Elided, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Positive)))
    }

    #[test]
    fn double_down_elided_elided() {
        let bonds = vec![
            Bond::new(BondKind::Double, 0),
            Bond::new(BondKind::Down, 1),
            Bond::new(BondKind::Elided, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn double_down_up() {
        let bonds = vec![
            Bond::new(BondKind::Double, 0),
            Bond::new(BondKind::Down, 1),
            Bond::new(BondKind::Up, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Positive)))
    }

    #[test]
    fn double_down_up_elided() {
        let bonds = vec![
            Bond::new(BondKind::Double, 0),
            Bond::new(BondKind::Down, 1),
            Bond::new(BondKind::Up, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn double_down_single() {
        let bonds = vec![
            Bond::new(BondKind::Double, 0),
            Bond::new(BondKind::Down, 1),
            Bond::new(BondKind::Single, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn double_elided_down() {
        let bonds = vec![
            Bond::new(BondKind::Double, 0),
            Bond::new(BondKind::Elided, 1),
            Bond::new(BondKind::Down, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Negative)))
    }

    #[test]
    fn double_elided_elided() {
        let bonds = vec![
            Bond::new(BondKind::Double, 0),
            Bond::new(BondKind::Elided, 1),
            Bond::new(BondKind::Elided, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn double_elided_down_elided() {
        let bonds = vec![
            Bond::new(BondKind::Double, 0),
            Bond::new(BondKind::Elided, 1),
            Bond::new(BondKind::Down, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn double_up_down() {
        let bonds = vec![
            Bond::new(BondKind::Double, 0),
            Bond::new(BondKind::Up, 1),
            Bond::new(BondKind::Down, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Negative)))
    }

    #[test]
    fn double_up_down_elided() {
        let bonds = vec![
            Bond::new(BondKind::Double, 0),
            Bond::new(BondKind::Up, 1),
            Bond::new(BondKind::Down, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn double_up_elided() {
        let bonds = vec![
            Bond::new(BondKind::Double, 0),
            Bond::new(BondKind::Up, 1),
            Bond::new(BondKind::Elided, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Negative)))
    }

    #[test]
    fn double_up_elided_elided() {
        let bonds = vec![
            Bond::new(BondKind::Double, 0),
            Bond::new(BondKind::Up, 1),
            Bond::new(BondKind::Elided, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn double_up_single() {
        let bonds = vec![
            Bond::new(BondKind::Double, 0),
            Bond::new(BondKind::Up, 1),
            Bond::new(BondKind::Single, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn double_single() {
        let bonds = vec![
            Bond::new(BondKind::Double, 0),
            Bond::new(BondKind::Single, 1),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn double_single_up() {
        let bonds = vec![
            Bond::new(BondKind::Double, 0),
            Bond::new(BondKind::Single, 1),
            Bond::new(BondKind::Up, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    // ----

    #[test]
    fn up_double() {
        let bonds = vec![Bond::new(BondKind::Up, 0), Bond::new(BondKind::Double, 1)];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Positive)))
    }

    #[test]
    fn up_double_elided() {
        let bonds = vec![
            Bond::new(BondKind::Up, 0),
            Bond::new(BondKind::Double, 1),
            Bond::new(BondKind::Elided, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Positive)))
    }

    #[test]
    fn up_double_elided_elided() {
        let bonds = vec![
            Bond::new(BondKind::Up, 0),
            Bond::new(BondKind::Double, 1),
            Bond::new(BondKind::Elided, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn up_double_single() {
        let bonds = vec![
            Bond::new(BondKind::Up, 0),
            Bond::new(BondKind::Double, 1),
            Bond::new(BondKind::Single, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn up_double_down() {
        let bonds = vec![
            Bond::new(BondKind::Up, 0),
            Bond::new(BondKind::Double, 1),
            Bond::new(BondKind::Down, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Positive)))
    }

    #[test]
    fn up_double_down_elided() {
        let bonds = vec![
            Bond::new(BondKind::Up, 0),
            Bond::new(BondKind::Double, 1),
            Bond::new(BondKind::Down, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn up_down() {
        let bonds = vec![Bond::new(BondKind::Up, 0), Bond::new(BondKind::Down, 1)];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn up_down_elided() {
        let bonds = vec![
            Bond::new(BondKind::Up, 0),
            Bond::new(BondKind::Down, 1),
            Bond::new(BondKind::Elided, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn up_down_double() {
        let bonds = vec![
            Bond::new(BondKind::Up, 0),
            Bond::new(BondKind::Down, 1),
            Bond::new(BondKind::Double, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Negative)))
    }

    #[test]
    fn up_down_double_elided() {
        let bonds = vec![
            Bond::new(BondKind::Up, 0),
            Bond::new(BondKind::Down, 1),
            Bond::new(BondKind::Double, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn up_down_single() {
        let bonds = vec![
            Bond::new(BondKind::Up, 0),
            Bond::new(BondKind::Down, 1),
            Bond::new(BondKind::Single, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn up_elided() {
        let bonds = vec![Bond::new(BondKind::Up, 0), Bond::new(BondKind::Elided, 1)];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn up_elided_double() {
        let bonds = vec![
            Bond::new(BondKind::Up, 0),
            Bond::new(BondKind::Elided, 1),
            Bond::new(BondKind::Double, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Negative)))
    }

    #[test]
    fn up_elided_double_elided() {
        let bonds = vec![
            Bond::new(BondKind::Up, 0),
            Bond::new(BondKind::Elided, 1),
            Bond::new(BondKind::Double, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn up_elided_single() {
        let bonds = vec![
            Bond::new(BondKind::Up, 0),
            Bond::new(BondKind::Elided, 1),
            Bond::new(BondKind::Single, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn up_single() {
        let bonds = vec![Bond::new(BondKind::Up, 0), Bond::new(BondKind::Single, 1)];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    // ----

    #[test]
    fn down_double() {
        let bonds = vec![Bond::new(BondKind::Down, 0), Bond::new(BondKind::Double, 1)];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Negative)))
    }

    #[test]
    fn down_double_elided() {
        let bonds = vec![
            Bond::new(BondKind::Down, 0),
            Bond::new(BondKind::Double, 1),
            Bond::new(BondKind::Elided, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Negative)))
    }

    #[test]
    fn down_double_elided_elided() {
        let bonds = vec![
            Bond::new(BondKind::Down, 0),
            Bond::new(BondKind::Double, 1),
            Bond::new(BondKind::Elided, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn down_double_up() {
        let bonds = vec![
            Bond::new(BondKind::Down, 0),
            Bond::new(BondKind::Double, 1),
            Bond::new(BondKind::Up, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Negative)))
    }

    #[test]
    fn down_double_up_elided() {
        let bonds = vec![
            Bond::new(BondKind::Down, 0),
            Bond::new(BondKind::Double, 1),
            Bond::new(BondKind::Up, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn down_double_single() {
        let bonds = vec![
            Bond::new(BondKind::Down, 0),
            Bond::new(BondKind::Double, 1),
            Bond::new(BondKind::Single, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn down_up() {
        let bonds = vec![Bond::new(BondKind::Down, 0), Bond::new(BondKind::Up, 1)];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn down_up_elided() {
        let bonds = vec![
            Bond::new(BondKind::Down, 0),
            Bond::new(BondKind::Up, 1),
            Bond::new(BondKind::Elided, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn down_up_double() {
        let bonds = vec![
            Bond::new(BondKind::Down, 0),
            Bond::new(BondKind::Up, 1),
            Bond::new(BondKind::Double, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Positive)))
    }

    #[test]
    fn down_up_double_elided() {
        let bonds = vec![
            Bond::new(BondKind::Down, 0),
            Bond::new(BondKind::Up, 1),
            Bond::new(BondKind::Double, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn down_up_single() {
        let bonds = vec![
            Bond::new(BondKind::Down, 0),
            Bond::new(BondKind::Up, 1),
            Bond::new(BondKind::Single, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn down_elided() {
        let bonds = vec![Bond::new(BondKind::Down, 0), Bond::new(BondKind::Elided, 1)];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn down_elided_double() {
        let bonds = vec![
            Bond::new(BondKind::Down, 0),
            Bond::new(BondKind::Elided, 1),
            Bond::new(BondKind::Double, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Positive)))
    }

    #[test]
    fn down_elided_double_elided() {
        let bonds = vec![
            Bond::new(BondKind::Down, 0),
            Bond::new(BondKind::Elided, 1),
            Bond::new(BondKind::Double, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn down_elided_single() {
        let bonds = vec![
            Bond::new(BondKind::Down, 0),
            Bond::new(BondKind::Elided, 1),
            Bond::new(BondKind::Single, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn down_single() {
        let bonds = vec![Bond::new(BondKind::Down, 0), Bond::new(BondKind::Single, 1)];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    // ---

    #[test]
    fn elided_up() {
        let bonds = vec![Bond::new(BondKind::Elided, 0), Bond::new(BondKind::Up, 1)];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn elided_up_double() {
        let bonds = vec![
            Bond::new(BondKind::Elided, 0),
            Bond::new(BondKind::Up, 1),
            Bond::new(BondKind::Double, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Positive)))
    }

    #[test]
    fn elided_up_double_elided() {
        let bonds = vec![
            Bond::new(BondKind::Elided, 0),
            Bond::new(BondKind::Up, 1),
            Bond::new(BondKind::Double, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn elided_up_single() {
        let bonds = vec![
            Bond::new(BondKind::Elided, 0),
            Bond::new(BondKind::Up, 1),
            Bond::new(BondKind::Single, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn elided_down() {
        let bonds = vec![Bond::new(BondKind::Elided, 0), Bond::new(BondKind::Down, 1)];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn elided_down_double() {
        let bonds = vec![
            Bond::new(BondKind::Elided, 0),
            Bond::new(BondKind::Down, 1),
            Bond::new(BondKind::Double, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Negative)))
    }

    #[test]
    fn elided_down_double_elided() {
        let bonds = vec![
            Bond::new(BondKind::Elided, 0),
            Bond::new(BondKind::Down, 1),
            Bond::new(BondKind::Double, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn elided_down_single() {
        let bonds = vec![
            Bond::new(BondKind::Elided, 0),
            Bond::new(BondKind::Down, 1),
            Bond::new(BondKind::Single, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn elided_double() {
        let bonds = vec![
            Bond::new(BondKind::Elided, 0),
            Bond::new(BondKind::Double, 1),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn elided_double_up() {
        let bonds = vec![
            Bond::new(BondKind::Elided, 0),
            Bond::new(BondKind::Double, 1),
            Bond::new(BondKind::Up, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Negative)))
    }

    #[test]
    fn elided_double_up_elided() {
        let bonds = vec![
            Bond::new(BondKind::Elided, 0),
            Bond::new(BondKind::Double, 1),
            Bond::new(BondKind::Up, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn elided_double_down() {
        let bonds = vec![
            Bond::new(BondKind::Elided, 0),
            Bond::new(BondKind::Double, 1),
            Bond::new(BondKind::Down, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(Some(Parity::Positive)))
    }

    #[test]
    fn elided_double_down_elided() {
        let bonds = vec![
            Bond::new(BondKind::Elided, 0),
            Bond::new(BondKind::Double, 1),
            Bond::new(BondKind::Down, 2),
            Bond::new(BondKind::Elided, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn elided_double_single() {
        let bonds = vec![
            Bond::new(BondKind::Elided, 0),
            Bond::new(BondKind::Double, 1),
            Bond::new(BondKind::Single, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn elided_single() {
        let bonds = vec![
            Bond::new(BondKind::Elided, 0),
            Bond::new(BondKind::Single, 1),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }

    #[test]
    fn elided_single_up() {
        let bonds = vec![
            Bond::new(BondKind::Elided, 0),
            Bond::new(BondKind::Single, 1),
            Bond::new(BondKind::Up, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn elided_single_down() {
        let bonds = vec![
            Bond::new(BondKind::Elided, 0),
            Bond::new(BondKind::Single, 1),
            Bond::new(BondKind::Down, 2),
        ];

        assert_eq!(trigonal_parity(&bonds), Err(()))
    }

    #[test]
    fn elided_single_single_down() {
        let bonds = vec![
            Bond::new(BondKind::Elided, 0),
            Bond::new(BondKind::Single, 1),
            Bond::new(BondKind::Single, 2),
            Bond::new(BondKind::Down, 3),
        ];

        assert_eq!(trigonal_parity(&bonds), Ok(None))
    }
}
