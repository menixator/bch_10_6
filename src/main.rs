const MOD_N: i64 = 11;
const CORRECTION_MIN: usize = 1;
const CORRECTION_MAX: usize = 10;
const NUM_DIGITS: usize = 10;
const RADIX: u32 = 10;

#[derive(PartialEq, Debug)]
pub struct Correction {
    position: usize,
    magnitude: usize,
}

#[derive(PartialEq, Debug, Clone, Copy)]
pub struct PQR(i64, i64, i64);

#[derive(PartialEq, Debug, Clone, Copy)]
pub struct Syndromes(i64, i64, i64, i64);

#[derive(PartialEq, Debug)]
pub enum BCHError {
    InvalidFormat,
    SingleError {
        corrected: String,
        syn: Syndromes,
        pqr: PQR,
        error: Correction,
    },
    DoubleError {
        corrected: String,
        syn: Syndromes,
        pqr: PQR,
        errors: [Correction; 2],
    },
    TripleError {
        syn: Syndromes,
        pqr: PQR,
    },
}

type Input = [u8; NUM_DIGITS];

pub struct SQRTError;
fn sqrt(input: i64) -> Result<i64, SQRTError> {
    let input = input.rem_euclid(MOD_N);
    for i in 1..=MOD_N {
        if (i.pow(2)).rem_euclid(MOD_N) == input {
            return Ok(i);
        }
    }
    Err(SQRTError)
}

pub fn inverse(a: i64, n: i64) -> i64 {
    let mut t: i64 = 0;
    let mut newt: i64 = 1;
    let mut r: i64 = n;
    let mut newr: i64 = a;
    let mut q: i64;
    let mut temp: i64;

    while newr != 0 {
        q = r / newr; /* integer division */
        temp = newt; /* remember newt    */
        newt = t - q * newt;
        t = temp;
        temp = newr; /* remember newr    */
        newr = r - q * newr;
        r = temp;
    }
    if r > 1 {
        return -1;
    } /* not invertible */
    if t < 0 {
        t = t + n;
    } /* change to positive */
    return t;
}

pub fn calculate_syn(input: &Input, pass: u8) -> i64 {
    let mut sum = 0;
    for (index, num) in input.iter().enumerate() {
        sum += (*num as i64) * ((index + 1) as i64).pow(pass as u32)
    }
    sum
}

pub enum FixError {
    OutOfBounds,
    CorrectedValueTooBig,
}

pub fn bch_fix(input: Input, corrections: &[&Correction]) -> Result<String, FixError> {
    let mut fixbuffer: Input = [0; NUM_DIGITS];
    fixbuffer.copy_from_slice(&input[0..NUM_DIGITS]);

    for correction in corrections.iter() {
        if correction.position < CORRECTION_MIN || correction.position > CORRECTION_MAX {
            return Err(FixError::OutOfBounds);
        }

        let current_value: i8 = fixbuffer[correction.position - 1] as i8;
        let corrected: i8 = current_value - (correction.magnitude as i8);
        let new_value: u8 = (corrected as i64).rem_euclid(MOD_N) as u8;

        if (new_value as usize) >= CORRECTION_MAX {
            return Err(FixError::CorrectedValueTooBig);
        }

        fixbuffer[correction.position - 1] = new_value;
    }
    let corrected: String = fixbuffer
        .iter()
        .map(|digit| {
            let digit_as_char = std::char::from_digit(*digit as u32, RADIX);
            digit_as_char.unwrap()
        })
        .collect();

    Ok(corrected)
}
pub fn bch_decode(input: usize) -> Result<Syndromes, BCHError> {
    let input = {
        let input = format!("{:010}", input);
        if input.len() != NUM_DIGITS {
            return Err(BCHError::InvalidFormat);
        }

        let input: Vec<_> = input
            .chars()
            .map(|digit| {
                digit
                    .to_digit(RADIX)
                    .ok_or(BCHError::InvalidFormat)
                    .map(|val| val as u8)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut input_slice: Input = [0; NUM_DIGITS];
        input_slice.copy_from_slice(&input[0..NUM_DIGITS]);
        input_slice
    };
    let s1: i64 = calculate_syn(&input, 0).rem_euclid(MOD_N);
    let s2: i64 = calculate_syn(&input, 1).rem_euclid(MOD_N);
    let s3: i64 = calculate_syn(&input, 2).rem_euclid(MOD_N);
    let s4: i64 = calculate_syn(&input, 3).rem_euclid(MOD_N);

    if s1 == 0 && s2 == 0 && s3 == 0 && s4 == 0 {
        return Ok(Syndromes(s1, s2, s3, s4));
    }
    let p: i64 = (s2.pow(2) - (s1 * s3)).rem_euclid(MOD_N);
    let q: i64 = ((s1 * s4) - (s2 * s3)).rem_euclid(MOD_N);
    let r: i64 = (s3.pow(2) - (s2 * s4)).rem_euclid(MOD_N);
    let syn = Syndromes(s1, s2, s3, s4);
    let pqr = PQR(p, q, r);

    if p == 0 && q == 0 && r == 0 {
        let magnitude = s1;

        let pos = (s2 * inverse(s1.rem_euclid(MOD_N), MOD_N)).rem_euclid(MOD_N);
        let error = Correction {
            position: pos as usize,
            magnitude: magnitude as usize,
        };

        let corrected =
            bch_fix(input, &[&error]).map_err(|_| BCHError::TripleError { syn, pqr })?;
        // Check if correction will make the number greater than 10
        return Err(BCHError::SingleError {
            corrected,
            syn,
            pqr,
            error,
        });
    } else {
        let pol: i64 = sqrt((q.pow(2) - (4 * p * r)).rem_euclid(MOD_N))
            .map_err(|_| BCHError::TripleError { syn, pqr })?;
        // position
        let i = ((-q + pol) * (inverse(2 * p, MOD_N))).rem_euclid(MOD_N);
        let j = ((-q - pol) * (inverse(2 * p, MOD_N))).rem_euclid(MOD_N);

        // magnitudes
        let b = ((i * s1 - s2) * inverse((i - j).rem_euclid(MOD_N), MOD_N).rem_euclid(MOD_N))
            .rem_euclid(MOD_N);
        let a = (s1 - b).rem_euclid(MOD_N);

        let error1 = Correction {
            position: i as usize,
            magnitude: a as usize,
        };
        let error2 = Correction {
            position: j as usize,
            magnitude: b as usize,
        };

        let corrected =
            bch_fix(input, &[&error1, &error2]).map_err(|_| BCHError::TripleError { syn, pqr })?;

        return Err(BCHError::DoubleError {
            corrected,
            syn,
            pqr,
            errors: [error1, error2],
        });
    }
}
#[test]
fn test_no_error() {
    assert_eq!(bch_decode(3745195876), Ok(Syndromes(0, 0, 0, 0)));
}

#[test]
fn test_single_error() {
    assert_eq!(
        bch_decode(3945195876),
        Err(BCHError::SingleError {
            corrected: String::from("3745195876"),
            syn: Syndromes(2, 4, 8, 5),
            pqr: PQR(0, 0, 0),
            error: Correction {
                position: 2,
                magnitude: 2
            }
        })
    );
}

#[test]
fn test_double_err_1() {
    assert_eq!(
        bch_decode(3715195076),
        Err(
            //test
            //input(3715195076)
            //output(3745195876)
            BCHError::DoubleError {
                corrected: String::from("3745195876"),
                syn: Syndromes(0, 4, 0, 3),
                pqr: PQR(5, 0, 10),
                errors: [
                    //i,a
                    Correction {
                        position: 8,
                        magnitude: 3
                    },
                    //j,b
                    Correction {
                        position: 3,
                        magnitude: 8
                    }
                ]
            }
        )
    );
}

#[test]
fn test_double_err_2() {
    //test
    //input(0743195876)
    //output(3745195876)
    assert_eq!(
        bch_decode(0743195876),
        Err(BCHError::DoubleError {
            corrected: String::from("3745195876"),
            syn: Syndromes(6, 0, 9, 1),
            pqr: PQR(1, 6, 4),
            errors: [
                //i,a
                Correction {
                    position: 4,
                    magnitude: 9
                },
                //j,b
                Correction {
                    position: 1,
                    magnitude: 8
                }
            ]
        })
    )
}

#[test]
fn test_double_err_3() {
    //test
    //input(3745195840)
    //output(3745195876)
    assert_eq!(
        bch_decode(3745195840),
        Err(BCHError::DoubleError {
            corrected: String::from("3745195876"),
            syn: Syndromes(2, 1, 4, 8),
            pqr: PQR(4, 1, 8),
            errors: [
                //i,a
                Correction {
                    position: 10,
                    magnitude: 5
                },
                //j,b
                Correction {
                    position: 9,
                    magnitude: 8
                }
            ]
        })
    );
}

#[test]
fn test_double_err_4() {
    //test
    //input(8745105876)
    //output(3745195876)
    assert_eq!(
        bch_decode(8745105876),
        Err(BCHError::DoubleError {
            corrected: String::from("3745195876"),
            syn: Syndromes(7, 6, 0, 8),
            pqr: PQR(3, 1, 7),
            errors: [
                //i,a
                Correction {
                    position: 6,
                    magnitude: 2
                },
                //j,b
                Correction {
                    position: 1,
                    magnitude: 5
                }
            ]
        })
    );
}

#[test]
fn test_double_err_5() {
    //test
    //input(3745102876)
    //output(3745195876)
    assert_eq!(
        bch_decode(3745102876),
        Err(BCHError::DoubleError {
            corrected: String::from("3745195876"),
            syn: Syndromes(10, 2, 2, 8),
            pqr: PQR(6, 10, 10),
            errors: [
                //i,a
                Correction {
                    position: 6,
                    magnitude: 2
                },
                //j,b
                Correction {
                    position: 7,
                    magnitude: 8
                }
            ]
        })
    );
}

#[test]
fn test_double_err_6() {
    //test
    //input(1145195876)
    //output(3745195876)
    assert_eq!(
        bch_decode(1145195876),
        Err(BCHError::DoubleError {
            corrected: String::from("3745195876"),
            syn: Syndromes(3, 8, 7, 5),
            pqr: PQR(10, 3, 9),
            errors: [
                //i,a
                Correction {
                    position: 1,
                    magnitude: 9
                },
                //j,b
                Correction {
                    position: 2,
                    magnitude: 5
                }
            ]
        })
    );
}

#[test]
fn test_double_err_7() {
    //test
    //input(3745191976)
    //output(3745195876)
    assert_eq!(
        bch_decode(3745191976),
        Err(BCHError::DoubleError {
            corrected: String::from("3745195876"),
            syn: Syndromes(8, 2, 0, 9),
            pqr: PQR(4, 6, 4),
            errors: [
                //i,a
                Correction {
                    position: 8,
                    magnitude: 1
                },
                //j,b
                Correction {
                    position: 7,
                    magnitude: 7
                }
            ]
        })
    );
}

#[test]
fn test_double_err_8() {
    //test
    //input(3745190872)
    //output(3745195876)
    assert_eq!(
        bch_decode(3745190872),
        Err(BCHError::DoubleError {
            corrected: String::from("3745195876"),
            syn: Syndromes(2, 2, 4, 5),
            pqr: PQR(7, 2, 6),
            errors: [
                // i,a
                Correction {
                    position: 7,
                    magnitude: 6
                },
                // j,b
                Correction {
                    position: 10,
                    magnitude: 7
                },
            ]
        })
    );
}

#[test]
fn test_triple_err_1() {
    assert_eq!(
        bch_decode(2745795878),
        Err(BCHError::TripleError {
            syn: Syndromes(7, 5, 8, 10),
            pqr: PQR(2, 8, 3)
        })
    )
}

#[test]
fn test_triple_err_2() {
    assert_eq!(
        bch_decode(3742102896),
        Err(BCHError::TripleError {
            syn: Syndromes(9, 8, 6, 9),
            pqr: PQR(10, 0, 8)
        })
    )
}

#[test]
fn test_triple_err_3() {
    assert_eq!(
        bch_decode(1115195876),
        Err(BCHError::TripleError {
            syn: Syndromes(0, 10, 2, 1),
            pqr: PQR(1, 2, 5)
        })
    )
}

#[test]
fn test_triple_err_4() {
    assert_eq!(
        bch_decode(3121195876),
        Err(BCHError::TripleError {
            syn: Syndromes(10, 10, 4, 5),
            pqr: PQR(5, 10, 10)
        })
    )
}

fn main() {
    for i in 0..(10_usize.pow(10)) {
        println!("{} -> {:?}", i, bch_decode(i));
    }
}
