use crate::limbs::MASK_120;
use crate::bignum::BigNum;
use crate::mul::{mul120_split, schoolbook, schoolbook_bn};

fn bn<let N: u32>(xs: [u128; N]) -> BigNum<N> { BigNum::from_limbs(xs) }

#[test]
fn mul120_basic_patterns() {
    let (lo, hi) = mul120_split(1u128, 1u128);
    assert(lo == 1u128);
    assert(hi == 0u128);

    let (lo2, hi2) = mul120_split(MASK_120, 1u128);
    assert(lo2 == MASK_120);
    assert(hi2 == 0u128);

    let two60: u128 = 1u128 << 60u128;
    let (lo3, hi3) = mul120_split(two60, two60);
    assert(lo3 == 0u128);
    assert(hi3 == 1u128);
}

#[test]
fn schoolbook_zero_and_one() {
    let a1: [u128; 1] = [0u128];
    let b1: [u128; 1] = [MASK_120];
    let p1 = schoolbook::<1, 2>(a1, b1);
    assert(p1[0] == 0u128);
    assert(p1[1] == 0u128);

    let a2: [u128; 1] = [1u128];
    let b2: [u128; 1] = [1u128];
    let p2 = schoolbook::<1, 2>(a2, b2);
    assert(p2[0] == 1u128);
    assert(p2[1] == 0u128);
}

#[test]
fn schoolbook_simple_carries() {
    let a: [u128; 1] = [MASK_120];
    let b: [u128; 1] = [2u128];
    let p = schoolbook::<1, 2>(a, b);
    assert(p[0] == (MASK_120 - 1u128)); // (2^120-2)
    assert(p[1] == 1u128);
}

#[test]
fn schoolbook_cross_limb_power() {
    let a = bn::<2>([0u128, 1u128]);
    let b = bn::<2>([0u128, 1u128]);
    let p = schoolbook_bn::<2, 4>(a, b);
    assert(p[0] == 0u128);
    assert(p[1] == 0u128);
    assert(p[2] == 1u128);
    assert(p[3] == 0u128);
}

#[test]
fn schoolbook_commutative_symmetry() {
    let a = bn::<2>([5u128, 7u128]);
    let b = bn::<2>([11u128, 3u128]);

    let p_ab = schoolbook_bn::<2, 4>(a, b);
    let p_ba = schoolbook_bn::<2, 4>(b, a);

    assert(p_ab[0] == p_ba[0]);
    assert(p_ab[1] == p_ba[1]);
    assert(p_ab[2] == p_ba[2]);
    assert(p_ab[3] == p_ba[3]);
}

#[test]
fn schoolbook_mixed_values_with_masks() {
    let a = bn::<3>([MASK_120, 0u128, 1u128]);
    let b = bn::<3>([2u128, MASK_120, 0u128]);
    let p = schoolbook_bn::<3, 6>(a, b);

    for i in 0..6 {
        assert((p[i] & (!MASK_120)) == 0u128);
    }

    let nz0: u128 = (p[0] != 0u128) as u128;
    let nz1: u128 = (p[1] != 0u128) as u128;
    let nz2: u128 = (p[2] != 0u128) as u128;
    let nz3: u128 = (p[3] != 0u128) as u128;
    let nz4: u128 = (p[4] != 0u128) as u128;
    let nz5: u128 = (p[5] != 0u128) as u128;
    let sum: u128 = nz0 + nz1 + nz2 + nz3 + nz4 + nz5;
    assert(sum >= 1u128);
}

