use crate::bignum::BigNum;
use crate::limbs::MASK_120;
use crate::cmp::{bignum_ct_ge, bignum_ct_eq, bignum_ct_cmp3};
use crate::addsub::{bignum_add_no_mod, bignum_sub_no_mod};

// helper
fn bn<let N: u32>(xs: [u128; N]) -> BigNum<N> { BigNum::from_limbs(xs) }

#[test]
fn cmp_eq_basic() {
    let a = bn::<3>([1u128, 2u128, 3u128]);
    let b = bn::<3>([1u128, 2u128, 3u128]);
    let c = bn::<3>([1u128, 2u128, 4u128]);

    assert(bignum_ct_eq(a, b) == 1u128);
    assert(bignum_ct_eq(a, c) == 0u128);
}

#[test]
fn cmp_ge_basic() {
    let x = bn::<2>([5u128, 0u128]);
    let y = bn::<2>([4u128, 0u128]);
    let z = bn::<2>([5u128, 1u128]);

    assert(bignum_ct_ge(x, y) == 1u128);
    assert(bignum_ct_ge(y, x) == 0u128);
    assert(bignum_ct_ge(z, x) == 1u128);
}

#[test]
fn cmp3_trichotomy() {
    let a = bn::<2>([7u128, 0u128]);
    let b = bn::<2>([7u128, 0u128]);
    let c = bn::<2>([6u128, 0u128]);
    let d = bn::<2>([7u128, 1u128]);

    assert(bignum_ct_cmp3(a, b) == 0);   // equal
    assert(bignum_ct_cmp3(c, a) == -1);  // less
    assert(bignum_ct_cmp3(d, a) == 1);   // greater
}

#[test]
fn add_no_mod_carry_propagates() {
    let a = bn::<3>([MASK_120, MASK_120, 0u128]);
    let b = bn::<3>([1u128, 0u128, 0u128]);
    let (s, carry) = bignum_add_no_mod(a, b);

    assert(s.limbs[0] == 0u128);
    assert(s.limbs[1] == 0u128);
    assert(s.limbs[2] == 1u128); 
    assert(carry == 0u128);
}

#[test]
fn sub_no_mod_borrow_propagates() {
    let a = bn::<3>([0u128, 0u128, 1u128]);
    let b = bn::<3>([1u128, 0u128, 0u128]);
    let (d, borrow) = bignum_sub_no_mod(a, b);

    assert(d.limbs[0] == MASK_120);
    assert(d.limbs[1] == MASK_120);
    assert(d.limbs[2] == 0u128);
    assert(borrow == 0u128);
}

#[test]
fn add_sub_inverse_identity() {
    let x = bn::<4>([11u128, 22u128, 33u128, 44u128]);
    let y = bn::<4>([5u128,  7u128,  0u128,  1u128]);

    let (xy, br1) = bignum_sub_no_mod(x, y);
    let (xy_plus_y, c1) = bignum_add_no_mod(xy, y);

    assert(xy_plus_y.limbs == x.limbs);
    assert(br1 <= 1u128);
    assert(c1  <= 1u128);
}

