use crate::bignum::BigNum;
use crate::limbs::{MASK_120, adc_120};
use crate::mul::schoolbook_bn;
use crate::montgomery::{MontCtx, mont_mul, to_mont, from_mont};
use crate::params::{test_p_N2, test_r2_N2, test_ninv};
use crate::cmp::{bignum_ct_eq, bignum_ct_ge};
use crate::addsub::bignum_sub_no_mod;

fn bn2(a0: u128, a1: u128) -> BigNum<2> {
    BigNum::<2>::from_limbs([a0 & MASK_120, a1 & MASK_120])
}

fn naive_reduce_special_p(t: [u128; 4]) -> BigNum<2> {
    let t0 = t[0] & MASK_120;
    let t1 = t[1] & MASK_120;
    let t2 = t[2] & MASK_120;
    let t3 = t[3] & MASK_120;

    // fold high into low: (t0 + t2, t1 + t3)
    let (s0, _c0) = adc_120(t0, t2, 0u128);
    let (s1, _c1) = adc_120(t1, t3, 0u128);

    let mut r = BigNum::<2>::from_limbs([s0, s1]);

    let p = test_p_N2();
    let ge = bignum_ct_ge(r, p);
    let (diff, _br) = bignum_sub_no_mod(r, p);
    r = BigNum::<2>::ct_select(r, diff, ge);
    r
}

#[test]
fn mont_roundtrip_identity() {
    let ctx0 = MontCtx::<2, 4> { p: test_p_N2(), r2: test_r2_N2(), ninv: test_ninv() };

    let xs = [
        bn2(0u128, 0u128),
        bn2(1u128, 0u128),
        bn2(123u128, 0u128),
        bn2(MASK_120 - 5u128, 7u128),
    ];

    let mut ctx = ctx0;
    for k in 0..xs.len() {
        let x = xs[k];
        let (ctx2, xm) = to_mont::<2, 4>(ctx, x);
        let (ctx3, back) = from_mont::<2, 4>(ctx2, xm);
        ctx = ctx3; // carry the context forward
        assert(bignum_ct_eq(back, x) == 1u128);
    }
}

#[test]
fn mont_mul_matches_naive_for_special_p() {
    let ctx0 = MontCtx::<2, 4> { p: test_p_N2(), r2: test_r2_N2(), ninv: test_ninv() };
    let a = bn2(5u128, 7u128);
    let b = bn2(11u128, 3u128);

    // raw product and naive reduction
    let t = schoolbook_bn::<2, 4>(a, b);
    let naive = naive_reduce_special_p(t);

    // through Montgomery
    let (ctx1, am) = to_mont::<2, 4>(ctx0, a);
    let (ctx2, bm) = to_mont::<2, 4>(ctx1, b);
    let (ctx3, prod_m) = mont_mul::<2, 4>(ctx2, am, bm);
    let (_ctx4, prod) = from_mont::<2, 4>(ctx3, prod_m);

    assert(bignum_ct_eq(prod, naive) == 1u128);
}

#[test]
fn mont_one_behaves_as_one() {
    let ctx0 = MontCtx::<2, 4> { p: test_p_N2(), r2: test_r2_N2(), ninv: test_ninv() };

    let one = bn2(1u128, 0u128);
    let x   = bn2(12345u128, 678u128);

    let (ctx1, xm) = to_mont::<2, 4>(ctx0, x);
    let (ctx2, om) = to_mont::<2, 4>(ctx1, one);

    let (ctx3, y)  = mont_mul::<2, 4>(ctx2, xm, om);
    let (_ctx4, back) = from_mont::<2, 4>(ctx3, y);

    assert(bignum_ct_eq(back, x) == 1u128);
}

