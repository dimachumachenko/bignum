use crate::bignum::BigNum;
use crate::limbs::MASK_120;
use crate::montgomery::MontCtx;
use crate::params::{test_p_N2, test_r2_N2, test_ninv};
use crate::field::{
    field_from_big, field_to_big, field_add, field_sub,
    field_mul, field_square, field_neg
};
use crate::cmp::bignum_ct_eq;

fn bn2(a0: u128, a1: u128) -> BigNum<2> {
    BigNum::<2>::from_limbs([a0 & MASK_120, a1 & MASK_120])
}

#[test]
fn field_add_sub_roundtrip() {
    let ctx0 = MontCtx::<2, 4> { p: test_p_N2(), r2: test_r2_N2(), ninv: test_ninv() };
    let a_bn = bn2(123u128, 456u128);
    let b_bn = bn2(789u128, 10u128);

    let (ctx1, a) = field_from_big::<2, 4>(ctx0, a_bn);
    let (ctx2, b) = field_from_big::<2, 4>(ctx1, b_bn);

    let (ctx3, s) = field_add::<2, 4>(ctx2, a, b);
    let (ctx4, d) = field_sub::<2, 4>(ctx3, s, b);
    let (_ctx5, back) = field_to_big::<2, 4>(ctx4, d);

    assert(bignum_ct_eq(back, a_bn) == 1u128);
}

#[test]
fn field_mul_square_consistency_no_inverse() {
    let ctx0 = MontCtx::<2, 4> { p: test_p_N2(), r2: test_r2_N2(), ninv: test_ninv() };

    let (ctx1, a) = field_from_big::<2, 4>(ctx0, bn2(5u128, 7u128));
    let (ctx2, b) = field_from_big::<2, 4>(ctx1, bn2(11u128, 3u128));

    let (ctx3, a2)   = field_square::<2, 4>(ctx2, a);
    let (ctx4, a2b)  = field_mul::<2, 4>(ctx3, a2, b);

    let (ctx5, ab)   = field_mul::<2, 4>(ctx4, a, b);
    let (ctx6, aab)  = field_mul::<2, 4>(ctx5, a, ab);

    let (_ctx7, a2b_big) = field_to_big::<2, 4>(ctx6, a2b);
    let (_ctx8, aab_big) = field_to_big::<2, 4>(ctx6, aab);
    assert(bignum_ct_eq(a2b_big, aab_big) == 1u128);

    let (ctx9, sum)  = field_add::<2, 4>(ctx6, a, b);
    let (ctx10, sum2)= field_square::<2, 4>(ctx9, sum);

    let (ctx11, b2)  = field_square::<2, 4>(ctx10, b);
    let (ctx12, ab2) = field_add::<2, 4>(ctx11, ab, ab); // 2ab
    let (ctx13, rhs1)= field_add::<2, 4>(ctx12, a2, ab2);
    let (_ctx14, rhs) = field_add::<2, 4>(ctx13, rhs1, b2);

    let (_ctx15, lhs_big) = field_to_big::<2, 4>(ctx13, sum2);
    let (_ctx16, rhs_big) = field_to_big::<2, 4>(ctx13, rhs);
    assert(bignum_ct_eq(lhs_big, rhs_big) == 1u128);
}

#[test]
fn field_neg_property() {
    let ctx0 = MontCtx::<2, 4> { p: test_p_N2(), r2: test_r2_N2(), ninv: test_ninv() };
    let (ctx1, a) = field_from_big::<2, 4>(ctx0, bn2(42u128, 1u128));
    let (ctx2, neg_a) = field_neg::<2, 4>(ctx1, a);
    let (_ctx3, s) = field_add::<2, 4>(ctx2, a, neg_a);
    let (_ctx4, s_big) = field_to_big::<2, 4>(_ctx3, s);

    let zero = bn2(0u128, 0u128);
    assert(bignum_ct_eq(s_big, zero) == 1u128);
}

