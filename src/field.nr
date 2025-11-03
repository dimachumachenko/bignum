use crate::bignum::BigNum;
use crate::limbs::MASK_120;
use crate::addsub::{bignum_add_no_mod, bignum_sub_no_mod};
use crate::cmp::{bignum_ct_ge, bignum_ct_eq};
use crate::montgomery::{MontCtx, to_mont, from_mont, mont_mul};

pub struct FieldElt<let N: u32, let P: u32> {
    pub mont: BigNum<N>,
}

pub fn field_from_big<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>, x: BigNum<N>
) -> (MontCtx<N, P>, FieldElt<N, P>) {
    let (ctx2, xm) = to_mont::<N, P>(ctx, x);
    (ctx2, FieldElt::<N, P> { mont: xm })
}

pub fn field_to_big<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>, a: FieldElt<N, P>
) -> (MontCtx<N, P>, BigNum<N>) {
    let (ctx2, x) = from_mont::<N, P>(ctx, a.mont);
    (ctx2, x)
}

pub fn field_add<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>, a: FieldElt<N, P>, b: FieldElt<N, P>
) -> (MontCtx<N, P>, FieldElt<N, P>) {
    let (sum_bn, _c) = bignum_add_no_mod::<N>(a.mont, b.mont);
    let ge = bignum_ct_ge(sum_bn, ctx.p);
    let (diff, _br) = bignum_sub_no_mod::<N>(sum_bn, ctx.p);
    let red = BigNum::<N>::ct_select(sum_bn, diff, ge);
    (ctx, FieldElt::<N, P> { mont: red })
}

pub fn field_sub<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>, a: FieldElt<N, P>, b: FieldElt<N, P>
) -> (MontCtx<N, P>, FieldElt<N, P>) {
    let (diff, borrow) = bignum_sub_no_mod::<N>(a.mont, b.mont);
    let (fixed, _c) = bignum_add_no_mod::<N>(
        diff,
        BigNum::<N>::ct_select(BigNum::<N>::from_limbs([0u128; N]), ctx.p, borrow),
    );
    (ctx, FieldElt::<N, P> { mont: fixed })
}

pub fn field_mul<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>, a: FieldElt<N, P>, b: FieldElt<N, P>
) -> (MontCtx<N, P>, FieldElt<N, P>) {
    let (ctx2, m) = mont_mul::<N, P>(ctx, a.mont, b.mont);
    (ctx2, FieldElt::<N, P> { mont: m })
}

pub fn field_square<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>, a: FieldElt<N, P>
) -> (MontCtx<N, P>, FieldElt<N, P>) {
    field_mul::<N, P>(ctx, a, a)
}

pub fn field_neg<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>, a: FieldElt<N, P>
) -> (MontCtx<N, P>, FieldElt<N, P>) {
    let is_zero = bignum_ct_eq(a.mont, BigNum::<N>::from_limbs([0u128; N]));
    let (tmp, _br) = bignum_sub_no_mod::<N>(ctx.p, a.mont);
    let neg = BigNum::<N>::ct_select(a.mont, tmp, ((1u128 - is_zero) & 1u128));
    (ctx, FieldElt::<N, P> { mont: neg })
}

pub fn field_pow_bnexp<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>, base: FieldElt<N, P>, exp: BigNum<N>
) -> (MontCtx<N, P>, FieldElt<N, P>) {
    // acc = 1 in Montgomery domain
    let one_bn = BigNum::<N>::from_limbs({
        let mut arr: [u128; N] = [0u128; N];
        arr[0] = 1u128;
        arr
    });
    let (mut ctx1, acc0) = to_mont::<N, P>(ctx, one_bn);
    let mut acc = FieldElt::<N, P> { mont: acc0 };
    let mut base_cur = base;

    let limb_bits: u32 = crate::limbs::LIMB_BITS;
    for i in 0..N {
        let ei = exp.limbs[i] & MASK_120;
        for bitpos in 0..limb_bits {
            let bit: u128 = (ei >> (bitpos as u128)) & 1u128;

            let (ctx2, prod) = field_mul::<N, P>(ctx1, acc, base_cur);
            let selected = BigNum::<N>::ct_select(acc.mont, prod.mont, bit);
            acc = FieldElt::<N, P> { mont: selected };

            let (ctx3, sq) = field_square::<N, P>(ctx2, base_cur);
            base_cur = sq;
            ctx1 = ctx3;
        }
    }
    (ctx1, acc)
}

pub fn field_inv_fermat<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>, a: FieldElt<N, P>, p_minus_2: BigNum<N>
) -> (MontCtx<N, P>, FieldElt<N, P>) {
    field_pow_bnexp::<N, P>(ctx, a, p_minus_2)
}

