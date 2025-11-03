use crate::bignum::BigNum;
use crate::limbs::{MASK_120, adc_120};
use crate::mul::{mul120_split, schoolbook_bn};

pub struct MontCtx<let N: u32, let P: u32> {
    pub p: BigNum<N>,   
    pub r2: BigNum<N>,  
    pub ninv: u128,    
}

pub fn mont_reduce<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>,
    t_in: [u128; P],
) -> (MontCtx<N, P>, BigNum<N>) {
    let mut t: [u128; P] = t_in;

    for i in 0..N {
        // m = (t[i] * ninv) mod B
        let m0: u128 = (t[i] & MASK_120) * (ctx.ninv & MASK_120);
        let m: u128 = m0 & MASK_120;

        // t[i..i+N+1] += m * p
        let mut carry: u128 = 0u128;
        for j in 0..N {
            let (lo, hi) = mul120_split(m, ctx.p.limbs[j]);

            let idx0: u32 = i + j;
            if idx0 < P {
                let (s0, c0) = adc_120(t[idx0], lo, carry);
                t[idx0] = s0;

                let idx1: u32 = idx0 + 1u32;
                if idx1 < P {
                    let (s1, c1) = adc_120(t[idx1], hi, c0);
                    t[idx1] = s1;
                    carry = c1;
                } else {
                    carry = 0u128;
                }
            }
        }

        // propagate carry
        let mut k: u32 = i + N + 1u32;
        for _it in 0..P {
            if carry == 0u128 {
                // nothing
            } else {
                if k < P {
                    let (sk, ck) = adc_120(t[k], 0u128, carry);
                    t[k] = sk;
                    carry = ck;
                    k += 1u32;
                }
            }
        }
    }

    let mut res_limbs: [u128; N] = [0u128; N];
    for i in 0..N {
        res_limbs[i] = t[i + N] & MASK_120;
    }
    let mut res = BigNum::<N>::from_limbs(res_limbs);

    let ge = crate::cmp::bignum_ct_ge(res, ctx.p);
    let (diff, _br) = crate::addsub::bignum_sub_no_mod(res, ctx.p);
    res = BigNum::<N>::ct_select(res, diff, ge);

    (ctx, res)
}

pub fn mont_mul<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>,
    a: BigNum<N>,
    b: BigNum<N>,
) -> (MontCtx<N, P>, BigNum<N>) {
    let t = schoolbook_bn::<N, P>(a, b);
    mont_reduce::<N, P>(ctx, t)
}

pub fn to_mont<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>,
    x: BigNum<N>,
) -> (MontCtx<N, P>, BigNum<N>) {
    mont_mul::<N, P>(ctx, x, ctx.r2)
}

pub fn from_mont<let N: u32, let P: u32>(
    ctx: MontCtx<N, P>,
    x: BigNum<N>,
) -> (MontCtx<N, P>, BigNum<N>) {
    // multiply by 1
    let one = BigNum::<N>::from_limbs({
        let mut arr: [u128; N] = [0u128; N];
        arr[0] = 1u128;
        arr
    });
    mont_mul::<N, P>(ctx, x, one)
}

