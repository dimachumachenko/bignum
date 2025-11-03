use crate::limbs::MASK_120;
use crate::bignum::BigNum;
use crate::limbs::adc_120;

pub fn mul120_split(a: u128, b: u128) -> (u128, u128) {
    let a120: u128 = a & MASK_120;
    let b120: u128 = b & MASK_120;

    let mask60: u128 = (1u128 << 60u128) - 1u128;

    let a0: u128 = a120 & mask60;
    let a1: u128 = (a120 >> 60u128) & mask60;
    let b0: u128 = b120 & mask60;
    let b1: u128 = (b120 >> 60u128) & mask60;

    let p00: u128 = a0 * b0;
    let p01: u128 = a0 * b1;
    let p10: u128 = a1 * b0;
    let p11: u128 = a1 * b1;

    let mid: u128 = p01 + p10;
    let mid_low: u128 = mid & mask60;
    let mid_high: u128 = mid >> 60u128;

    let lo_raw: u128 = p00 + (mid_low << 60u128);
    let carry_lo: u128 = lo_raw >> 120u128;
    let lo: u128 = lo_raw & MASK_120;

    let hi: u128 = (p11 + mid_high + carry_lo) & MASK_120;
    (lo, hi)
}

pub fn schoolbook<let N: u32, let P: u32>(a: [u128; N], b: [u128; N]) -> [u128; P] {
    let mut out: [u128; P] = [0u128; P];

    for i in 0..N {
        for j in 0..N {
            let (lo, hi) = mul120_split(a[i], b[j]);

            let idx0: u32 = i + j;
            if idx0 < P {
                let (s0, c0) = adc_120(out[idx0], lo, 0u128);
                out[idx0] = s0;

                let idx1: u32 = idx0 + 1u32;
                let mut carry: u128 = c0;

                if idx1 < P {
                    let (s1, c1) = adc_120(out[idx1], hi, carry);
                    out[idx1] = s1;
                    carry = c1;
                } else {
                    carry = 0u128;
                }

                let mut k: u32 = idx1 + 1u32;
                for _t in 0..P {
                    if carry == 0u128 {
                    } else {
                        if k < P {
                            let (sk, ck) = adc_120(out[k], 0u128, carry);
                            out[k] = sk;
                            carry = ck;
                            k += 1u32;
                        }
                    }
                }
            }
        }
    }
    out
}

pub fn schoolbook_bn<let N: u32, let P: u32>(a: BigNum<N>, b: BigNum<N>) -> [u128; P] {
    schoolbook::<N, P>(a.limbs, b.limbs)
}

