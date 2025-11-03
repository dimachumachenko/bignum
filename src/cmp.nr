use crate::limbs::MASK_120;
use crate::bignum::BigNum;

pub fn bignum_ct_ge<let N: u32>(a: BigNum<N>, b: BigNum<N>) -> u128 {
    let mut gt: u128 = 0u128;
    let mut lt: u128 = 0u128;

    for i in 0..N {
        let j: u32 = (N - 1u32) - i;
        let ai: u128 = a.limbs[j] & MASK_120;
        let bi: u128 = b.limbs[j] & MASK_120;

        let ai_gt_bi: u128 = ((ai > bi) as u128) & 1u128;
        let ai_lt_bi: u128 = ((ai < bi) as u128) & 1u128;

        let undecided: u128 = (((gt | lt) == 0u128) as u128) & 1u128;
        gt |= ai_gt_bi & undecided;
        lt |= ai_lt_bi & undecided;
    }
    (((lt == 0u128) as u128) & 1u128)
}

pub fn bignum_ct_eq<let N: u32>(a: BigNum<N>, b: BigNum<N>) -> u128 {
    let mut acc: u128 = 0u128;
    for i in 0..N {
        let ai = a.limbs[i] & MASK_120;
        let bi = b.limbs[i] & MASK_120;
        acc |= ai ^ bi;
    }
    (((acc == 0u128) as u128) & 1u128)
}

pub fn bignum_ct_cmp3<let N: u32>(a: BigNum<N>, b: BigNum<N>) -> i32 {
    let mut gt: u128 = 0u128;
    let mut lt: u128 = 0u128;

    for i in 0..N {
        let j: u32 = (N - 1u32) - i;
        let ai: u128 = a.limbs[j] & MASK_120;
        let bi: u128 = b.limbs[j] & MASK_120;

        let ai_gt_bi: u128 = ((ai > bi) as u128) & 1u128;
        let ai_lt_bi: u128 = ((ai < bi) as u128) & 1u128;

        let undecided: u128 = (((gt | lt) == 0u128) as u128) & 1u128;
        gt |= ai_gt_bi & undecided;
        lt |= ai_lt_bi & undecided;
    }
    let pos: i32 = (gt as i32) & (! (lt as i32));
    let neg: i32 = (lt as i32) & (! (gt as i32));
    pos - neg
}

