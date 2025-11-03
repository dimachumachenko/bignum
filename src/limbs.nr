pub global LIMB_BITS: u32 = 120;
pub global MASK_120: u128 = (1u128 << (LIMB_BITS as u128)) - 1u128;

pub fn normalize_limb(x: u128) -> u128 {
	x& MASK_120
}

pub fn adc_120(a: u128, b: u128, cin: u128) -> (u128, u128) {
    let aa: u128 = a & MASK_120;
    let bb: u128 = b & MASK_120;
    let cc: u128 = cin & 1u128;

    let s: u128 = aa + bb + cc;

    let lo: u128 = s & MASK_120;

    let carry: u128 = s >> (LIMB_BITS as u128);

    (lo, carry)
}

pub fn sbb_120(a: u128, b: u128, bin: u128) -> (u128, u128) {
    let aa: u128 = a & MASK_120;
    let bb: u128 = (b & MASK_120) + (bin & 1u128);

    if aa >= bb {
        (aa - bb, 0u128)
    } else {
        ((aa + (1u128 << (LIMB_BITS as u128))) - bb, 1u128)
    }
}

pub fn add_no_mod<let N: u32>(a: [u128; N], b: [u128; N]) -> ([u128; N], u128) {
    let mut out: [u128; N] = [0u128; N];
    let mut carry: u128 = 0u128;
    for i in 0..N {
        let (lo, c) = adc_120(a[i], b[i], carry);
        out[i] = lo;
        carry = c;
    }
    (out, carry)
}

pub fn sub_no_mod<let N: u32>(a: [u128; N], b: [u128; N]) -> ([u128; N], u128) {
    let mut out: [u128; N] = [0u128; N];
    let mut borrow: u128 = 0u128;
    for i in 0..N {
        let (lo, br) = sbb_120(a[i], b[i], borrow);
        out[i] = lo;
        borrow = br;
    }
    (out, borrow)
}

pub fn cmp_limbs<let N: u32>(a: [u128; N], b: [u128; N]) -> i32 {
    let mut decided: bool = false;
    let mut res: i32 = 0;

    for i in 0..N {
        let j: u32 = (N - 1u32) - i;
        let ai: u128 = a[j] & MASK_120;
        let bi: u128 = b[j] & MASK_120;

        let gt_i: bool = ai > bi;
        let lt_i: bool = ai < bi;
        let neq_i: bool = gt_i | lt_i;             
        let pick:  bool = (!decided) & neq_i;

        let gt_val: i32 = if gt_i { 1 } else { 0 };
        let lt_val: i32 = if lt_i { 1 } else { 0 };
        let cand: i32  = gt_val - lt_val;        

        res     = if pick { cand } else { res };
        decided = decided | pick;              
    }
    res
}

pub fn normalize_limbs<let N: u32>(x: [u128; N]) -> [u128; N] {
    let mut out: [u128; N] = [0u128; N];
    for i in 0..N {
        out[i] = x[i] & MASK_120;
    }
    out
}

